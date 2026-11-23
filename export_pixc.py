##################################################
## Authors: Merritt Harlan (USGS), Alexander Hernes (USGS)
## Last Updated: November 2026
#-------------------------------------------------
# Usage (from the repo root): python src/export_pixc.py [--config config.yaml] [--lakes NAME ...]
#                              [--time-limit-min 100] [--max-granules N]
#-------------------------------------------------
# DESCRIPTION: Finds SWOT PIXC granules over the monitored lakes that have not been
# processed yet, downloads each granule once (even when it covers several lakes),
# and saves the pixels assigned to each lake as data/pixc/<lake>/<granule>.parquet.
# Every (lake, granule) is logged in data/processed_granules.csv so later runs skip
# it, including granules that had no pixels. Needs Earthdata credentials: a .netrc
# locally, or EARTHDATA_USERNAME / EARTHDATA_PASSWORD in CI.
#-------------------------------------------------
# DISCLAIMER: This software is preliminary or provisional and is subject to revision.
# It is being provided to meet the need for timely best science. The software has not
# received final approval by the U.S. Geological Survey (USGS). No warranty, expressed
# or implied, is made by the USGS or the U.S. Government as to the functionality of the
# software and related material nor shall the fact of release constitute any such
# warranty. The software is provided on the condition that neither the USGS nor the U.S.
# Government shall be held liable for any damages resulting from the authorized or
# unauthorized use of the software. This software is provided 'AS IS'
##################################################

import argparse
import os
import shutil
import time
import warnings
from datetime import datetime, timezone
from itertools import islice
from pathlib import Path

import earthaccess
import h5py
import numpy as np
import pandas as pd
import scipy.ndimage
import shapely
from pyproj import Transformer

# earthaccess calls its own soon-to-change DataGranule.size(); the warning is not ours to fix
warnings.filterwarnings("ignore", category=FutureWarning, module="earthaccess")
from common import (DONE_STATUSES, lake_pixel_dir, load_config, load_lakes, manifest_record,
                    parse_granule_name, pass_matches, qc_mask, read_manifest, write_manifest)

# PIXC variables the assignment, QC and aggregation steps need
REQUIRED_VARS = [
    "height", "geoid", "solid_earth_tide", "load_tide_fes", "pole_tide",
    "classification", "geolocation_qual", "sig0_qual", "interferogram_qual",
    "classification_qual", "cross_track", "phase_noise_std", "dheight_dphase",
    "dlatitude_dphase", "dlongitude_dphase", "eff_num_rare_looks", "eff_num_medium_looks",
    "interferogram", "power_plus_y", "power_minus_y",
]


def batch_iter(iterable, size):
    it = iter(iterable)
    while batch := list(islice(it, size)):
        yield batch


# -------------------------------------------------------------------------
# Pixel assignment
#   Polygon method: pixel's PIXC geolocation falls inside the lake polygon.
#   Refined method (adapted from the SWOT Pixel Assignment Example notebook, B. Williams, JPL):
#     1) segment detected water (classification > 2) over the full tile in the slant plane
#     2) find the dominant segment among open-water pixels (class 4) inside the polygon
#     3) refine geolocation by moving each pixel along its height-sensitivity vector to
#        that segment's median height
#     4) keep detected water (class > 2, not dark water) inside a dilated polygon, in the
#        dominant segment, without height outliers
# -------------------------------------------------------------------------

def dilate_polygon(geom, buffer_m):
    """Buffer a lat/lon polygon by buffer_m metres in a local azimuthal equidistant
    projection, so the buffer is the same distance in every direction."""
    c = geom.centroid
    aeqd = f"+proj=aeqd +lat_0={c.y} +lon_0={c.x} +datum=WGS84 +units=m"
    fwd = Transformer.from_crs("EPSG:4326", aeqd, always_xy=True)
    inv = Transformer.from_crs(aeqd, "EPSG:4326", always_xy=True)

    def apply(t):
        return lambda xy: np.column_stack(t.transform(xy[:, 0], xy[:, 1]))
    return shapely.transform(shapely.transform(geom, apply(fwd)).buffer(buffer_m), apply(inv))


def segment_slant_plane(classification, azimuth_index, range_index):
    """2D label image of connected detected water (classification > 2); 0 = not water."""
    mask = np.zeros((azimuth_index.max() + 1, range_index.max() + 1), dtype=bool)
    water = classification > 2
    mask[azimuth_index[water], range_index[water]] = True
    labels, _ = scipy.ndimage.label(mask)
    return labels


def slant_window_mask(lat, lon, azimuth_index, range_index, bounds, slant_buff=10):
    """Pixels inside the slant-plane rectangle covering every pixel geolocated within
    bounds = (minx, miny, maxx, maxy), padded by slant_buff. None if nothing is inside."""
    minx, miny, maxx, maxy = bounds
    in_box = (lon >= minx) & (lon <= maxx) & (lat >= miny) & (lat <= maxy)
    if not np.any(in_box):
        return None
    az0, az1 = azimuth_index[in_box].min() - slant_buff, azimuth_index[in_box].max() + slant_buff
    r0, r1 = range_index[in_box].min() - slant_buff, range_index[in_box].max() + slant_buff
    return ((azimuth_index >= az0) & (azimuth_index <= az1) &
            (range_index >= r0) & (range_index <= r1))


def _dominant_label(labels):
    labels = labels[labels > 0]
    if labels.size == 0:
        return None
    vals, counts = np.unique(labels, return_counts=True)
    return int(vals[np.argmax(counts)])


def contains(geom, lon, lat):
    """Vectorized point-in-polygon that ignores non-finite coordinates."""
    out = np.zeros(lat.shape, dtype=bool)
    ok = np.isfinite(lat) & np.isfinite(lon)
    if np.any(ok):
        out[ok] = shapely.contains_xy(geom, lon[ok], lat[ok])
    return out


def refined_pixel_assignment(lat, lon, height, wse, classification, seg_label,
                             dlat_dphase, dlon_dphase, dheight_dphase,
                             lake_geom, lake_geom_dilated, max_dh=50.0, iqr_k=3.0,
                             median_wse_bounds=None):
    """Returns dict: assigned (bool array), latitude_refined, longitude_refined, status."""
    n = lat.shape[0]
    out = {"assigned": np.zeros(n, dtype=bool),
           "latitude_refined": np.full(n, np.nan),
           "longitude_refined": np.full(n, np.nan),
           "status": "ok"}

    # Step 2: open water in the polygon -> dominant segment and its median height
    trusted = (classification == 4) & contains(lake_geom, lon, lat) & np.isfinite(height)
    if median_wse_bounds is not None:
        trusted &= (wse >= median_wse_bounds[0]) & (wse <= median_wse_bounds[1])
    dom = _dominant_label(seg_label[trusted])
    if dom is None:
        out["status"] = "no_open_water_in_polygon"
        return out
    med = np.nanmedian(height[trusted & (seg_label == dom)])

    # Step 3: refined geolocation
    dh = np.clip(height - med, -max_dh, max_dh)
    with np.errstate(divide="ignore", invalid="ignore"):
        lat_ref = lat - (dlat_dphase / dheight_dphase) * dh
        lon_ref = lon - (dlon_dphase / dheight_dphase) * dh
    out["latitude_refined"], out["longitude_refined"] = lat_ref, lon_ref

    # Step 4: detected water in the dilated polygon, dominant segment, no outliers
    det = (classification > 2) & (classification != 5)
    water_in = det & contains(lake_geom_dilated, lon_ref, lat_ref)
    dom2 = _dominant_label(seg_label[water_in])
    if dom2 is None:
        out["status"] = "no_water_in_dilated_polygon"
        return out
    water_seg_in = water_in & (seg_label == dom2)
    p25, med2, p75 = np.nanpercentile(height[water_seg_in], [25, 50, 75])
    iqr = p75 - p25
    outliers = np.abs(height - med2) > iqr_k * iqr if iqr > 0 else np.zeros(n, dtype=bool)
    out["assigned"] = water_seg_in & ~outliers
    return out


# -------------------------------------------------------------------------
# Reading a granule
# -------------------------------------------------------------------------

def read_masked(ds, mask, lo, hi):
    """Read one contiguous slice then mask in memory (much faster than h5py boolean
    indexing); float fill values become NaN."""
    a = ds[lo:hi][mask[lo:hi]]
    fill = ds.attrs.get("_FillValue")
    if fill is not None and a.dtype.kind == "f":
        a = a.copy()
        a[a == np.asarray(fill).ravel()[0]] = np.nan
    return a


def extract_lake(grp, full, lake, cfg):
    """Pixels from an open PIXC group for one lake. Returns (DataFrame or None, refined_status)."""
    ra, pf = cfg["refined_assignment"], cfg["pixel_filter"]
    lat, lon = full["latitude"], full["longitude"]
    geom = lake["geometry"]
    sw_lon, sw_lat, ne_lon, ne_lat = lake["sw_lon"], lake["sw_lat"], lake["ne_lon"], lake["ne_lat"]
    in_bbox = (lat >= sw_lat) & (lat <= ne_lat) & (lon >= sw_lon) & (lon <= ne_lon)

    if ra["enabled"]:
        minx, miny, maxx, maxy = lake["geometry_dilated"].bounds
        b = ra["bbox_buff_deg"]
        bounds = (min(minx, sw_lon) - b, min(miny, sw_lat) - b,
                  max(maxx, ne_lon) + b, max(maxy, ne_lat) + b)
        window = slant_window_mask(lat, lon, full["azimuth_index"], full["range_index"],
                                   bounds, slant_buff=ra["slant_buff_pix"])
        if window is None:
            return None, "no_pixels_in_window"
    else:
        window = in_bbox
    if not window.any():
        return None, "no_pixels_in_window"

    idx = np.flatnonzero(window)
    lo, hi = idx[0], idx[-1] + 1
    d = {v: read_masked(grp[v], window, lo, hi)
         for v in REQUIRED_VARS + list(cfg.get("extra_pixc_variables") or [])}
    ifg = d.pop("interferogram")
    d["interferogram_real"], d["interferogram_imag"] = ifg[:, 0], ifg[:, 1]
    d["latitude"], d["longitude"] = lat[window], lon[window]
    d["wse_m_egm2008"] = d["height"] - d["geoid"] - d["solid_earth_tide"] - d["load_tide_fes"] - d["pole_tide"]
    d["wse_m_navd88"] = d["wse_m_egm2008"] + lake["NAVD88-EGM2008_m"]
    d["in_polygon"] = in_bbox[window] & contains(geom, d["longitude"], d["latitude"])

    n = len(d["latitude"])
    if ra["enabled"]:
        seg_label = full["segmentation"]()[full["azimuth_index"][window], full["range_index"][window]]
        bounds = (lake["elev_min"], lake["elev_max"]) if ra["use_elev_bounds_for_median"] else None
        r = refined_pixel_assignment(
            lat=d["latitude"], lon=d["longitude"], height=d["height"], wse=d["wse_m_navd88"],
            classification=d["classification"], seg_label=seg_label,
            dlat_dphase=d["dlatitude_dphase"], dlon_dphase=d["dlongitude_dphase"],
            dheight_dphase=d["dheight_dphase"], lake_geom=geom,
            lake_geom_dilated=lake["geometry_dilated"], max_dh=ra["max_dh_m"],
            iqr_k=ra["iqr_k"], median_wse_bounds=bounds)
        d["refined_assigned"] = r["assigned"]
        d["latitude_refined"], d["longitude_refined"] = r["latitude_refined"], r["longitude_refined"]
        d["segmentation_label"] = seg_label.astype(np.int32)
        refined_status = r["status"]
    else:
        d["refined_assigned"] = np.zeros(n, dtype=bool)
        d["latitude_refined"] = d["longitude_refined"] = np.full(n, np.nan)
        d["segmentation_label"] = np.full(n, -1, dtype=np.int32)
        refined_status = "not_run"

    ct = np.abs(d["cross_track"])
    keep = (np.isin(d["classification"], pf["classification_vals"]) &
            (ct >= pf["cross_track_min_m"]) & (ct <= pf["cross_track_max_m"]) &
            (d["in_polygon"] | d["refined_assigned"]))
    if not keep.any():
        return None, refined_status

    df = pd.DataFrame(d).loc[keep].reset_index(drop=True)
    qc = qc_mask(df, lake, cfg)
    df["in_polygon_qc"] = df["in_polygon"] & qc
    df["refined_qc"] = df["refined_assigned"] & qc
    return df, refined_status


def process_granule(path, granule, lakes, cfg):
    """Open a granule once and extract every lake it covers.
    Returns a list of (lake, DataFrame or None, refined_status, error message or '')."""
    info = parse_granule_name(granule)
    results = []
    with h5py.File(path, "r") as f:
        grp = f["pixel_cloud"]
        full = {"latitude": grp["latitude"][:], "longitude": grp["longitude"][:]}
        if cfg["refined_assignment"]["enabled"]:
            full["azimuth_index"] = grp["azimuth_index"][:].astype(np.int64)
            full["range_index"] = grp["range_index"][:].astype(np.int64)
            cache = {}

            def segmentation():  # whole-tile segmentation, computed once and only if needed
                if "labels" not in cache:
                    cache["labels"] = segment_slant_plane(
                        grp["classification"][:], full["azimuth_index"], full["range_index"])
                return cache["labels"]
            full["segmentation"] = segmentation
        looks_to_efflooks = float(np.squeeze(grp.attrs.get("looks_to_efflooks", np.nan)))

        for lake in lakes:
            try:
                df, status = extract_lake(grp, full, lake, cfg)
            except Exception as e:  # keep going with the other lakes
                results.append((lake, None, "", f"{type(e).__name__}: {e}"))
                continue
            if df is not None:
                meta = {"lake_name": str(lake["name"]), "granule": granule, **info,
                        "refined_status": status, "looks_to_efflooks": looks_to_efflooks}
                for i, (k, v) in enumerate(meta.items()):
                    df.insert(i, k, v)
            results.append((lake, df, status, ""))
    return results


# -------------------------------------------------------------------------
# Search and run
# -------------------------------------------------------------------------

def granule_name(g):
    for link in g.data_links():
        base = link.rsplit("/", 1)[-1]
        if base.endswith(".nc"):
            return base[:-3]
    return g["umm"]["GranuleUR"]


def download(granules, tmp_dir):
    """Download granules; if the batch fails (e.g. a server error on one file), retry them
    one at a time so a single bad file doesn't stop the run. Returns the local paths that
    succeeded; missing ones are logged as errors by the caller and retried next run."""
    try:
        return earthaccess.download(granules, local_path=str(tmp_dir), threads=len(granules)) or []
    except Exception as e:
        if len(granules) == 1:
            print(f"  Download failed: {granule_name(granules[0])}: {e}")
            return []
        print(f"  Batch download failed ({type(e).__name__}); retrying one at a time")
    files = []
    for g in granules:
        try:
            files += earthaccess.download([g], local_path=str(tmp_dir), threads=1) or []
        except Exception as e:
            print(f"  Download failed: {granule_name(g)}: {e}")
    return files


def find_new_granules(lakes, cfg, done):
    """{granule_name: {'granule': DataGranule, 'lakes': [...]}} for unprocessed granules
    in season and on each lake's passes."""
    s = cfg["search"]
    end = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S")
    todo = {}
    for lake in lakes:
        results = earthaccess.search_data(
            short_name=s["collection"], temporal=(s["start_date"], end),
            bounding_box=(lake["sw_lon"], lake["sw_lat"], lake["ne_lon"], lake["ne_lat"]))
        skipped = {"unparsed name": 0, "older version": 0, "off season": 0,
                   "other pass": 0, "done": 0}
        example = None
        # NASA sometimes reprocesses a granule (..._PID0_01 -> ..._PID0_02); only the
        # newest version of each cycle/pass/tile is downloaded. Names sort by CRID, counter.
        newest = {}
        for g in results:
            n = granule_name(g)
            info = parse_granule_name(n)
            if info is not None:
                key = (info["cycle_number"], info["pass_number"], info["tile_name"])
                newest[key] = max(newest.get(key, n), n)
        for g in results:
            name = granule_name(g)
            example = example or name
            info = parse_granule_name(name)
            if info is None:
                skipped["unparsed name"] += 1
            elif newest[(info["cycle_number"], info["pass_number"], info["tile_name"])] != name:
                skipped["older version"] += 1
            elif info["datetime_utc"].month not in s["season_months"]:
                skipped["off season"] += 1
            elif not pass_matches(lake["passes"], name, info["pass_number"]):
                skipped["other pass"] += 1
            elif (str(lake["name"]), name) in done:
                skipped["done"] += 1
            else:
                todo.setdefault(name, {"granule": g, "lakes": [], "time": info["datetime_utc"]})
                todo[name]["lakes"].append(lake)
        n_new = len(results) - sum(skipped.values())
        why = ", ".join(f"{v} {k}" for k, v in skipped.items() if v)
        print(f"{lake['name']}: {len(results)} granules found, {n_new} new"
              + (f" (skipped: {why})" if why else ""))
        if results and n_new == 0 and skipped["done"] == 0:
            print(f"    passes in GeoJSON: {lake['passes']}; example granule: {example}")
    return todo


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--config", default="config.yaml")
    ap.add_argument("--lakes", nargs="*", help="Only process these lake names")
    ap.add_argument("--time-limit-min", type=float, default=None,
                    help="Stop starting new downloads after this many minutes")
    ap.add_argument("--max-granules", type=int, default=None)
    args = ap.parse_args()
    t0 = time.time()

    cfg = load_config(args.config)
    data_root, tmp_dir = Path(cfg["data_root"]), Path(cfg["tmp_dir"])
    lakes = load_lakes(cfg["lakes_geojson"])
    if args.lakes:
        lakes = [l for l in lakes if l["name"] in args.lakes]
    for lake in lakes:
        lake["geometry_dilated"] = dilate_polygon(lake["geometry"], cfg["refined_assignment"]["buffer_m"])

    manifest = read_manifest(data_root)
    done = set(zip(manifest.loc[manifest["status"].isin(DONE_STATUSES), "lake_name"],
                   manifest.loc[manifest["status"].isin(DONE_STATUSES), "granule"]))

    auth = earthaccess.login(strategy="environment" if os.getenv("CI") else "all")
    if not auth.authenticated:
        raise SystemExit("Earthdata login failed. Set EARTHDATA_USERNAME / EARTHDATA_PASSWORD "
                         "or a .netrc (account: https://urs.earthdata.nasa.gov/home).")

    todo = find_new_granules(lakes, cfg, done)
    names = sorted(todo, key=lambda n: todo[n]["time"])  # oldest first
    if args.max_granules is not None:
        names = names[:args.max_granules]
    print(f"\n{len(names)} granules to download")

    version = str(cfg.get("processing_version", ""))
    n_done = 0
    try:
        for batch in batch_iter(names, cfg["search"]["download_batch_size"]):
            if args.time_limit_min and (time.time() - t0) / 60 > args.time_limit_min:
                print(f"Time limit reached; {len(names) - n_done} granules left for the next run")
                break
            tmp_dir.mkdir(parents=True, exist_ok=True)
            files = download([todo[n]["granule"] for n in batch], tmp_dir)
            by_name = {Path(p).stem: Path(p) for p in files or []}
            records = []
            for name in batch:
                path = by_name.get(name)
                if path is None or not path.exists():
                    records += [manifest_record(l["name"], name, "error", message="download failed",
                                                processing_version=version) for l in todo[name]["lakes"]]
                    continue
                try:
                    results = process_granule(path, name, todo[name]["lakes"], cfg)
                except Exception as e:
                    results = [(l, None, "", f"{type(e).__name__}: {e}") for l in todo[name]["lakes"]]
                for lake, df, refined_status, err in results:
                    if err:
                        print(f"  ERROR {lake['name']} {name}: {err}")
                        records.append(manifest_record(lake["name"], name, "error", message=err,
                                                       processing_version=version))
                    elif df is None:
                        records.append(manifest_record(lake["name"], name, "no_pixels",
                                                       refined_status=refined_status,
                                                       processing_version=version))
                    else:
                        out = lake_pixel_dir(data_root, lake["name"]) / f"{name}.parquet"
                        out.parent.mkdir(parents=True, exist_ok=True)
                        df.to_parquet(out, index=False, compression="zstd")
                        records.append(manifest_record(lake["name"], name, "ok", n_pixels=len(df),
                                                       refined_status=refined_status,
                                                       processing_version=version))
                        print(f"  {lake['name']} {name}: {len(df):,} pixels")
                path.unlink(missing_ok=True)
                n_done += 1
            # Save progress after every batch so a timeout loses at most one batch
            manifest = write_manifest(pd.concat([manifest, pd.DataFrame(records)], ignore_index=True),
                                      data_root)
    finally:
        shutil.rmtree(tmp_dir, ignore_errors=True)

    print(f"\nProcessed {n_done} granules in {(time.time() - t0) / 60:.1f} min")


if __name__ == "__main__":
    main()
