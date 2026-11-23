"""Shared helpers for the SWOT lake pipeline: config, lakes, granule names, manifest, QC."""

import json
import re
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from shapely.geometry import shape

GRANULE_RE = re.compile(
    r"SWOT_L2_HR_PIXC_(?P<cycle>\d{3})_(?P<pass>\d{3})_(?P<tile>\d{3}[LR])_"
    r"(?P<start>\d{8}T\d{6})_(?P<end>\d{8}T\d{6})_(?P<crid>\w{4})_(?P<counter>\d{2})"
)

MANIFEST_NAME = "processed_granules.csv"
MANIFEST_COLS = ["lake_name", "granule", "datetime_utc", "status", "n_pixels",
                 "refined_status", "processing_version", "processed_utc", "message"]
DONE_STATUSES = {"ok", "no_pixels"}  # anything else (errors) is retried on the next run


def load_config(path="config.yaml"):
    with open(path) as f:
        return yaml.safe_load(f)


def parse_granule_name(name):
    """Cycle, pass, tile, start time and CRID from a PIXC granule name, or None."""
    m = GRANULE_RE.search(name)
    if not m:
        return None
    d = m.groupdict()
    return {
        "cycle_number": int(d["cycle"]),
        "pass_number": int(d["pass"]),
        "tile_name": d["tile"],
        "datetime_utc": pd.to_datetime(d["start"], format="%Y%m%dT%H%M%S"),
        "crid": d["crid"],
    }


def load_lakes(geojson_path):
    """Lakes as a list of dicts: the GeoJSON properties plus 'geometry' (shapely, EPSG:4326),
    'passes' (list of strings, empty = any pass) and 'properties' (the untouched originals)."""
    with open(geojson_path) as f:
        gj = json.load(f)
    lakes = []
    for feat in gj["features"]:
        props = feat["properties"]
        lake = dict(props)
        lake["properties"] = props
        lake["geometry"] = shape(feat["geometry"])
        passes = props.get("passes") or []
        if isinstance(passes, (str, int)):
            passes = [passes]
        lake["passes"] = [str(p).strip() for p in passes if str(p).strip()]
        lakes.append(lake)
    return lakes


def pass_matches(lake_passes, granule, pass_number):
    """True if the granule is on one of the lake's passes. Plain numbers ("14", "014", 14)
    are compared to the pass number; anything else (e.g. "014_123L") is matched as text
    in the granule name, like the original script. No passes listed = any pass."""
    if not lake_passes:
        return True
    for p in lake_passes:
        if p.isdigit():
            if int(p) == pass_number:
                return True
        elif p in granule:
            return True
    return False


def safe_name(name):
    return re.sub(r"[^A-Za-z0-9._-]+", "_", str(name)).strip("_")


def lake_pixel_dir(data_root, lake_name):
    return Path(data_root) / "pixc" / safe_name(lake_name)


def read_manifest(data_root):
    path = Path(data_root) / MANIFEST_NAME
    if not path.exists():
        return pd.DataFrame(columns=MANIFEST_COLS)
    return pd.read_csv(path, dtype={"granule": str, "lake_name": str})


def write_manifest(manifest, data_root):
    path = Path(data_root) / MANIFEST_NAME
    path.parent.mkdir(parents=True, exist_ok=True)
    manifest = (manifest[MANIFEST_COLS]
                .drop_duplicates(["lake_name", "granule"], keep="last")
                .sort_values(["lake_name", "datetime_utc", "granule"]))
    manifest.to_csv(path, index=False)
    return manifest


def manifest_record(lake_name, granule, status, n_pixels=0, refined_status="", message="",
                    processing_version=""):
    info = parse_granule_name(granule) or {}
    dt = info.get("datetime_utc")
    return {
        "lake_name": lake_name,
        "granule": granule,
        "datetime_utc": dt.strftime("%Y-%m-%d %H:%M:%S") if dt is not None else "",
        "status": status,
        "n_pixels": int(n_pixels),
        "refined_status": refined_status,
        "processing_version": processing_version,
        "processed_utc": datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S"),
        "message": str(message)[:300],
    }


def qc_mask(d, lake, cfg):
    """Quality, classification and elevation criteria shared by both assignment methods.
    `d` is a dict of arrays or a DataFrame; bounds are compared in NAVD88."""
    q = cfg["qc"]
    m = np.isin(np.asarray(d["classification"]), q["water_classification"])
    for col in ["geolocation_qual", "sig0_qual", "interferogram_qual", "classification_qual"]:
        m &= np.asarray(d[col]) <= q[f"{col}_max"]
    wse = np.asarray(d["wse_m_navd88"])
    m &= (wse >= lake["elev_min"]) & (wse <= lake["elev_max"])
    return m
