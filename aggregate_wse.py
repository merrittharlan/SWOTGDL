##################################################
## Authors: Merritt Harlan (USGS), Alexander Hernes (USGS)
## Last Updated: November 2026
#-------------------------------------------------
# Usage: python src/aggregate_wse.py [--config config.yaml]   (run from the repo root)
#-------------------------------------------------
# DESCRIPTION: Reads the per-granule pixel files in data/pixc/ and writes
# data/aggregated/wse_timeseries.parquet: one row per lake and overpass, with WSE
# statistics for each pixel set.
# Pixel sets: polygon_all (pixels in the lake polygon), polygon_qc (polygon + QC) and
# refined_qc (refined assignment + QC). Statistics per set: median, sd, weighted_mean
# (inverse-variance weighted mean), uncert_std and uncert_multilook (RiverObs
# aggregate.py, NaN if RiverObs isn't installed), and n_pix.
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
import re
from pathlib import Path

import numpy as np
import pandas as pd

from common import load_config

# RiverObs aggregate.py (optional) for uncertainties that account for correlated pixels:
#   git clone https://github.com/SWOTAlgorithms/RiverObs && export PYTHONPATH=RiverObs/src
try:
    import SWOTWater.aggregate as swot_agg
except ImportError:
    swot_agg = None
    print("SWOTWater (RiverObs) not found; uncert_* columns will be NaN.")

KEYS = ["lake_name", "datetime_utc", "granule", "cycle_number", "pass_number", "tile_name", "crid"]
PIXEL_SETS = {"polygon_all": "in_polygon", "polygon_qc": "in_polygon_qc", "refined_qc": "refined_qc"}
MULTILOOK_COLS = ["interferogram_real", "interferogram_imag", "power_plus_y", "power_minus_y",
                  "dlatitude_dphase", "dlongitude_dphase", "looks_to_efflooks"]


def aggregate_uncertainties(g):
    """(uncert_std, uncert_multilook) from RiverObs aggregate.py: the sample-std estimate
    scaled for correlated pixels (Eq. 14, Williams 2018) and the multilooked-coherence
    estimate (Eq. 7)."""
    if swot_agg is None or g.empty:
        return np.nan, np.nan
    h = g["wse_m_navd88"].to_numpy(float)
    height_std = (g["phase_noise_std"] * g["dheight_dphase"]).to_numpy(float)
    good = np.isfinite(h) & np.isfinite(height_std) & (height_std != 0)
    if not good.any():
        return np.nan, np.nan
    height_std = np.where(good, height_std, 1.0)
    rare = g["eff_num_rare_looks"].to_numpy(float)
    med = g["eff_num_medium_looks"].to_numpy(float)
    uncert_std = float(swot_agg.height_uncert_std(h, good, rare, med, height_std=height_std,
                                                  method="weight"))
    uncert_ml = np.nan
    if all(c in g.columns for c in MULTILOOK_COLS) and g[MULTILOOK_COLS].notna().all().all():
        res = swot_agg.height_with_uncerts(
            h, good, rare, med,
            g["interferogram_real"].to_numpy(float) + 1j * g["interferogram_imag"].to_numpy(float),
            g["power_plus_y"].to_numpy(float), g["power_minus_y"].to_numpy(float),
            float(g["looks_to_efflooks"].iloc[0]), g["dheight_dphase"].to_numpy(float),
            g["dlatitude_dphase"].to_numpy(float), g["dlongitude_dphase"].to_numpy(float),
            height_std=height_std, method="weight")
        uncert_ml = float(res[2])
    return uncert_std, uncert_ml


# Columns read from each pixel file (others are skipped to keep memory low)
READ_COLS = (KEYS + ["refined_status", "wse_m_navd88", "phase_noise_std", "dheight_dphase",
                     "eff_num_rare_looks", "eff_num_medium_looks"]
             + list(PIXEL_SETS.values()) + MULTILOOK_COLS)


def set_stats(sub, name):
    """WSE statistics for one overpass and one pixel set."""
    p = f"wse_m_navd88_{name}"
    out = {f"n_pix_{name}": len(sub)}
    if sub.empty:
        return out | {f"{p}_{k}": np.nan for k in
                      ["median", "sd", "weighted_mean", "uncert_std", "uncert_multilook"]}
    wse = sub["wse_m_navd88"].astype(float)
    w = 1.0 / (sub["phase_noise_std"].astype(float) * sub["dheight_dphase"].astype(float)) ** 2
    ok = np.isfinite(w) & (w > 0) & wse.notna()
    uncert_std, uncert_ml = aggregate_uncertainties(sub)
    out.update({
        f"{p}_median": wse.median(),
        f"{p}_sd": wse.std(),
        f"{p}_weighted_mean": (float((w[ok] * wse[ok]).sum() / w[ok].sum())
                               if ok.any() else np.nan),
        f"{p}_uncert_std": uncert_std,
        f"{p}_uncert_multilook": uncert_ml,
    })
    return out


def summarize_file(path):
    """One row of the time series from one pixel file (one lake, one overpass)."""
    import pyarrow.parquet as pq
    names = set(pq.read_schema(path).names)
    df = pd.read_parquet(path, columns=[c for c in READ_COLS if c in names])
    if df.empty:
        return None
    row = {k: df[k].iloc[0] for k in KEYS}
    row["datetime_utc"] = pd.Timestamp(row["datetime_utc"])
    for name, col in PIXEL_SETS.items():
        mask = df[col].fillna(False).astype(bool) if col in df else pd.Series(False, index=df.index)
        row.update(set_stats(df[mask], name))
    row["refined_status"] = df["refined_status"].iloc[0] if "refined_status" in df else None
    return row


def latest_versions(files):
    """NASA sometimes reprocesses a granule, giving a second file for the same overpass
    (e.g. ..._PID0_01 and ..._PID0_02). Keep only the newest version of each
    cycle/pass/tile; file names sort by CRID then counter."""
    newest = {}
    for f in sorted(files):
        m = re.search(r"PIXC_(\d{3}_\d{3}_\d{3}[LR])_", f.name)
        newest[m.group(1) if m else f.name] = f
    return sorted(newest.values())


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--config", default="config.yaml")
    args = ap.parse_args()
    cfg = load_config(args.config)
    data_root = Path(cfg["data_root"])

    rows = []
    for lake_dir in sorted(p for p in (data_root / "pixc").glob("*") if p.is_dir()):
        all_files = sorted(lake_dir.glob("*.parquet"))
        files = latest_versions(all_files)
        if not files:
            continue
        skipped = len(all_files) - len(files)
        print(f"{lake_dir.name}: {len(files)} overpasses"
              + (f" ({skipped} older granule versions skipped)" if skipped else ""))
        for f in files:
            try:
                row = summarize_file(f)
            except Exception as e:
                print(f"  could not summarize {f.name}: {type(e).__name__}: {e}")
                continue
            if row is not None:
                rows.append(row)
    if not rows:
        print("No pixel files found.")
        return

    ts = pd.DataFrame(rows)
    for name in PIXEL_SETS:
        ts[f"n_pix_{name}"] = ts[f"n_pix_{name}"].fillna(0).astype(int)
    ts["wse_m_navd88_refined_minus_polygon_qc_median"] = (
        ts["wse_m_navd88_refined_qc_median"] - ts["wse_m_navd88_polygon_qc_median"])
    ts = ts.sort_values(["lake_name", "datetime_utc"]).reset_index(drop=True)

    out_dir = data_root / "aggregated"
    out_dir.mkdir(parents=True, exist_ok=True)
    ts.to_parquet(out_dir / "wse_timeseries.parquet", index=False, compression="zstd")
    print(f"Wrote {len(ts):,} overpasses for {ts['lake_name'].nunique()} lakes to {out_dir}")


if __name__ == "__main__":
    main()
