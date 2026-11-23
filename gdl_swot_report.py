#!/usr/bin/env python3
"""
Glacier-dammed lake (GDL) SWOT monitoring report.

Builds a single, self-contained HTML page that lets forecasters quickly inspect
SWOT water-surface elevation (WSE) time series for any number of lakes.

Python port of Snow_Skilak_SWOT_GDL_Monitoring.Rmd, with:
  * lakes defined in a TOML config file (add a lake = add a [[lakes]] block)
  * optional download of input files from a ScienceBase item
  * a "forecaster overview" table and per-lake key numbers
  * a table of contents (sticky sidebar) listing every lake
  * a SWOT overpass table per lake (used, rejected, and predicted passes)

Examples
--------
    # From this repo's pipeline (what the GitLab job runs, from the repo root):
    python src/gdl_swot_report.py --config report.toml --all-lakes \\
        --timeseries data/aggregated/wse_timeseries.parquet --pixc-root data/pixc \\
        -o public/index.html

    python gdl_swot_report.py --config lakes.toml -o gdl_report.html
    python gdl_swot_report.py --config lakes.toml --sciencebase-item 1234abcd...
    python gdl_swot_report.py --config lakes.toml --lakes Snow --no-3d
    python gdl_swot_report.py --config lakes.toml --list-lakes
    python gdl_swot_report.py --config lakes.toml --publish   # build + push to GitHub Pages

Requires: Python >= 3.11, pandas, numpy, plotly, pyarrow.
"""

from __future__ import annotations

import argparse
import datetime as dt
import html
import json
import logging
import re
import shutil
import subprocess
import sys
import tomllib
import urllib.request
import zipfile
from pathlib import Path
from zoneinfo import ZoneInfo

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.offline import get_plotlyjs, get_plotlyjs_version

log = logging.getLogger("gdl_swot")

__version__ = "1.2.0"

# SWOT science orbit: every pass repeats about every 20.86 days.
SWOT_REPEAT_DAYS = 20.864

# Colorblind-safe palette (Okabe-Ito) used for years without a configured color.
PALETTE = ["#0072B2", "#E69F00", "#009E73", "#CC79A7", "#56B4E9",
           "#D55E00", "#999999", "#000000"]

# WSE column plotted and summarized everywhere: filtered (NAVD88) median.
WSE_COL = "wse_m_NAVD88_filtered_med"
SD_COL = "wse_m_NAVD88_filtered_sd"

STATUS_USED = "Used"
STATUS_OUTLIER = "Outlier (Hampel)"
STATUS_LOWPIX = "Too few pixels"
STATUS_NOWSE = "No WSE"
STATUS_PRED = "Predicted"

DEFAULT_REPORT = {
    "title": "Alaska Glacier-Dammed Lake SWOT Monitoring",
    "author": "U.S. Geological Survey",
    "display_timezone": "America/Anchorage",
    # Link shown in every lake's notes; a lake can override it with gdl_url.
    "gdl_url": "https://www.weather.gov/aprfc/gdlData?6",
    "gdl_link_text": "Glacial Dammed Lake Data",
    # How the quality assessment is cited on the page.
    "assessment_source": "Harlan et al., Monitoring Alaskan Glacier-Dammed Lakes "
                         "with the SWOT Satellite (manuscript in review)",
    "basemap": "Imagery",         # "Imagery" or "Topographic" (Esri)
}

DEFAULT_DATA = {
    "search_dirs": ["."],
    # One or more file names or glob patterns. Every match in every search
    # folder is read and combined (duplicates dropped). If nothing matches,
    # any CSV whose header has the required columns is used.
    "aggregated_csv": ["SWOT_lakes_aggregated_wse*.csv"],
}

DEFAULT_PUBLISH = {
    # GitHub repo served by GitHub Pages. --publish writes the report into a
    # local clone of it as `file`, commits, and pushes; the repo's Pages
    # workflow then redeploys the site.
    "repo_url": "https://github.com/merrittharlan/SWOTGDL.git",
    "repo_dir": "SWOTGDL",        # local clone (relative to the config file)
    "branch": "main",
    "file": "index.html",
    "commit_message": "Update SWOT GDL report ({date})",
}

DEFAULT_SB = {
    "item_id": "",
    "download_dir": "sciencebase_cache",
    "include_children": "auto",   # true / false / "auto" (if the item has children)
    "max_small_file_mb": 200,     # non-zip files up to this size are always fetched
}

# Column names the script uses. Map a different name in your files to these
# with a [columns] table in the config, e.g.  wse_m_NAVD88_filtered_med = "wse_fmed".
AGG_REQUIRED = ["lake_name", "datetime_UTC", "n_pix_filtered",
                "wse_m_NAVD88_filtered_med", "wse_m_NAVD88_filtered_sd"]

DEFAULT_SETTINGS = {
    "season_start_month": 6,
    "season_end_month": 11,
    "years": [],                 # empty = every year in the data
    "min_pixels": 10,
    "pixel_thresholds": [10, 25, 50, 100],
    "hampel_k": 3,
    "hampel_t0": 3.0,
    "min_obs_for_hampel": 5,
    "max_gap_days": 60,
    "trend_window_days": 30,
    "pre_release_window_days": 30,
    "stale_days": 21,
    "predict_days_ahead": 45,
    "max_3d_points": 5000,
    "max_3d_passes": 8,
    "pixc_min_year": 0,          # 0 = current year only
    "year_colors": {},
    # Only glacier-dammed lakes with these SWOT quality-assessment classes are
    # reported (--all-classes turns this off). Lakes not in the assessment
    # (e.g. gaged validation lakes) are left out unless a [[lakes]] entry sets
    # `classification`.
    "include_classes": ["Likely successful", "Potentially successful"],
}

# ---- SWOT quality assessment of NWS-monitored glacier-dammed lakes -------------
# From Table 3 of Harlan et al., "Monitoring Alaskan Glacier-Dammed Lakes with the
# SWOT Satellite" (manuscript in review). Only lakes listed here are glacier-dammed
# lakes (GDLs); the gaged validation lakes (Bradley, Eklutna, ...) are not.
# Keys are the manuscript names; matching is loose (see name_key) plus ALIASES.
QA_LIKELY = "Likely successful"
QA_POTENTIAL = "Potentially successful"
QA_FAILURE = "Failure"
QA_NODATA = "No data"
QA_RANK = {QA_LIKELY: 0, QA_POTENTIAL: 1, QA_FAILURE: 2, QA_NODATA: 3}

# name: (lat, lon, area_km2, pixel_sd_m, mean_gap_days, pass_diff_m,
#        pass_diff_from_unfiltered, mean_filtered_pixels, tri, classification)
GDL_ASSESSMENT = {
    "Bear":           (60.06721, -149.6106, 0.78, 17.7, 62.3, 102, True, 27, 9.06, QA_FAILURE),
    "Berg":           (60.40759, -143.8059, 17.8, 10.9, 6.7, 15.3, False, 4163, 4.00, QA_POTENTIAL),
    "Big Lake":       (61.60144, -146.6005, 3.4, 14.7, 12.7, 16.9, False, 257, 8.91, QA_POTENTIAL),
    "Blockade":       (60.95570, -152.4536, 11.7, 1.86, 51.8, 19.5, False, 1, 11.7, QA_FAILURE),
    "Frustration":    (61.48923, -151.9932, 2.29, 8.8, 10, 32.2, False, 391, 6.36, QA_POTENTIAL),
    "Hidden":         (61.56026, -143.0949, 0.66, 5.1, 28.1, 7.4, False, 9, 11.9, QA_FAILURE),
    "Iceberg":        (61.59540, -146.4780, 0.98, 10.7, 11.2, 28.4, False, 118, 10.7, QA_POTENTIAL),
    "Iceburg":        (60.75793, -142.9497, 0.90, 3.4, 9.4, 20.1, False, 1154, 4.79, QA_POTENTIAL),
    "Lake Erie":      (61.56808, -142.8838, 0.09, 0.8, 97.3, 61.5, True, 58, 8.41, QA_FAILURE),
    "Lake No Lake":   (58.88099, -133.8652, 1.82, 5.4, 14.2, 6.0, False, 20, 11.0, QA_LIKELY),
    "Mendenhall":     (58.45931, -134.4956, 0.50, None, None, 107, True, 0, 15.4, QA_NODATA),
    "Moose":          (60.47737, -149.0237, 0.06, 3.1, 34.9, 42.0, False, 37, 7.52, QA_FAILURE),
    "Nelchina North": (61.69307, -146.9274, 0.3, 7.2, 12.8, 42.7, False, 127, 9.11, QA_POTENTIAL),
    "Nelchina South": (61.66418, -146.8723, 1.57, 12.3, 10.0, 61.5, False, 752, 8.95, QA_POTENTIAL),
    "Sherman":        (60.56100, -145.1636, 2.29, 6.2, 21.9, 42.5, False, 76, 6.92, QA_FAILURE),
    "Skilak":         (60.20015, -149.9439, 0.37, 4.5, 7.7, 9.4, False, 186, 7.30, QA_LIKELY),
    "Snow":           (60.48265, -148.9259, 1.21, 7.7, 7, 9.1, False, 230, 10.3, QA_LIKELY),
    "Strandline":     (61.49589, -151.9020, 5.02, 10.0, 10.5, 44.3, False, 187, 8.61, QA_POTENTIAL),
    "Summit":         (56.18932, -130.0730, 1.29, 13.2, 26.7, 103, False, 7, 9.14, QA_FAILURE),
    "Tazlione":       (61.59944, -146.6089, 3.52, 12.9, 10.2, 36.0, False, 605, 8.94, QA_POTENTIAL),
    "Tulsequah":      (58.82559, -133.8114, 0.38, 11.0, 19.2, 35.7, False, 68, 7.74, QA_FAILURE),
    "Valdez":         (61.21912, -146.1371, 0.39, None, None, 25.8, True, 0, 8.69, QA_NODATA),
}
# Other spellings a data file might use for the same lake.
GDL_ALIASES = {
    "Mendenhall": ["Suicide Basin", "Suicide"],
    "Tazlione": ["Tazlina"],
    "Nelchina North": ["North Nelchina"],
    "Nelchina South": ["South Nelchina"],
    "Lake No Lake": ["No Lake"],
    "Big Lake": ["Big"],
}
# Table 1 thresholds: (key, label, unit, test, threshold text)
QA_CRITERIA = [
    ("area_km2", "Lake area", "km²", lambda v: v < 0.0625, "< 0.0625 km²"),
    ("gap_days", "Filtered temporal frequency (mean gap)", "days", lambda v: v > 15, "> 15 days → failure"),
    ("pixel_sd_m", "Filtered pixel standard deviation", "m", lambda v: v > 10, "> 10 m"),
    ("pass_diff_m", "Pass difference (between SWOT passes)", "m", lambda v: v > 10, "> 10 m"),
    ("pixels", "Average filtered pixel count", "px", lambda v: v < 10, "< 10 pixels"),
    ("tri", "ArcticDEM Terrain Roughness Index", "", lambda v: v > 13.9, "> 13.9"),
]
QA_SHORT = {"area_km2": "small area", "gap_days": "sparse observations",
            "pixel_sd_m": "high pixel SD", "pass_diff_m": "large pass difference",
            "pixels": "few pixels", "tri": "rough terrain"}
DATA_RELEASE_URL = "https://doi.org/10.5066/P1UR7DHM"
GDL_MAIN_URL = "https://www.weather.gov/aprfc/gdlMain"

# Esri basemap tiles (free use with attribution).
ESRI_TILES = {
    "Imagery": ("https://server.arcgisonline.com/ArcGIS/rest/services/World_Imagery/"
                "MapServer/tile/{z}/{y}/{x}",
                "Tiles &copy; Esri &mdash; Source: Esri, Maxar, Earthstar Geographics, "
                "and the GIS User Community"),
    "Topographic": ("https://server.arcgisonline.com/ArcGIS/rest/services/World_Topo_Map/"
                    "MapServer/tile/{z}/{y}/{x}",
                    "Tiles &copy; Esri &mdash; Esri, HERE, Garmin, FAO, NOAA, USGS, "
                    "and the GIS User Community"),
}
ESRI_LABELS = ("https://server.arcgisonline.com/ArcGIS/rest/services/Reference/"
               "World_Boundaries_and_Places/MapServer/tile/{z}/{y}/{x}")
LEAFLET_VERSION = "1.9.4"

PIXC_NAME_RE = re.compile(r"SWOT_L2_HR_PIXC_(\d{3})_(\d{3})_(\d{3}[LR])")
SB_ID_RE = re.compile(r"[0-9a-f]{24}")
SB_BASE = "https://www.sciencebase.gov/catalog"


# =============================================================================
# CONFIG
# =============================================================================

def load_config(path: Path) -> dict:
    with open(path, "rb") as f:
        raw = tomllib.load(f)

    cfg = {
        "report": {**DEFAULT_REPORT, **raw.get("report", {})},
        "data": {**DEFAULT_DATA, **raw.get("data", {})},
        "sciencebase": {**DEFAULT_SB, **raw.get("sciencebase", {})},
        "publish": {**DEFAULT_PUBLISH, **raw.get("publish", {})},
        "settings": {**DEFAULT_SETTINGS, **raw.get("settings", {})},
        "lakes": [],
        "base_dir": path.parent.resolve(),
        "columns": {str(k): str(v) for k, v in raw.get("columns", {}).items()},
    }
    agg = cfg["data"]["aggregated_csv"]
    cfg["data"]["aggregated_csv"] = [agg] if isinstance(agg, str) else list(agg)
    cfg["settings"]["year_colors"] = {
        int(k): v for k, v in cfg["settings"]["year_colors"].items()
    }

    for lk in raw.get("lakes", []):
        if "name" not in lk:
            raise ValueError(f"A [[lakes]] entry in {path} has no 'name'.")
        name = str(lk["name"])
        cfg["lakes"].append({
            "name": name,
            "display_name": lk.get("display_name", default_display_name(name)),
            "pixc_parquet": lk.get("pixc_parquet", f"{name}_complete.parquet"),
            "outline": lk.get("outline", f"{name}_lake.geojson"),
            "sciencebase_files": list(lk.get("sciencebase_files", [])),
            "notes": lk.get("notes", ""),
            "gdl_url": lk.get("gdl_url", ""),
            # Override / supply the quality-assessment class and spellings.
            "classification": lk.get("classification", ""),
            "aliases": [str(x) for x in lk.get("aliases", [])],
            "releases": [_to_date(d) for d in lk.get("releases", [])],
            # Outburst crest dates (e.g. from the NWS APRFC GDL event tables).
            "crests": [_to_date(d) for d in lk.get("crests", [])],
            "cap": [
                {"date": _to_date(c["date"]),
                 "elevation_ft": float(c["elevation_ft"])}
                for c in lk.get("cap", [])
            ],
        })
    return cfg


def default_display_name(name: str) -> str:
    """'Snow' -> 'Snow Lake', 'Big_Lake' -> 'Big Lake', 'Lake_No_Lake' -> 'Lake No Lake'."""
    n = str(name).replace("_", " ").strip()
    return n if "lake" in n.lower() else f"{n} Lake"


def default_lake(name: str) -> dict:
    return {
        "name": name, "display_name": default_display_name(name),
        "pixc_parquet": f"{name}_complete.parquet",
        "outline": f"{name}_lake.geojson", "sciencebase_files": [],
        "notes": "", "gdl_url": "", "classification": "", "aliases": [],
        "releases": [], "crests": [], "cap": [],
    }


def lake_assessment(lake: dict) -> dict | None:
    """Quality-assessment record for a lake (None if it is not an assessed GDL).

    A `classification` set in the lake's config entry overrides (or supplies)
    the class; metrics are only known for lakes in GDL_ASSESSMENT."""
    keys = {name_key(n) for n in [lake["name"], lake["display_name"], *lake.get("aliases", [])]}
    rec = None
    for nm, row in GDL_ASSESSMENT.items():
        if keys & {name_key(n) for n in [nm, *GDL_ALIASES.get(nm, [])]}:
            lat, lon, area, sd, gap, pdiff, p_unf, pix, tri, cls = row
            rec = {"name": nm, "lat": lat, "lon": lon, "area_km2": area,
                   "pixel_sd_m": sd, "gap_days": gap, "pass_diff_m": pdiff,
                   "pass_diff_unfiltered": p_unf, "pixels": pix, "tri": tri,
                   "classification": cls}
            break
    override = str(lake.get("classification") or "").strip()
    if override:
        match = next((c for c in QA_RANK if c.lower() == override.lower()), override)
        rec = {**(rec or {"name": lake["name"]}), "classification": match}
    if rec is not None:
        rec["flags"] = [k for k, _, _, test, _ in QA_CRITERIA
                        if rec.get(k) is not None and test(rec[k])]
    return rec


def _to_date(v) -> dt.date:
    if isinstance(v, dt.datetime):
        return v.date()
    if isinstance(v, dt.date):
        return v
    return dt.date.fromisoformat(str(v))


def slug(s: str) -> str:
    return re.sub(r"[^a-z0-9]+", "-", s.lower()).strip("-") or "lake"


def name_key(s: str) -> str:
    """'Snow Lake', 'snow_lake', 'SNOW' -> 'snow' (for loose name matching)."""
    k = re.sub(r"[^a-z0-9]", "", str(s).lower())
    k = re.sub(r"(glacierlake|glacialdammedlake|lake)$", "", k)
    return k or re.sub(r"[^a-z0-9]", "", str(s).lower())


def name_in(lake: str, text: str) -> bool:
    """True if the lake name appears in a file name/path as its own token."""
    toks = [t for t in re.split(r"[^a-z0-9]+", str(text).lower()) if t]
    k = name_key(lake)
    joined = re.sub(r"[^a-z0-9]", "", str(lake).lower())
    # match single token ("snow") or consecutive tokens ("suicide","basin")
    for i in range(len(toks)):
        acc = ""
        for t in toks[i:i + 4]:
            acc += t
            if acc in (k, joined):
                return True
    return False


# =============================================================================
# FILES: local search + ScienceBase download
# =============================================================================

class FileFinder:
    """Finds input files in an ordered list of folders (searched recursively)."""

    SKIP = {"__pycache__", ".git"}

    def __init__(self, dirs):
        seen, self.dirs = set(), []
        for d in dirs:
            if d and Path(d).is_dir() and Path(d).resolve() not in seen:
                seen.add(Path(d).resolve())
                self.dirs.append(Path(d))
        self._files = None

    def files(self) -> list[Path]:
        if self._files is None:
            out, seen = [], set()
            for d in self.dirs:
                for f in sorted(d.rglob("*")):
                    if (f.is_file() and not self.SKIP & set(f.parts)
                            and not f.name.endswith(".part")
                            and f.resolve() not in seen):
                        seen.add(f.resolve())
                        out.append(f)
            self._files = out
        return self._files

    def find(self, name: str | None) -> Path | None:
        hits = self.find_all(name)
        return hits[0] if hits else None

    def find_all(self, pattern: str | None) -> list[Path]:
        """Every file matching a name or glob, in search-folder order."""
        if not pattern:
            return []
        p = Path(pattern)
        if p.is_absolute():
            return sorted(p.parent.glob(p.name))
        out = []
        for d in self.dirs:
            direct = sorted(d.glob(pattern))
            out += [f for f in direct if f.is_file()]
        pat = p.name.lower()
        import fnmatch
        out += [f for f in self.files() if fnmatch.fnmatch(f.name.lower(), pat)]
        return list(dict.fromkeys(out))

    def lake_files(self, lake: str, exts: tuple[str, ...]) -> list[Path]:
        """Files with one of `exts` whose name or parent folder names the lake."""
        out = []
        for f in self.files():
            if f.suffix.lower() not in exts:
                continue
            rel = f.as_posix()
            if name_in(lake, f.stem) or any(name_in(lake, part) for part in f.parent.parts[-3:]):
                out.append(f)
        return out


def csv_header(path: Path) -> list[str]:
    try:
        with open(path, newline="", encoding="utf-8-sig", errors="replace") as fh:
            return [h.strip().strip('"') for h in fh.readline().split(",")]
    except OSError:
        return []


def parse_sb_id(text: str) -> str:
    m = SB_ID_RE.search(text.lower())
    if not m:
        raise ValueError(f"Could not find a ScienceBase item id in '{text}'.")
    return m.group(0)


def _http_json(url: str) -> dict:
    req = urllib.request.Request(url, headers={
        "User-Agent": f"gdl-swot-report/{__version__}",
        "Accept": "application/json",
    })
    with urllib.request.urlopen(req, timeout=60) as r:
        return json.load(r)


def _sb_files(item: dict) -> list[dict]:
    files = list(item.get("files") or [])
    for facet in item.get("facets") or []:
        files += facet.get("files") or []
    out, seen = [], set()
    for f in files:
        url = f.get("downloadUri") or f.get("url")
        if f.get("name") and url and f["name"] not in seen:
            seen.add(f["name"])
            out.append({"name": f["name"], "url": url, "size": f.get("size"),
                        "item": item.get("id", ""), "item_title": item.get("title", "")})
    return out


def sciencebase_catalog(item_id: str, include_children) -> tuple[dict, list[dict]]:
    """Return (item json, list of files across the item and, optionally, children)."""
    item = _http_json(f"{SB_BASE}/item/{item_id}?format=json")
    items = [item]
    want_kids = (include_children is True or
                 (str(include_children).lower() == "auto" and item.get("hasChildren")))
    if want_kids:
        kids = _http_json(f"{SB_BASE}/items?parentId={item_id}"
                          "&format=json&max=1000&fields=id,title")
        for k in kids.get("items", []):
            items.append(_http_json(f"{SB_BASE}/item/{k['id']}?format=json"))
    entries = []
    for it in items:
        entries += _sb_files(it)
    return item, entries


def select_sb_files(entries: list[dict], lakes: list[dict], all_lakes: bool,
                    want_pixc: bool, max_small_mb: float) -> tuple[list[dict], list[str]]:
    """Pick what to download.

    * every non-zip file up to `max_small_mb` (the master CSV, GeoJSON, metadata)
    * zip packages only when pixel clouds are wanted, and then only the ones whose
      name contains a configured lake name (or that a lake lists explicitly in
      `sciencebase_files`). --all-lakes takes every zip.
    """
    msgs = []
    picked = []
    explicit = {n.lower() for lk in lakes for n in lk.get("sciencebase_files", [])}
    zips = [e for e in entries if e["name"].lower().endswith(".zip")]
    for e in entries:
        nm = e["name"].lower()
        size_mb = (int(e["size"]) / 1e6) if e.get("size") else 0
        if nm in explicit:
            picked.append(e)
        elif not nm.endswith(".zip") and size_mb <= max_small_mb:
            picked.append(e)
    if want_pixc:
        for e in zips:
            if e in picked:
                continue
            if all_lakes or any(name_in(lk["name"], e["name"]) or
                                name_in(lk["display_name"], e["name"]) for lk in lakes):
                picked.append(e)
        for lk in lakes:
            has = any(name_in(lk["name"], e["name"]) or name_in(lk["display_name"], e["name"])
                      or e["name"].lower() in {n.lower() for n in lk.get("sciencebase_files", [])}
                      for e in picked if e["name"].lower().endswith(".zip"))
            if zips and not has and not all_lakes:
                msgs.append(f"No ScienceBase zip matched lake '{lk['name']}'. Set "
                            f"sciencebase_files = [\"<zip name>\"] for it in the config. "
                            f"Zips on the item: {', '.join(z['name'] for z in zips[:12])}"
                            + (" …" if len(zips) > 12 else ""))
    return picked, msgs


def download_sciencebase(entries: list[dict], root: Path, refresh: bool) -> list[Path]:
    """Download entries into root/. Files already present with the same size are
    skipped unless `refresh`. Zip archives are extracted to root/<zip stem>/."""
    root.mkdir(parents=True, exist_ok=True)
    got = []
    for e in entries:
        target = root / Path(e["name"]).name
        if (target.exists() and not refresh and
                (not e.get("size") or target.stat().st_size == int(e["size"]))):
            log.info("  cached  %s", e["name"])
        else:
            mb = f" ({int(e['size']) / 1e6:.1f} MB)" if e.get("size") else ""
            log.info("  get     %s%s", e["name"], mb)
            tmp = target.with_name(target.name + ".part")
            req = urllib.request.Request(
                e["url"], headers={"User-Agent": f"gdl-swot-report/{__version__}"})
            with urllib.request.urlopen(req, timeout=900) as r, open(tmp, "wb") as f:
                while chunk := r.read(1 << 20):
                    f.write(chunk)
            tmp.replace(target)
        if target.suffix.lower() == ".zip":
            out_dir = root / target.stem
            marker = out_dir / ".extracted"
            if refresh or not marker.exists() or marker.stat().st_mtime < target.stat().st_mtime:
                log.info("  unzip   %s", target.name)
                with zipfile.ZipFile(target) as z:
                    z.extractall(out_dir)
                marker.touch()
        got.append(target)
    (root / "_downloaded.json").write_text(json.dumps({
        "downloaded_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "files": [e["name"] for e in entries],
    }, indent=2))
    return got


def print_sb_listing(item: dict, entries: list[dict]) -> None:
    print(f"{item.get('title', '')}\n{SB_BASE}/item/{item.get('id', '')}")
    dates = {d.get("type"): d.get("dateString") for d in item.get("dates", []) or []}
    if dates:
        print("  " + ", ".join(f"{k}: {v}" for k, v in dates.items() if v))
    cur = None
    for e in entries:
        if e["item"] != cur:
            cur = e["item"]
            print(f"\n[{e['item']}] {e['item_title']}")
        size = f"{int(e['size']) / 1e6:9.1f} MB" if e.get("size") else "        ? MB"
        print(f"  {size}  {e['name']}")


# =============================================================================
# DATA PROCESSING
# =============================================================================

def _rename(df: pd.DataFrame, columns: dict) -> pd.DataFrame:
    inv = {v: k for k, v in columns.items() if v in df.columns and k not in df.columns}
    return df.rename(columns=inv)


def find_aggregated(finder: FileFinder, patterns: list[str], columns: dict,
                    sb_dir: Path | None = None) -> list[Path]:
    """CSV files to combine: every name/glob match, plus any CSV with the right
    header inside the ScienceBase download (whatever the release calls it).
    ScienceBase files come first so local files win on duplicate passes."""
    need = set(AGG_REQUIRED)
    def ok(f):
        hdr = set(csv_header(f))
        hdr |= {k for k, v in columns.items() if v in hdr}
        return need <= hdr
    hits = []
    for pat in patterns:
        hits += [f for f in finder.find_all(pat) if f.suffix.lower() == ".csv"]
    hits = list(dict.fromkeys(hits))
    good = [f for f in hits if ok(f)]
    for f in set(hits) - set(good):
        log.warning("Skipping %s: header lacks %s", f, sorted(need - set(csv_header(f))))
    if sb_dir is not None and sb_dir.is_dir():
        sb = [f for f in sorted(sb_dir.rglob("*.csv")) if ok(f)]
        good = list(dict.fromkeys(sb + [f for f in good if sb_dir not in f.parents]))
    if good:
        return good
    # Nothing matched by name: take any CSV with the right header.
    auto = [f for f in finder.files() if f.suffix.lower() == ".csv" and ok(f)]
    if auto:
        log.info("Aggregated CSV found by header: %s", ", ".join(map(str, auto)))
    return auto


def read_aggregated_many(paths: list[Path], columns: dict) -> pd.DataFrame:
    """Read and combine aggregated CSVs; later files win on duplicate passes."""
    frames = []
    for p in paths:
        d = _rename(pd.read_csv(p), columns)
        d["_source"] = p.name
        frames.append(d)
    df = read_aggregated(pd.concat(frames, ignore_index=True))
    if len(paths) > 1:
        # Same lake (loose name match) and same overpass minute = same pass.
        k = [df["lake_name"].map(name_key), df["time"].dt.round("min")]
        n0 = len(df)
        df = df[~pd.DataFrame({"a": k[0], "b": k[1]}).duplicated(keep="last")]
        log.info("Combined %d CSVs: %d rows (%d duplicate passes dropped)",
                 len(paths), len(df), n0 - len(df))
    return df.sort_values(["lake_name", "time"]).reset_index(drop=True)


def unify_lake_names(agg: pd.DataFrame, lakes: list[dict]) -> pd.DataFrame:
    """Give every spelling of a lake one name ('Snow Lake', 'Snow' -> config name)."""
    canon = {name_key(lk["name"]): lk["name"] for lk in lakes}
    for n in agg["lake_name"].value_counts().index:
        canon.setdefault(name_key(n), n)
    agg["lake_name"] = agg["lake_name"].map(lambda n: canon[name_key(n)])
    return agg


def read_aggregated(df: pd.DataFrame) -> pd.DataFrame:
    missing = set(AGG_REQUIRED) - set(df.columns)
    if missing:
        raise ValueError(f"Aggregated data is missing columns {sorted(missing)}; "
                         f"columns present: {sorted(df.columns)}. Map names with a "
                         "[columns] table in the config.")

    df["lake_name"] = df["lake_name"].astype(str)
    df["time"] = (pd.to_datetime(df["datetime_UTC"], utc=True,
                                 errors="coerce", format="mixed")
                  .dt.tz_convert(None))
    df = df.dropna(subset=["time"]).copy()

    url = df["pixc_url"].astype(str) if "pixc_url" in df else pd.Series("", index=df.index)
    parts = url.str.extract(PIXC_NAME_RE)
    df["cycle"] = parts[0]
    df["pass"] = parts[1]
    df["tile"] = parts[2]
    # Fallback used by the original R script: characters 124-127 of the URL.
    no_pass = df["pass"].isna()
    df.loc[no_pass, "pass"] = url[no_pass].str.slice(123, 127).replace("", np.nan)
    return df.sort_values(["lake_name", "time"]).reset_index(drop=True)


def hampel_flags(x: np.ndarray, k: int, t0: float) -> np.ndarray:
    """Centered Hampel filter; True where a value is an outlier.

    The window shrinks at the ends of the series. Windows whose MAD is zero
    never flag (avoids flagging every value in a perfectly flat window).
    """
    n = len(x)
    flags = np.zeros(n, dtype=bool)
    for i in range(n):
        w = x[max(0, i - k): i + k + 1]
        med = np.nanmedian(w)
        mad = 1.4826 * np.nanmedian(np.abs(w - med))
        if mad > 0 and abs(x[i] - med) > t0 * mad:
            flags[i] = True
    return flags


def classify_lake(raw: pd.DataFrame, s: dict) -> tuple[pd.DataFrame, list[str]]:
    """Return every in-window pass for one lake with a status column."""
    notes = []
    d = raw.copy()
    m = d["time"].dt.month
    d = d[(m >= s["season_start_month"]) & (m <= s["season_end_month"])]
    if s["years"]:
        d = d[d["time"].dt.year.isin(s["years"])]
    d = d.sort_values("time").copy()

    d["status"] = STATUS_USED
    no_wse = d[WSE_COL].isna()
    low = ~no_wse & (d["n_pix_filtered"] < s["min_pixels"])
    d.loc[no_wse, "status"] = STATUS_NOWSE
    d.loc[low, "status"] = STATUS_LOWPIX

    cand = d["status"] == STATUS_USED
    if cand.sum() >= s["min_obs_for_hampel"]:
        x = d.loc[cand, "wse_m_NAVD88_filtered_med"].to_numpy(float)
        flags = hampel_flags(x, int(s["hampel_k"]), float(s["hampel_t0"]))
        d.loc[d.index[cand.to_numpy()][flags], "status"] = STATUS_OUTLIER
    elif cand.sum() > 0:
        notes.append(f"Only {int(cand.sum())} valid observations — outlier "
                     "filter not applied.")

    d["date"] = d["time"].dt.normalize()
    d["year"] = d["time"].dt.year
    d["doy"] = d["time"].dt.dayofyear
    return d, notes


def adjust_cap(cap: list[dict], used: pd.DataFrame) -> pd.DataFrame:
    """CAP elevations shifted by the annual mean CAP-to-nearest-SWOT offset."""
    if not cap:
        return pd.DataFrame(columns=["date", "year", "doy", "cap_m", "cap_adj_m"])
    c = pd.DataFrame(cap)
    c["date"] = pd.to_datetime(c["date"])
    c["year"] = c["date"].dt.year
    c["doy"] = c["date"].dt.dayofyear
    c["cap_m"] = c["elevation_ft"] * 0.3048

    offsets = []
    for _, r in c.iterrows():
        same = used[used["year"] == r["year"]]
        if same.empty:
            offsets.append(np.nan)
            continue
        i = (same["date"] - r["date"]).abs().idxmin()
        offsets.append(same.at[i, WSE_COL] - r["cap_m"])
    c["offset"] = offsets
    c["mean_offset"] = c.groupby("year")["offset"].transform("mean")
    c["cap_adj_m"] = c["cap_m"] + c["mean_offset"]
    return c


def predict_overpasses(raw: pd.DataFrame, now: pd.Timestamp,
                       days_ahead: int) -> pd.DataFrame:
    """Project future overpasses from the last time each pass number was seen."""
    rows = []
    d = raw.dropna(subset=["pass"])
    if d.empty:
        return pd.DataFrame(columns=["time", "pass", "tile"])
    horizon = now + pd.Timedelta(days=days_ahead)
    rep = pd.Timedelta(days=SWOT_REPEAT_DAYS)
    for p, g in d.groupby("pass"):
        last = g["time"].max()
        tiles = ", ".join(sorted(set(g.loc[g["time"] == last, "tile"].dropna())))
        n = max(1, int(np.ceil((now - last) / rep)))
        t = last + n * rep
        while t <= horizon:
            rows.append({"time": t, "pass": p, "tile": tiles})
            t += rep
    return pd.DataFrame(rows).sort_values("time") if rows else \
        pd.DataFrame(columns=["time", "pass", "tile"])


def lake_summary(used: pd.DataFrame, lake: dict, s: dict,
                 now: pd.Timestamp, predicted: pd.DataFrame) -> dict:
    out = {"n_used": len(used)}
    out["next_pass"] = predicted.iloc[0] if len(predicted) else None
    if used.empty:
        return out

    wse = WSE_COL
    last = used.iloc[-1]
    out["latest"] = last
    out["age_days"] = (now - last["time"]) / pd.Timedelta(days=1)

    if len(used) > 1:
        prev = used.iloc[-2]
        days = (last["time"] - prev["time"]) / pd.Timedelta(days=1)
        out["delta_prev"] = last[wse] - prev[wse]
        out["delta_prev_days"] = days
        out["prev_time"] = prev["time"]

    win = used[used["time"] >= last["time"] - pd.Timedelta(days=s["trend_window_days"])]
    if len(win) >= 3:
        t = (win["time"] - win["time"].min()) / pd.Timedelta(days=1)
        out["trend_m_per_day"] = float(np.polyfit(t, win[wse], 1)[0])
        out["trend_n"] = len(win)

    yr = int(last["year"])
    this_year = used[used["year"] == yr]
    pk = this_year.loc[this_year[wse].idxmax()]
    out["season_peak"] = (pk[wse], pk["time"])

    prior = used[used["year"] == yr - 1].sort_values("doy")
    if len(prior) >= 2 and prior["doy"].min() <= last["doy"] <= prior["doy"].max():
        ref = np.interp(last["doy"], prior["doy"], prior[wse])
        out["vs_prior_year"] = (last[wse] - ref, yr - 1)

    past_rel = sorted(r for r in lake["releases"] if pd.Timestamp(r) <= now)
    for rel in reversed(past_rel):
        rel_t = pd.Timestamp(rel)
        pre = used[(used["time"] <= rel_t + pd.Timedelta(days=1)) &
                   (used["time"] >= rel_t - pd.Timedelta(days=s["pre_release_window_days"]))]
        if len(pre):
            out["vs_release"] = (last[wse] - pre[wse].max(), pre[wse].max(), rel)
            break
    if past_rel:
        out["last_release"] = past_rel[-1]
    return out


def _feature_names_lake(f: dict, lake: str) -> bool:
    props = f.get("properties") or {}
    return any(isinstance(v, str) and name_key(v) == name_key(lake) for v in props.values())


def find_outline(finder: FileFinder, lake: dict) -> tuple[Path | None, str | None]:
    """Return (geojson path, lake to select inside it or None for whole file)."""
    p = finder.find(lake["outline"])
    if p:
        return p, None
    own = [f for f in finder.lake_files(lake["name"], (".geojson", ".json"))
           if not f.name.startswith("_")]
    if own:
        # If the file holds several lakes, use only the features that name this one
        try:
            feats = json.loads(own[0].read_text()).get("features") or []
        except (ValueError, OSError, AttributeError):
            feats = []
        named = any(_feature_names_lake(ft, lake["name"]) for ft in feats)
        return own[0], (lake["name"] if named else None)
    # A shared GeoJSON (e.g. the release's all-lakes file): pick features by name.
    for f in finder.files():
        if f.suffix.lower() not in (".geojson", ".json") or f.name.startswith("_"):
            continue
        if f.stat().st_size > 200e6:
            continue
        try:
            gj = json.loads(f.read_text())
        except (ValueError, OSError):
            continue
        feats = gj.get("features") if isinstance(gj, dict) else None
        if feats and any(_feature_names_lake(ft, lake["name"]) and
                         (ft.get("geometry") or {}).get("type") in ("Polygon", "MultiPolygon")
                         for ft in feats):
            return f, lake["name"]
    return None, None


def read_outline(path: Path | None, select: str | None = None) -> list[tuple[list, list]]:
    """Return polygon rings as (lon list, lat list). GeoJSON is lon/lat WGS84.
    With `select`, only features whose properties name that lake are used."""
    if path is None:
        return []
    gj = json.loads(path.read_text())
    crs = (gj.get("crs") or {}).get("properties", {}).get("name", "")
    if crs and not any(c in crs for c in ("4326", "CRS84")):
        log.warning("%s declares CRS %s; GeoJSON outlines are expected in "
                    "lon/lat (EPSG:4326) and are plotted as-is.", path.name, crs)
    feats = gj["features"] if gj.get("type") == "FeatureCollection" else [gj]
    if select:
        feats = [f for f in feats if _feature_names_lake(f, select)]
    rings = []
    for f in feats:
        g = f.get("geometry", f)
        if g is None:
            continue
        polys = [g["coordinates"]] if g["type"] == "Polygon" else \
            g["coordinates"] if g["type"] == "MultiPolygon" else []
        for poly in polys:
            for ring in poly:
                rings.append(([p[0] for p in ring], [p[1] for p in ring]))
    return rings


PIXC_COLS = ["datetime_UTC", "filtered_boolean", "elevation_m_NAVD88",
             "latitude_WGS84", "longitude_WGS84", "pixel_area_sqm",
             "water_classification", "cross_track_dist_m"]


PIXC_REQUIRED = {"datetime_UTC", "elevation_m_NAVD88", "latitude_WGS84", "longitude_WGS84"}


def _pixc_cols(path: Path, columns: dict) -> list[str] | None:
    """Source column names to read from a pixel-cloud file, or None if not one."""
    if path.suffix.lower() == ".parquet":
        import pyarrow.parquet as pq
        try:
            names = set(pq.read_schema(path).names)
        except Exception:
            return None
    else:
        names = set(csv_header(path))
    std = {c: (c if c in names else columns.get(c)) for c in PIXC_COLS}
    std = {k: v for k, v in std.items() if v in names}
    return list(std.values()) if PIXC_REQUIRED <= set(std) else None


def find_pixc(finder: FileFinder, lake: dict, columns: dict) -> list[Path]:
    p = finder.find(lake["pixc_parquet"])
    if p:
        return [p]
    cands = finder.lake_files(lake["name"], (".parquet", ".csv"))
    return [f for f in cands if _pixc_cols(f, columns) is not None]


def read_pixc(paths: list[Path], s: dict, now: pd.Timestamp, columns: dict) -> pd.DataFrame:
    if not paths:
        return pd.DataFrame()
    import pyarrow.parquet as pq
    frames = []
    for path in paths:
        cols = _pixc_cols(path, columns)
        if cols is None:
            continue
        if path.suffix.lower() == ".parquet":
            fb = columns.get("filtered_boolean", "filtered_boolean")
            filters = [(fb, "==", True)] if fb in cols else None
            d = pq.read_table(path, columns=cols, filters=filters).to_pandas()
        else:
            d = pd.read_csv(path, usecols=cols)
        frames.append(_rename(d, columns))
    if not frames:
        return pd.DataFrame()
    p = pd.concat(frames, ignore_index=True)
    if "filtered_boolean" in p:
        fb = p["filtered_boolean"]
        if fb.dtype == object:
            fb = fb.astype(str).str.lower().isin(["true", "1", "t", "yes"])
        p = p[fb.astype(bool)]
    p = p.dropna(subset=["elevation_m_NAVD88", "latitude_WGS84", "longitude_WGS84"])
    p["time"] = (pd.to_datetime(p["datetime_UTC"], utc=True, errors="coerce",
                                format="mixed").dt.tz_convert(None))
    p = p.dropna(subset=["time"])
    m = p["time"].dt.month
    min_year = s["pixc_min_year"] or now.year
    p = p[(m >= s["season_start_month"]) & (m <= s["season_end_month"])
          & (p["time"].dt.year >= min_year)]
    p["date"] = p["time"].dt.normalize()
    return p


# ---- inputs from this repo's pipeline (export_pixc.py / aggregate_wse.py) ----
# Pixel set -> boolean column in the per-granule pixel files
WSE_SETS = {"polygon_qc": "in_polygon_qc", "refined_qc": "refined_qc", "polygon_all": "in_polygon"}


def read_timeseries(path: Path, wse_set: str) -> pd.DataFrame:
    """data/aggregated/wse_timeseries.parquet, renamed to the columns this report uses."""
    d = pd.read_parquet(path)
    p = f"wse_m_navd88_{wse_set}"
    df = pd.DataFrame({
        "lake_name": d["lake_name"].astype(str),
        "datetime_UTC": d["datetime_utc"],
        "n_pix_filtered": d[f"n_pix_{wse_set}"],
        WSE_COL: d[f"{p}_median"],
        SD_COL: d[f"{p}_sd"],
        "pixc_url": d["granule"].astype(str),   # cycle/pass/tile are parsed from the name
    })
    return read_aggregated(df)


def pipeline_pixc_files(root: Path, lake: str, s: dict, now: pd.Timestamp) -> list[Path]:
    """Per-granule pixel files for one lake (data/pixc/<lake>/), from pixc_min_year on.
    The year comes from the granule name, so older files are never opened."""
    dirs = [d for d in Path(root).iterdir() if d.is_dir() and name_key(d.name) == name_key(lake)]
    min_year = s["pixc_min_year"] or now.year
    newest: dict[str, Path] = {}   # newest version of each cycle/pass/tile (names sort by version)
    for d in dirs:
        for f in sorted(d.glob("*.parquet")):
            m = re.search(r"_(\d{4})\d{4}T\d{6}_", f.name)
            if m is None or int(m.group(1)) >= min_year:
                k = re.search(r"PIXC_(\d{3}_\d{3}_\d{3}[LR])_", f.name)
                newest[k.group(1) if k else f.name] = f
    return sorted(newest.values())


def read_pipeline_pixc(paths: list[Path], s: dict, wse_set: str) -> pd.DataFrame:
    """Pixels of the chosen set, renamed to the columns fig_pixc() uses. The refined set
    is drawn at its refined geolocation."""
    if not paths:
        return pd.DataFrame()
    import pyarrow.parquet as pq
    mask = WSE_SETS[wse_set]
    lat, lon = (("latitude_refined", "longitude_refined") if wse_set == "refined_qc"
                else ("latitude", "longitude"))
    rename = {"wse_m_navd88": "elevation_m_NAVD88", lat: "latitude_WGS84",
              lon: "longitude_WGS84", "pixel_area": "pixel_area_sqm",
              "classification": "water_classification", "cross_track": "cross_track_dist_m"}
    frames = []
    for f in paths:
        names = set(pq.read_schema(f).names)
        cols = [c for c in ["datetime_utc", mask, *rename] if c in names]
        if mask not in cols:
            continue
        frames.append(pq.read_table(f, columns=cols, filters=[(mask, "==", True)]).to_pandas())
    if not frames:
        return pd.DataFrame()
    p = pd.concat(frames, ignore_index=True).rename(columns=rename)
    p = p.dropna(subset=["elevation_m_NAVD88", "latitude_WGS84", "longitude_WGS84"])
    p["time"] = pd.to_datetime(p["datetime_utc"])
    m = p["time"].dt.month
    p = p[(m >= s["season_start_month"]) & (m <= s["season_end_month"])].copy()
    p["date"] = p["time"].dt.normalize()
    return p


# =============================================================================
# FIGURES
# =============================================================================

def year_color(year: int, s: dict) -> str:
    yc = s["year_colors"]
    if year in yc:
        return yc[year]
    return PALETTE[year % len(PALETTE)]


def _rgba(hex_color: str, a: float) -> str:
    h = hex_color.lstrip("#")
    r, g, b = (int(h[i:i + 2], 16) for i in (0, 2, 4))
    return f"rgba({r},{g},{b},{a})"


def _threshold_labels(s: dict) -> list[tuple[str, int]]:
    out = []
    for t in sorted(set(int(t) for t in s["pixel_thresholds"])):
        if t <= s["min_pixels"]:
            continue
        out.append((f"≥{t} pixels", t))
    return [(f"All (≥{s['min_pixels']} px)", int(s["min_pixels"]))] + out


BASE_LAYOUT = dict(
    template="plotly_white",
    font=dict(family="Arial, Helvetica, sans-serif", size=12, color="#222"),
    hovermode="closest",
    margin=dict(l=60, r=20, t=40, b=40),
    hoverlabel=dict(font_size=12),
)


def _threshold_menu(labels, groups, n_traces, always):
    buttons = []
    for i, (lab, _) in enumerate(labels):
        vis = [False] * n_traces
        for j in groups[i]:
            vis[j] = True
        for j in always:
            vis[j] = True
        buttons.append(dict(method="restyle", label=lab, args=[{"visible": vis}]))
    return [dict(type="dropdown", direction="down", x=1, xanchor="right",
                 y=1.13, yanchor="bottom", buttons=buttons, showactive=True,
                 bgcolor="white", bordercolor="#cbd2d9", font=dict(size=11))]


CREST_STYLE = dict(color="#3d3d3d", dash="dashdot", width=1.5)


def add_event_lines(fig, xs, labels, lo, hi, name, line, legendgroup=None,
                    showlegend=True) -> list[int]:
    """Vertical lines drawn as one hoverable trace (x = positions). Returns trace ids.
    meta='event' tells the page's auto-y code to ignore them."""
    if not xs:
        return []
    x, y, txt = [], [], []
    for xv, lab in zip(xs, labels):
        x += [xv, xv, None]
        y += [lo, hi, None]
        txt += [lab, lab, ""]
    fig.add_trace(go.Scatter(
        x=x, y=y, mode="lines", name=name, line=line, meta="event",
        legendgroup=legendgroup or name, showlegend=showlegend,
        text=txt, hovertemplate="%{text}<extra></extra>"))
    return [len(fig.data) - 1]


def _event_span(used):
    lo, hi = used[WSE_COL].min(), used[WSE_COL].max()
    pad = max(5.0, 0.5 * (hi - lo))
    return lo - pad, hi + pad


def fig_timeseries(used, cap, lake, s, now) -> go.Figure:
    wse, sd = WSE_COL, SD_COL
    fig = go.Figure()
    labels = _threshold_labels(s)
    groups = []

    for ti, (_, thr) in enumerate(labels):
        d = used[used["n_pix_filtered"] >= thr].sort_values("time").copy()
        gap = d["time"].diff() / pd.Timedelta(days=1)
        d["seg"] = ((gap.isna()) | (gap > s["max_gap_days"])
                    | (d["year"].diff() != 0)).cumsum()
        idx = []
        shown_years = set()
        # Faint dotted line across gaps longer than max_gap_days within a season,
        # so no observation looks orphaned.
        prev = d.shift(1)
        jumps = d[(d["seg"] != prev["seg"]) & (d["year"] == prev["year"])].index
        for j in jumps:
            k = d.index.get_loc(j)
            a_, b_ = d.iloc[k - 1], d.iloc[k]
            col = year_color(int(b_["year"]), s)
            fig.add_trace(go.Scatter(
                x=[a_["time"], b_["time"]], y=[a_[wse], b_[wse]], mode="lines",
                line=dict(color=col, width=1, dash="dot"), hoverinfo="skip",
                showlegend=False, legendgroup=str(int(b_["year"])), visible=ti == 0,
                name=f"gap > {s['max_gap_days']} d"))
            idx.append(len(fig.data) - 1)
        for _, g in d.groupby("seg"):
            yr = int(g["year"].iloc[0])
            col = year_color(yr, s)
            sdv = g[sd].fillna(0).to_numpy()
            y = g[wse].to_numpy()
            x = list(g["time"])
            fig.add_trace(go.Scatter(
                x=x + x[::-1], y=list(y + sdv) + list((y - sdv)[::-1]),
                mode="lines",   # Plotly adds markers to short traces otherwise
                fill="toself", fillcolor=_rgba(col, 0.15), line=dict(width=0),
                hoverinfo="skip", showlegend=False, legendgroup=str(yr),
                visible=ti == 0, name=f"{yr} ±1 SD"))
            idx.append(len(fig.data) - 1)
            fig.add_trace(go.Scatter(
                x=x, y=y, mode="lines+markers", name=str(yr),
                legendgroup=str(yr), showlegend=yr not in shown_years,
                line=dict(color=col, width=2), marker=dict(size=7, color=col),
                error_y=dict(type="data", array=sdv, visible=True, color=col,
                             thickness=1.2, width=3),
                customdata=np.stack([
                    g["time"].dt.strftime("%Y-%m-%d %H:%M"), g[sd].round(2),
                    g["n_pix_filtered"], g["pass"].fillna("—")], axis=-1),
                hovertemplate=("<b>%{customdata[0]} UTC</b><br>"
                               "WSE %{y:.2f} m ± %{customdata[1]} m<br>"
                               "%{customdata[2]} pixels · pass %{customdata[3]}"
                               "<extra></extra>"),
                visible=ti == 0))
            idx.append(len(fig.data) - 1)
            shown_years.add(yr)
        groups.append(idx)

    always = []
    cap = cap.dropna(subset=["cap_adj_m"]) if len(cap) else cap
    for yr, g in (cap.groupby("year") if len(cap) else []):
        col = year_color(int(yr), s)
        fig.add_trace(go.Scatter(
            x=g["date"], y=g["cap_adj_m"], mode="markers", name=f"CAP {yr}",
            marker=dict(size=11, symbol="diamond", color=col,
                        line=dict(color="black", width=1.5)),
            customdata=np.stack([g["date"].dt.strftime("%Y-%m-%d"),
                                 g["elevation_ft"]], axis=-1),
            hovertemplate=("<b>CAP %{customdata[0]}</b><br>Adjusted WSE "
                           "%{y:.2f} m<br>Reported %{customdata[1]} ft"
                           "<extra></extra>")))
        always.append(len(fig.data) - 1)

    for rel in lake["releases"]:
        fig.add_vline(x=pd.Timestamp(rel).timestamp() * 1000, line_dash="dash",
                      line_width=1, line_color="#52606d",
                      annotation_text=f"release {rel:%Y-%m-%d}",
                      annotation_font_size=10, annotation_position="top left")
    if lake.get("crests") and len(used):
        lo, hi = _event_span(used)
        always += add_event_lines(
            fig, [pd.Timestamp(c) for c in lake["crests"]],
            [f"Crest {c:%Y-%m-%d}" for c in lake["crests"]], lo, hi, "Crest", CREST_STYLE)
    fig.add_vline(x=now.timestamp() * 1000, line_dash="dot", line_width=1,
                  line_color="#9aa5b1", annotation_text="today",
                  annotation_font_size=10, annotation_position="top right")

    # Default view: every year in the data (use the range buttons to zoom in).
    xrange = None

    fig.update_layout(
        **{**BASE_LAYOUT, "margin": dict(l=60, r=20, t=95, b=30)}, height=500,
        updatemenus=_threshold_menu(labels, groups, len(fig.data), always),
        legend=dict(orientation="h", x=0, y=1.2, yanchor="bottom",
                    font=dict(size=11)),
        yaxis=dict(title="WSE (m NAVD88)", tickformat=".1f"),
        xaxis=dict(
            range=xrange,
            rangeslider=dict(visible=True, thickness=0.08),
            rangeselector=dict(
                x=0, y=1.02, yanchor="bottom", font=dict(size=11),
                bgcolor="#e8edf3", activecolor="#b9cde5",
                buttons=[dict(count=30, label="30 d", step="day", stepmode="backward"),
                         dict(count=90, label="90 d", step="day", stepmode="backward"),
                         dict(count=1, label="YTD", step="year", stepmode="todate"),
                         dict(step="all", label="All")]),
        ),
    )
    return fig


def fig_seasonal(used, cap, lake, s, now) -> go.Figure:
    wse, sd = WSE_COL, SD_COL
    fig = go.Figure()
    labels = _threshold_labels(s)
    groups = []
    years = sorted(used["year"].unique()) if len(used) else []
    cur = max(years) if years else now.year

    for ti, (_, thr) in enumerate(labels):
        d = used[used["n_pix_filtered"] >= thr]
        idx = []
        for yr in years:
            g = d[d["year"] == yr].sort_values("time")
            col = year_color(int(yr), s)
            is_cur = yr == cur
            fig.add_trace(go.Scatter(
                x=g["doy"], y=g[wse], mode="lines+markers", name=str(yr),
                legendgroup=str(yr), opacity=1 if is_cur else 0.75,
                line=dict(color=col, width=3 if is_cur else 1.5),
                marker=dict(size=8 if is_cur else 6, color=col),
                customdata=np.stack([g["time"].dt.strftime("%Y-%m-%d"),
                                     g[sd].round(2), g["n_pix_filtered"]], axis=-1)
                if len(g) else None,
                hovertemplate=(f"<b>{yr}</b> %{{customdata[0]}}<br>"
                               "WSE %{y:.2f} m ± %{customdata[1]} m<br>"
                               "%{customdata[2]} pixels<extra></extra>"),
                visible=ti == 0))
            idx.append(len(fig.data) - 1)
        groups.append(idx)

    always = []
    cap = cap.dropna(subset=["cap_adj_m"]) if len(cap) else cap
    for yr, g in (cap.groupby("year") if len(cap) else []):
        fig.add_trace(go.Scatter(
            x=g["doy"], y=g["cap_adj_m"], mode="markers", name=f"CAP {yr}",
            legendgroup=str(yr),
            marker=dict(size=11, symbol="diamond", color=year_color(int(yr), s),
                        line=dict(color="black", width=1.5)),
            customdata=g["date"].dt.strftime("%Y-%m-%d"),
            hovertemplate="<b>CAP %{customdata}</b><br>Adjusted WSE %{y:.2f} m"
                          "<extra></extra>"))
        always.append(len(fig.data) - 1)

    if len(used):
        lo, hi = used[wse].min(), used[wse].max()
        pad = max(0.5, 0.05 * (hi - lo))
        for rel in lake["releases"]:
            ts = pd.Timestamp(rel)
            fig.add_trace(go.Scatter(
                x=[ts.dayofyear] * 2, y=[lo - pad, hi + pad], mode="lines",
                name=f"{ts.year} release", legendgroup=str(ts.year),
                line=dict(color=year_color(ts.year, s), dash="dash", width=1.2),
                hovertemplate=f"{ts.year} release {ts:%Y-%m-%d}<extra></extra>"))
            always.append(len(fig.data) - 1)
        for c in lake.get("crests", []):
            ts = pd.Timestamp(c)
            always += add_event_lines(
                fig, [ts.dayofyear], [f"{ts.year} crest {ts:%Y-%m-%d}"], lo - pad, hi + pad,
                f"{ts.year} crest", dict(CREST_STYLE, color=year_color(ts.year, s)),
                legendgroup=str(ts.year))

    fig.add_vline(x=now.dayofyear, line_dash="dot", line_width=1,
                  line_color="#9aa5b1", annotation_text="today",
                  annotation_font_size=10)

    months = range(s["season_start_month"], s["season_end_month"] + 1)
    ref = [pd.Timestamp(2025, m, 1) for m in months]
    end = pd.Timestamp(2025, s["season_end_month"], 1) + pd.offsets.MonthEnd(0)
    fig.update_layout(
        **{**BASE_LAYOUT, "margin": dict(l=60, r=20, t=70, b=40)}, height=430,
        updatemenus=_threshold_menu(labels, groups, len(fig.data), always),
        legend=dict(orientation="h", x=0, y=1.13, yanchor="bottom",
                    font=dict(size=11)),
        yaxis=dict(title="WSE (m NAVD88)", tickformat=".1f"),
        xaxis=dict(title=None, tickvals=[r.dayofyear for r in ref],
                   ticktext=[r.strftime("%b") for r in ref],
                   range=[ref[0].dayofyear - 4, end.dayofyear + 4],
                   showgrid=True),
    )
    return fig


PASS_COLORS = ["#1f4e8c", "#d55e00", "#009e73", "#cc79a7", "#56b4e9",
               "#e69f00", "#6a3d9a", "#7f7f7f"]


def fig_by_pass(used, s, now, lake=None) -> go.Figure | None:
    """WSE coloured by SWOT pass number, to show between-pass offsets."""
    d = used.dropna(subset=["pass"])
    if d.empty:
        return None
    fig = go.Figure()
    for i, (p, g) in enumerate(sorted(d.groupby("pass"), key=lambda t: str(t[0]))):
        g = g.sort_values("time")
        col = PASS_COLORS[i % len(PASS_COLORS)]
        # break the line between seasons
        x, y, cd = [], [], []
        prev_year = None
        for r in g.itertuples():
            if prev_year is not None and r.year != prev_year:
                x.append(None); y.append(None); cd.append(["", "", ""])
            x.append(r.time); y.append(getattr(r, WSE_COL))
            cd.append([f"{r.time:%Y-%m-%d %H:%M}", int(r.n_pix_filtered),
                       f"{getattr(r, SD_COL):.2f}" if pd.notna(getattr(r, SD_COL)) else "—"])
            prev_year = r.year
        fig.add_trace(go.Scatter(
            x=x, y=y, mode="lines+markers", name=f"Pass {p} (n={len(g)})",
            line=dict(color=col, width=1.2), marker=dict(size=7, color=col),
            customdata=cd, connectgaps=False,
            hovertemplate=(f"<b>Pass {p}</b> · %{{customdata[0]}} UTC<br>"
                           "WSE %{y:.2f} m ± %{customdata[2]} m<br>"
                           "%{customdata[1]} pixels<extra></extra>")))
    if lake:
        lo, hi = _event_span(d)
        add_event_lines(fig, [pd.Timestamp(r) for r in lake.get("releases", [])],
                        [f"Release {r:%Y-%m-%d}" for r in lake.get("releases", [])], lo, hi,
                        "Release", dict(color="#52606d", dash="dash", width=1))
        add_event_lines(fig, [pd.Timestamp(c) for c in lake.get("crests", [])],
                        [f"Crest {c:%Y-%m-%d}" for c in lake.get("crests", [])], lo, hi,
                        "Crest", CREST_STYLE)
    last_year = int(d["year"].max())
    x0 = pd.Timestamp(last_year, s["season_start_month"], 1)
    x1 = pd.Timestamp(last_year, s["season_end_month"], 1) + pd.offsets.MonthEnd(0)
    fig.update_layout(
        **{**BASE_LAYOUT, "margin": dict(l=60, r=20, t=60, b=30)}, height=420,
        legend=dict(orientation="h", x=0, y=1.12, yanchor="bottom", font=dict(size=11)),
        yaxis=dict(title="WSE (m NAVD88)", tickformat=".1f"),
        xaxis=dict(range=[x0, x1], rangeslider=dict(visible=True, thickness=0.08)),
    )
    return fig


def fig_pixc(pixc: pd.DataFrame, rings, lake, s) -> go.Figure | None:
    if pixc.empty:
        return None
    passes = (pixc.groupby("date").size().rename("n").reset_index()
              .sort_values("date", ascending=False).head(s["max_3d_passes"]))
    ring_x, ring_y = [], []
    for xs, ys in rings:
        ring_x += xs + [None]
        ring_y += ys + [None]

    fig = go.Figure()
    per_pass = 2 if rings else 1
    rng = np.random.default_rng(100)
    for i, r in enumerate(passes.itertuples()):
        g = pixc[pixc["date"] == r.date]
        med = float(np.median(g["elevation_m_NAVD88"]))
        if len(g) > s["max_3d_points"]:
            g = g.iloc[rng.choice(len(g), s["max_3d_points"], replace=False)]
        cd_cols = [c for c in ("pixel_area_sqm", "water_classification",
                               "cross_track_dist_m") if c in g]
        cd = np.stack([g[c].round(1) if g[c].dtype.kind == "f" else g[c]
                       for c in cd_cols], axis=-1) if cd_cols else None
        hover = ("<b>%{z:.2f} m</b><br>%{y:.5f}°, %{x:.5f}°"
                 + "".join(f"<br>{c}: %{{customdata[{k}]}}"
                           for k, c in enumerate(cd_cols)) + "<extra></extra>")
        label = f"{r.date:%Y-%m-%d} — {r.n:,} px"
        fig.add_trace(go.Scatter3d(
            x=g["longitude_WGS84"].round(5), y=g["latitude_WGS84"].round(5),
            z=g["elevation_m_NAVD88"].round(2), mode="markers", name=label,
            customdata=cd, hovertemplate=hover, visible=i == 0,
            marker=dict(size=2.5, opacity=0.7, color=g["elevation_m_NAVD88"],
                        colorscale="Viridis", showscale=True,
                        colorbar=dict(title="WSE (m)", thickness=12, len=0.6))))
        if rings:
            fig.add_trace(go.Scatter3d(
                x=ring_x, y=ring_y, z=[None if v is None else med for v in ring_x],
                mode="lines", name="Lake outline (at median WSE)",
                line=dict(color="black", width=4), hoverinfo="skip",
                visible=i == 0, showlegend=i == 0))

    n = len(fig.data)
    buttons = []
    for i, r in enumerate(passes.itertuples()):
        vis = [False] * n
        for j in range(i * per_pass, i * per_pass + per_pass):
            vis[j] = True
        buttons.append(dict(method="restyle", args=[{"visible": vis}],
                            label=f"{r.date:%Y-%m-%d} — {r.n:,} px"))
    fig.update_layout(
        **{**BASE_LAYOUT, "margin": dict(l=0, r=0, t=50, b=0)},
        height=560, showlegend=False,
        updatemenus=[dict(type="dropdown", x=0, xanchor="left", y=1.08,
                          yanchor="bottom", buttons=buttons, bgcolor="white",
                          bordercolor="#cbd2d9", font=dict(size=11))],
        scene=dict(xaxis_title="Longitude (°)", yaxis_title="Latitude (°)",
                   zaxis_title="WSE (m NAVD88)", aspectmode="auto"),
    )
    return fig


# =============================================================================
# HTML
# =============================================================================

esc = html.escape


def fmt_m(v, signed=False, nd=2):
    if v is None or (isinstance(v, float) and np.isnan(v)):
        return "—"
    return f"{v:+.{nd}f} m" if signed else f"{v:.{nd}f} m"


def fmt_rate(v):
    if v is None:
        return "—"
    return f"{v * 100:+.1f} cm/day"


def arrow(v, eps=0.005):
    if v is None:
        return ""
    return "▲" if v > eps else "▼" if v < -eps else "▶"


class Html:
    def __init__(self):
        self.fig_specs: list[tuple[str, str]] = []

    def plot(self, fig_id: str, fig: go.Figure, cls="plot", autoy=False) -> str:
        js = fig.to_json().replace("</", "<\\/")
        self.fig_specs.append((fig_id, js))
        h = fig.layout.height or 450
        ay = ' data-autoy="1"' if autoy else ""
        return (f'<div class="{cls}" data-fig="{fig_id}"{ay} style="min-height:{h}px">'
                f'<div class="plot-loading">Loading chart…</div></div>')


def table_html(tid, headers, rows, status_filter=None, csv_name="table.csv",
               caption=""):
    """rows: list of (cells, row_attrs). cells: list of (display, sort_value)."""
    head = "".join(f'<th scope="col">{esc(h)}<span class="sort"></span></th>'
                   for h in headers)
    body = []
    for cells, attrs in rows:
        a = " ".join(f'{k}="{esc(str(v))}"' for k, v in attrs.items())
        tds = "".join(f'<td data-v="{esc(str(sv))}">{disp}</td>' for disp, sv in cells)
        body.append(f"<tr {a}>{tds}</tr>")
    sel = ""
    if status_filter:
        opts = "".join(f'<option value="{esc(o)}">{esc(o)}</option>' for o in status_filter)
        sel = f'<select class="status" aria-label="Filter by status"><option value="">All statuses</option>{opts}</select>'
    return f"""
<div class="table-block">
  <div class="table-tools">
    <input class="q" type="search" placeholder="Search…" aria-label="Search table">
    {sel}
    <span class="count"></span>
    <button class="csv" data-name="{esc(csv_name)}">Download CSV</button>
  </div>
  {f'<p class="caption">{caption}</p>' if caption else ''}
  <div class="table-scroll"><table class="dt" id="{tid}">
    <thead><tr>{head}</tr></thead><tbody>{''.join(body)}</tbody>
  </table></div>
</div>"""


def overpass_table(lake_id, passes: pd.DataFrame, predicted: pd.DataFrame,
                   tz: ZoneInfo, name: str) -> str:
    headers = ["Date / time (UTC)", f"Local ({tz.key.split('/')[-1].replace('_', ' ')})",
               "Cycle", "Pass", "Tile", "Pixels", "WSE (m)", "± SD (m)", "Status"]
    rows = []
    for r in predicted.itertuples():
        loc = r.time.tz_localize("UTC").tz_convert(tz)
        rows.append(([
            (f"≈ {r.time:%Y-%m-%d %H:%M}", r.time.isoformat()),
            (f"≈ {loc:%Y-%m-%d %H:%M}", r.time.isoformat()),
            ("—", ""), (esc(str(r.pass_)), r.pass_), (esc(r.tile or "—"), r.tile or ""),
            ("", ""), ("", ""), ("", ""),
            (f'<span class="st st-pred">{STATUS_PRED}</span>', STATUS_PRED),
        ], {"data-status": STATUS_PRED, "class": "pred"}))
    for r in passes.sort_values("time", ascending=False).itertuples():
        loc = r.time.tz_localize("UTC").tz_convert(tz)
        st_cls = {STATUS_USED: "used", STATUS_OUTLIER: "out",
                  STATUS_LOWPIX: "low", STATUS_NOWSE: "nowse"}[r.status]
        w = getattr(r, WSE_COL)
        sd = r.wse_m_NAVD88_filtered_sd
        rows.append(([
            (f"{r.time:%Y-%m-%d %H:%M}", r.time.isoformat()),
            (f"{loc:%Y-%m-%d %H:%M}", r.time.isoformat()),
            (esc(str(r.cycle)) if pd.notna(r.cycle) else "—", r.cycle if pd.notna(r.cycle) else ""),
            (esc(str(r.pass_)) if pd.notna(r.pass_) else "—", r.pass_ if pd.notna(r.pass_) else ""),
            (esc(str(r.tile)) if pd.notna(r.tile) else "—", r.tile if pd.notna(r.tile) else ""),
            (f"{int(r.n_pix_filtered):,}" if pd.notna(r.n_pix_filtered) else "—",
             r.n_pix_filtered if pd.notna(r.n_pix_filtered) else -1),
            (f"{w:.2f}" if pd.notna(w) else "—", w if pd.notna(w) else ""),
            (f"{sd:.2f}" if pd.notna(sd) else "—", sd if pd.notna(sd) else ""),
            (f'<span class="st st-{st_cls}">{esc(r.status)}</span>', r.status),
        ], {"data-status": r.status}))
    statuses = [STATUS_PRED, STATUS_USED, STATUS_OUTLIER, STATUS_LOWPIX, STATUS_NOWSE]
    cap = ("Every SWOT pixel-cloud overpass of the lake inside the seasonal window, "
           "newest first. <b>Used</b> rows are plotted; other rows show why a pass "
           "was left out. <b>Predicted</b> rows project each pass number forward "
           f"on SWOT's ~{SWOT_REPEAT_DAYS:.1f}-day repeat orbit and are approximate.")
    return table_html(f"{lake_id}-passes", headers, rows, statuses,
                      f"{slug(name)}_swot_overpasses.csv", cap)


def age_text(age):
    """Days since the latest valid pass, as plain text (deliberately no colour:
    data age says nothing about how close a lake is to releasing)."""
    if age is None:
        return "No data"
    if age < 1:
        return "< 1 d ago"
    return f"{age:.0f} d ago"


# ---- quality assessment ------------------------------------------------------

def qa_badge(qa, short=False) -> str:
    if not qa:
        return '<span class="qa qa-other">Not assessed</span>'
    c = qa["classification"]
    kind = {QA_LIKELY: "likely", QA_POTENTIAL: "potential"}.get(c, "other")
    txt = {QA_LIKELY: "Likely", QA_POTENTIAL: "Potentially"}.get(c, c) if short else c
    return (f'<span class="qa qa-{kind}" title="SWOT quality assessment: {esc(c)}">'
            f'{esc(txt)}</span>')


def qa_value(qa, key) -> str:
    v = qa.get(key)
    if v is None:
        return "—"
    if key == "area_km2":
        return f"{v:g} km²"
    if key == "gap_days":
        return f"{v:g} d"
    if key in ("pixel_sd_m", "pass_diff_m"):
        star = "*" if key == "pass_diff_m" and qa.get("pass_diff_unfiltered") else ""
        return f"{v:g} m{star}"
    if key == "pixels":
        return f"{int(v):,}"
    return f"{v:.2f}" if isinstance(v, float) else str(v)


def qa_summary(qa) -> str:
    """One line for the lake header: class + which criteria were flagged."""
    if not qa:
        return ""
    if "area_km2" not in qa:   # class set in config only
        return f'<p class="qa-line">SWOT quality assessment: {qa_badge(qa)}</p>'
    flags = qa.get("flags", [])
    if flags:
        what = "; ".join(f"{QA_SHORT[k]} ({qa_value(qa, k)})" for k in flags)
        tail = f"Flagged: {what}."
    else:
        tail = "No quality criteria flagged."
    return (f'<p class="qa-line">SWOT quality assessment: {qa_badge(qa)} '
            f'<span class="qa-why">{esc(tail)}</span></p>')


QA_GUIDANCE = {
    QA_LIKELY: ("No quality criteria were flagged, so the SWOT time series is "
                "internally consistent. That does not mean SWOT samples often enough "
                "to catch every fill–drain cycle at this lake."),
    QA_POTENTIAL: ("One or more criteria were flagged. Treat single observations with "
                   "caution. When passes disagree by more than ~10 m, compare "
                   "observations from the same pass (see <i>By SWOT pass</i>) and "
                   "watch the filling trajectory rather than any one value."),
}


def qa_lake_table(qa, rep) -> str:
    if not qa:
        return '<p class="help">This lake is not in the SWOT quality assessment.</p>'
    if "area_km2" not in qa:
        return (f'<p>{qa_badge(qa)} — class set in the report configuration; '
                'no assessment metrics are available for this lake.</p>')
    rows = []
    for k, label, _, test, thr in QA_CRITERIA:
        v = qa.get(k)
        flagged = k in qa.get("flags", [])
        rows.append(f'<tr{" class=flag" if flagged else ""}><td>{esc(label)}</td>'
                    f'<td class="num">{qa_value(qa, k)}</td><td>{esc(thr)}</td>'
                    f'<td>{"Flagged" if flagged else ("—" if v is None else "OK")}</td></tr>')
    star = ('<p class="help">* Pass difference computed from all (unfiltered) pixels '
            'because fewer than two SWOT passes had filtered observations.</p>'
            if qa.get("pass_diff_unfiltered") else "")
    guide = QA_GUIDANCE.get(qa["classification"], "")
    return f"""
<p>{qa_badge(qa)} {guide}</p>
<table class="grid qa-table"><thead><tr><th>Criterion</th><th>Value</th>
<th>Flag if</th><th>Status</th></tr></thead><tbody>{''.join(rows)}</tbody></table>
{star}
<p class="help">Source: {esc(rep['assessment_source'])}. The assessment rates the whole
time series and uses only remotely sensed inputs (no in-situ data). It is a first
indicator of usability, not a check of any single observation.</p>"""


def qa_overview_table(lakes_out) -> str:
    headers = ["Lake", "Assessment", "Area (km²)", "Pixel SD (m)", "Mean gap (d)",
               "Pass diff. (m)", "Avg. pixels", "TRI", "Flagged criteria"]
    rows = []
    for lo in lakes_out:
        qa = lo["qa"] or {}
        def cell(k, nd=None):
            v = qa.get(k)
            if v is None:
                return ("—", "")
            disp = qa_value(qa, k).replace(" km²", "").replace(" m", "").replace(" d", "")
            if k in qa.get("flags", []):
                disp = f'<b class="flagged" title="exceeds threshold">{disp}</b>'
            return (disp, v)
        flags = ", ".join(QA_SHORT[k] for k in qa.get("flags", [])) or "none"
        rows.append(([
            (f'<a href="#{lo["id"]}">{esc(lo["lake"]["display_name"])}</a>',
             lo["lake"]["display_name"]),
            (qa_badge(qa or None), QA_RANK.get(qa.get("classification"), 9)),
            cell("area_km2"), cell("pixel_sd_m"), cell("gap_days"), cell("pass_diff_m"),
            cell("pixels"), cell("tri"),
            (esc(flags) if "area_km2" in qa else "—", flags),
        ], {}))
    cap = ("Bold values exceed the threshold for that criterion. "
           "* pass difference from unfiltered pixels.")
    return table_html("qa-all", headers, rows, None, "gdl_swot_quality_assessment.csv", cap)


# ---- KPIs and overview -------------------------------------------------------

def kpi(label, value, sub="", cls=""):
    return (f'<div class="kpi {cls}"><div class="kpi-label">{label}</div>'
            f'<div class="kpi-value">{value}</div>'
            f'<div class="kpi-sub">{sub}</div></div>')


def kpi_block(sm: dict, s: dict, tz: ZoneInfo) -> str:
    if "latest" not in sm:
        nxt = sm.get("next_pass")
        return ('<div class="kpis">' + kpi("Latest WSE", "—", "No valid observations") +
                (kpi("Next overpass", f"≈ {nxt['time']:%b %d}", f"pass {nxt['pass']}")
                 if nxt is not None else "") + "</div>")
    L = sm["latest"]
    age = age_text(sm["age_days"])
    out = [kpi("Latest WSE",
               fmt_m(L[WSE_COL]),
               f'<b>{L["time"]:%b %d, %Y %H:%M} UTC</b> ({age}) · '
               f'±{L[SD_COL]:.2f} m · {int(L["n_pix_filtered"]):,} px',
               "kpi-main")]
    if "delta_prev" in sm:
        dlt = sm["delta_prev"]
        out.append(kpi("Since previous pass",
                       f'{arrow(dlt)} {fmt_m(dlt, True)}',
                       f'over {sm["delta_prev_days"]:.0f} d '
                       f'({fmt_rate(dlt / max(sm["delta_prev_days"], 0.01))})'))
    if "trend_m_per_day" in sm:
        tr = sm["trend_m_per_day"]
        out.append(kpi(f'{s["trend_window_days"]}-day trend',
                       f'{arrow(tr / 10)} {fmt_rate(tr)}',
                       f'linear fit, {sm["trend_n"]} obs'))
    if "vs_prior_year" in sm:
        v, y = sm["vs_prior_year"]
        out.append(kpi(f"vs {y}, same date", f"{arrow(v)} {fmt_m(v, True)}",
                       "interpolated"))
    if "season_peak" in sm:
        v, t = sm["season_peak"]
        out.append(kpi(f'{L["year"]} season peak', fmt_m(v), f"{t:%b %d}"))
    if "vs_release" in sm:
        dv, peak, rel = sm["vs_release"]
        out.append(kpi("vs last pre-release peak", f"{arrow(dv)} {fmt_m(dv, True)}",
                       f"peak {peak:.2f} m before {rel:%Y-%m-%d} release"))
    nxt = sm.get("next_pass")
    if nxt is not None:
        loc = nxt["time"].tz_localize("UTC").tz_convert(tz)
        out.append(kpi("Next overpass (predicted)", f"≈ {nxt['time']:%b %d}",
                       f"pass {nxt['pass']} · ≈{loc:%H:%M} local"))
    return '<div class="kpis">' + "".join(out) + "</div>"


def overview_table(lakes_out, s) -> str:
    headers = ["Lake", "Assessment", "Latest obs (UTC)", "Age", "WSE (m)", "Δ prev pass",
               f"{s['trend_window_days']}-d trend", "vs prior yr",
               "vs pre-release peak", "Next overpass"]
    rows = []
    for lo in lakes_out:
        sm = lo["summary"]
        L = sm.get("latest")
        d = sm.get("delta_prev")
        tr = sm.get("trend_m_per_day")
        vy = sm.get("vs_prior_year", (None, None))[0]
        vr = sm.get("vs_release", (None,))[0]
        nxt = sm.get("next_pass")
        qa = lo["qa"]
        rows.append(([
            (f'<a href="#{lo["id"]}">{esc(lo["lake"]["display_name"])}</a>',
             lo["lake"]["display_name"]),
            (qa_badge(qa, short=True), QA_RANK.get((qa or {}).get("classification"), 9)),
            (f"{L['time']:%Y-%m-%d %H:%M}" if L is not None else "—",
             L["time"].isoformat() if L is not None else ""),
            (age_text(sm.get("age_days")), sm.get("age_days", 9e9)),
            (fmt_m(L[WSE_COL]).replace(" m", "") if L is not None else "—",
             L[WSE_COL] if L is not None else ""),
            (f"{arrow(d)} {fmt_m(d, True)}" if d is not None else "—", d if d is not None else ""),
            (f"{arrow(tr / 10 if tr else None)} {fmt_rate(tr)}" if tr is not None else "—",
             tr if tr is not None else ""),
            (f"{arrow(vy)} {fmt_m(vy, True)}" if vy is not None else "—", vy if vy is not None else ""),
            (f"{arrow(vr)} {fmt_m(vr, True)}" if vr is not None else "—", vr if vr is not None else ""),
            (f"≈ {nxt['time']:%b %d} · p{nxt['pass']}" if nxt is not None else "—",
             nxt["time"].isoformat() if nxt is not None else ""),
        ], {}))
    return table_html("overview-tbl", headers, rows, None, "gdl_overview.csv")


# ---- maps --------------------------------------------------------------------

def map_payload(lakes_out, rep) -> dict:
    lakes = []
    for lo in lakes_out:
        sm, qa = lo["summary"], lo["qa"] or {}
        L = sm.get("latest")
        rings = []
        for xs, ys in lo.get("rings", []):
            step = max(1, len(xs) // 1500)          # keep the page small
            rings.append([[round(y, 5), round(x, 5)] for x, y in
                          list(zip(xs, ys))[::step]])
        lat, lon = qa.get("lat"), qa.get("lon")
        if lat is None and rings:
            pts = [p for r in rings for p in r]
            lat = sum(p[0] for p in pts) / len(pts)
            lon = sum(p[1] for p in pts) / len(pts)
        latest = (f'{fmt_m(L[WSE_COL])} on {L["time"]:%b %d, %Y} '
                  f'({age_text(sm["age_days"])})') if L is not None else "No valid observations"
        popup = (f'<b>{esc(lo["lake"]["display_name"])}</b><br>{qa_badge(qa or None)}<br>'
                 f'Latest WSE: {latest}<br><a href="#{lo["id"]}">Go to lake &raquo;</a>')
        lakes.append({"id": lo["id"], "name": lo["lake"]["display_name"],
                      "cls": qa.get("classification", ""), "lat": lat, "lon": lon,
                      "rings": rings, "popup": popup})
    base = rep.get("basemap", "Imagery")
    return {"lakes": lakes, "basemap": base if base in ESRI_TILES else "Imagery",
            "tiles": ESRI_TILES, "labels": ESRI_LABELS}


CSS = """
:root{--ink:#222;--muted:#555f6b;--line:#c8d2de;--soft:#e8eef5;--navy:#0a3a7a;
--navy2:#1f4e8c;--side:230px}
*{box-sizing:border-box}
html{scroll-behavior:smooth;scroll-padding-top:8px}
body{margin:0;font:13px/1.45 Arial,Helvetica,sans-serif;color:var(--ink);background:#dfe6ee;
font-variant-numeric:tabular-nums}
a{color:#0b4ea2;text-decoration:none}a:hover{text-decoration:underline}
#page{max-width:1480px;margin:0 auto;background:#fff;box-shadow:0 0 6px rgba(0,0,0,.18)}
#banner{background:linear-gradient(90deg,#0b2a5b 0%,#14468f 55%,#4f86c9 100%);color:#fff;
padding:16px 22px 14px}
#banner .t1{font-size:12px;letter-spacing:.06em;text-transform:uppercase;opacity:.85}
#banner .t2{font-size:24px;font-weight:bold;margin-top:2px}
#banner .t3{font-size:12px;opacity:.85;margin-top:3px}
#topnav{background:linear-gradient(#3b5f8e,#22416b);display:flex;flex-wrap:wrap;
border-top:1px solid #6d8fbb;border-bottom:3px solid #0b2a5b}
#topnav a{color:#fff;font-weight:bold;font-size:12px;text-transform:uppercase;
padding:8px 16px;border-right:1px solid #5d7ea9}
#topnav a:hover{background:#4f77ab;text-decoration:none}
#topnav a.ext{margin-left:auto;border-left:1px solid #5d7ea9;border-right:0}
#crumb{font-size:11.5px;color:var(--muted);padding:6px 22px;border-bottom:1px solid var(--soft)}
#wrap{display:flex;align-items:flex-start}
#side{width:var(--side);flex:none;padding:14px 0 14px 14px;position:sticky;top:0;
max-height:100vh;overflow-y:auto}
.box{border:1px solid var(--line);margin-bottom:12px;background:#f6f9fc}
.box h4{margin:0;background:var(--navy2);color:#fff;font-size:12px;padding:6px 9px}
.box a{display:flex;justify-content:space-between;gap:6px;padding:5px 9px;
border-top:1px solid #dde5ee;color:#0b3d82;font-size:12.5px}
.box a:first-of-type{border-top:0}
.box a:hover{background:#e3ecf6;text-decoration:none}
.box a.active{background:#d3e1f1;font-weight:bold}
.box a small{color:var(--muted);font-weight:normal;white-space:nowrap}
.box .txt{padding:6px 9px;font-size:11.5px;color:var(--muted)}
main{flex:1;min-width:0;padding:14px 22px 40px}
section{margin-bottom:30px}
h2{font-size:19px;color:var(--navy);margin:0 0 8px;padding-bottom:4px;
border-bottom:1px solid var(--line);display:flex;align-items:baseline;gap:10px;flex-wrap:wrap}
h2 a.totop{margin-left:auto;font-size:11.5px;font-weight:normal}
h3{font-size:14px;color:var(--navy);margin:18px 0 4px}
.lake-notes{color:var(--muted);margin:0 0 6px;font-size:12.5px}
.lake-top{display:grid;grid-template-columns:minmax(0,1fr) 360px;gap:16px;align-items:start}
.qa-line{margin:2px 0 6px}.qa-why{color:var(--muted)}
.qa{display:inline-block;padding:1px 8px;border-radius:3px;font-size:11.5px;font-weight:bold;
white-space:nowrap;vertical-align:1px}
.qa-likely{background:var(--navy2);color:#fff;border:1px solid var(--navy2)}
.qa-potential{background:#fff;color:var(--navy2);border:1px solid var(--navy2)}
.qa-other{background:#eee;color:#444;border:1px solid #bbb}
.kpis{display:grid;grid-template-columns:repeat(auto-fill,minmax(170px,1fr));gap:8px;margin:6px 0}
.kpi{border:1px solid var(--line);padding:7px 10px;background:#fafbfd}
.kpi-main{grid-column:span 2;background:#f1f5fa}
.kpi-label{font-size:11px;color:var(--muted);text-transform:uppercase;letter-spacing:.03em}
.kpi-value{font-size:19px;font-weight:bold;margin:1px 0;white-space:nowrap}
.kpi-main .kpi-value{font-size:25px}
.kpi-sub{font-size:11.5px;color:var(--muted)}.kpi-sub b{color:var(--ink);font-weight:normal}
.map{border:1px solid var(--line);background:#e8eef5}
.map-lg{height:430px;margin:8px 0 14px}.map-sm{height:270px}
.map-off{display:flex;height:100%;align-items:center;justify-content:center;color:var(--muted);
padding:10px;text-align:center}
.map .leaflet-control-attribution{font-size:10px;line-height:1.3}
.legend{background:rgba(255,255,255,.92);padding:6px 8px;font-size:11.5px;line-height:1.6;
border:1px solid #999}
.legend i{display:inline-block;width:11px;height:11px;border-radius:50%;margin-right:5px;
vertical-align:-1px}
.warn{background:#fff8e1;border:1px solid #e5c76b;padding:6px 10px;margin:8px 0;font-size:12.5px}
.warn ul{margin:0;padding-left:18px}
.help{color:var(--muted);font-size:12.5px;margin:2px 0 6px}
.tabs{display:flex;flex-wrap:wrap;gap:2px;border-bottom:1px solid var(--line);margin-top:12px}
.tab{font:inherit;font-size:12.5px;font-weight:bold;color:#0b3d82;background:var(--soft);
border:1px solid var(--line);border-bottom:0;padding:6px 13px;cursor:pointer;margin-bottom:-1px}
.tab:hover{background:#dbe5f0}
.tab.active{background:#fff;color:var(--ink);border-top:3px solid var(--navy2);padding-top:4px}
.tabpanel{border:1px solid var(--line);border-top:0;padding:10px 12px}
.plot{position:relative;width:100%}
.plot-loading{position:absolute;inset:0;display:flex;align-items:center;justify-content:center;
color:var(--muted)}
.table-block{margin-top:6px}
.table-tools{display:flex;gap:8px;align-items:center;flex-wrap:wrap;margin-bottom:6px}
.table-tools input,.table-tools select{padding:4px 7px;border:1px solid #aab6c4;font:inherit}
.table-tools .count{color:var(--muted);font-size:12px}
.table-tools button{margin-left:auto;padding:4px 10px;border:1px solid #aab6c4;
background:var(--soft);cursor:pointer;font:inherit;font-size:12px}
.table-tools button:hover{background:#dbe5f0}
.caption{color:var(--muted);font-size:12px;margin:0 0 6px}
.table-scroll{max-height:440px;overflow:auto;border:1px solid var(--line)}
table.dt,table.grid{border-collapse:collapse;width:100%;font-size:12.5px}
table.dt th,table.grid th{position:sticky;top:0;background:#dde6f0;text-align:left;
padding:6px 9px;border:1px solid var(--line);white-space:nowrap;color:#1b2b40}
table.dt th{cursor:pointer;user-select:none}
table.dt td,table.grid td{padding:4px 9px;border:1px solid #dfe5ec;white-space:nowrap}
table.dt tbody tr:nth-child(even) td,table.grid tbody tr:nth-child(even) td{background:#f7f9fb}
table.dt tr:hover td{background:#eaf1f9}
table.dt tr.pred td{color:#4a5563;font-style:italic;background:#f4f4f4}
table.grid{width:auto;margin:6px 0}table.grid td.num{text-align:right}
table.qa-table tr.flag td{font-weight:bold}
b.flagged{border-bottom:1px dotted #333}
th .sort::after{content:" ↕";color:#9fb0c4}th.asc .sort::after{content:" ↑";color:var(--ink)}
th.desc .sort::after{content:" ↓";color:var(--ink)}
.st{display:inline-block;padding:0 6px;font-size:11.5px;font-style:normal;border:1px solid #bbb;
background:#fff}
.st-used{border-color:#7a9cc6;color:#1b3f6e}.st-out,.st-low,.st-nowse{color:#555}
.st-pred{color:#555;border-style:dashed}
.notes li{margin-bottom:5px}
code{background:#eef2f6;padding:1px 4px;font-size:12px}
footer{border-top:3px solid #0b2a5b;background:#eef2f6;padding:12px 22px;font-size:11.5px;
color:var(--muted)}
footer p{margin:3px 0}
@media (max-width:1000px){.lake-top{grid-template-columns:1fr}}
@media (max-width:820px){#wrap{display:block}#side{position:static;width:auto;max-height:none;
padding:10px}main{padding:10px}.kpi-main{grid-column:span 1}}
@media print{#side,#topnav,.tabs,.map{display:none}.tabpanel[hidden]{display:block}}
"""

JS = r"""
(function(){
const cfg={responsive:true,displaylogo:false,
  modeBarButtonsToRemove:['lasso2d','select2d']};
function render(el){
  if(el.dataset.done)return; el.dataset.done=1;
  const spec=JSON.parse(document.getElementById('fig-'+el.dataset.fig).textContent);
  el.innerHTML='';
  Plotly.newPlot(el,spec.data,spec.layout,
    Object.assign({},cfg,{toImageButtonOptions:{format:'png',scale:2,filename:el.dataset.fig}}))
  .then(()=>{
    if(!el.dataset.autoy)return;
    autoY(el);
    el.on('plotly_relayout',ev=>{
      if(Object.keys(ev).some(k=>k.startsWith('xaxis.range')||k==='xaxis.autorange'))autoY(el);});
    el.on('plotly_restyle',()=>autoY(el));
  });
}
// Fit the y axis to the points inside the visible date range.
function autoY(el){
  const xr=el._fullLayout.xaxis.range;if(!xr)return;
  const x0=new Date(xr[0]).getTime(),x1=new Date(xr[1]).getTime();
  let lo=Infinity,hi=-Infinity;
  el._fullData.forEach(t=>{
    if(t.visible!==true||t.fill==='toself'||t.meta==='event'||!t.x)return;
    for(let i=0;i<t.x.length;i++){
      if(t.x[i]==null)continue;
      const tx=new Date(t.x[i]).getTime(),y=t.y[i];
      if(tx>=x0&&tx<=x1&&y!=null&&isFinite(y)){
        const e=(t.error_y&&t.error_y.visible&&t.error_y.array)?(+t.error_y.array[i]||0):0;
        if(y-e<lo)lo=y-e;if(y+e>hi)hi=y+e;}
    }});
  if(!isFinite(lo))return;
  const pad=Math.max(0.5,(hi-lo)*0.1);
  Plotly.relayout(el,{'yaxis.range':[lo-pad,hi+pad]});
}
const io=new IntersectionObserver(es=>es.forEach(e=>{
  if(e.isIntersecting){render(e.target);io.unobserve(e.target);}}),{rootMargin:'300px'});
document.querySelectorAll('.plot[data-fig]').forEach(el=>io.observe(el));

// Tabs within each lake section.
document.querySelectorAll('.tabs').forEach(tl=>tl.addEventListener('click',e=>{
  const b=e.target.closest('.tab');if(!b)return;
  tl.querySelectorAll('.tab').forEach(x=>{const on=x===b;x.classList.toggle('active',on);
    x.setAttribute('aria-selected',on);document.getElementById(x.dataset.tab).hidden=!on;});
  document.getElementById(b.dataset.tab).querySelectorAll('.plot[data-done]')
    .forEach(el=>Plotly.Plots.resize(el));
}));

document.querySelectorAll('.table-block').forEach(block=>{
  const t=block.querySelector('table'),tb=t.tBodies[0],rows=[...tb.rows];
  const q=block.querySelector('input.q'),sel=block.querySelector('select.status'),
        cnt=block.querySelector('.count');
  function apply(){
    const s=(q.value||'').toLowerCase(),st=sel?sel.value:'';let n=0;
    rows.forEach(r=>{const ok=(!s||r.textContent.toLowerCase().includes(s))&&
      (!st||r.dataset.status===st);r.hidden=!ok;if(ok)n++;});
    cnt.textContent=n+' of '+rows.length+' rows';
  }
  q.addEventListener('input',apply); if(sel)sel.addEventListener('change',apply); apply();
  [...t.tHead.rows[0].cells].forEach((th,i)=>th.addEventListener('click',()=>{
    const asc=!th.classList.contains('asc');
    [...t.tHead.rows[0].cells].forEach(h=>h.classList.remove('asc','desc'));
    th.classList.add(asc?'asc':'desc');
    const val=r=>{const v=r.cells[i].dataset.v;const f=parseFloat(v);
      return v!==''&&!isNaN(f)&&/^-?[\d.eE+-]+$/.test(v)?f:v;};
    const sorted=[...tb.rows].sort((a,b)=>{const x=val(a),y=val(b);
      if(x===''&&y!=='')return 1;if(y===''&&x!=='')return -1;
      return (x<y?-1:x>y?1:0)*(asc?1:-1);});
    sorted.forEach(r=>tb.appendChild(r));
  }));
  block.querySelector('button.csv').addEventListener('click',e=>{
    const quote=s=>'"'+String(s).replace(/"/g,'""')+'"';
    const lines=[[...t.tHead.rows[0].cells].map(c=>quote(c.textContent.trim()))];
    [...tb.rows].filter(r=>!r.hidden).forEach(r=>lines.push(
      [...r.cells].map(c=>quote(c.textContent.trim()))));
    const blob=new Blob([lines.map(l=>l.join(',')).join('\n')],{type:'text/csv'});
    const a=document.createElement('a');a.href=URL.createObjectURL(blob);
    a.download=e.target.dataset.name;a.click();URL.revokeObjectURL(a.href);
  });
});

// Esri basemaps (Leaflet).
const MD=JSON.parse(document.getElementById('map-data').textContent);
const byId=Object.fromEntries(MD.lakes.map(k=>[k.id,k]));
function markerStyle(c){return c==='Likely successful'
  ?{radius:7,color:'#fff',weight:2,fillColor:'#1f4e8c',fillOpacity:1}
  :{radius:6,color:'#1f4e8c',weight:3,fillColor:'#fff',fillOpacity:1};}
function makeMap(el,lakes,opt){
  if(typeof L==='undefined'){el.innerHTML='<div class="map-off">Map unavailable: the map '+
    'library could not be loaded (no internet connection?).</div>';return;}
  const bl={};for(const [k,v] of Object.entries(MD.tiles))
    bl[k]=L.tileLayer(v[0],{attribution:v[1],maxZoom:18});
  const labels=L.tileLayer(MD.labels,{maxZoom:18});
  const m=L.map(el,{scrollWheelZoom:false,layers:[bl[MD.basemap]]});
  if(MD.basemap==='Imagery')labels.addTo(m);
  L.control.layers(bl,{'Place names':labels}).addTo(m);
  L.control.scale().addTo(m);
  const b=L.latLngBounds([]);
  lakes.forEach(k=>{
    k.rings.forEach(r=>{const p=L.polygon(r,{color:'#3fd0ff',weight:2,fillOpacity:.12})
      .addTo(m);if(opt.single)b.extend(p.getBounds());});
    if(k.lat==null)return;
    if(!opt.single||!k.rings.length){
      const mk=L.circleMarker([k.lat,k.lon],markerStyle(k.cls)).addTo(m);
      if(!opt.single){mk.bindPopup(k.popup);mk.bindTooltip(k.name);}
    }
    if(!opt.single||!k.rings.length)b.extend([k.lat,k.lon]);
  });
  if(b.isValid())m.fitBounds(b.pad(opt.single?0.25:0.08),{maxZoom:opt.single?13:9});
  else m.setView([61,-147],5);
  if(!opt.single){
    const lg=L.control({position:'bottomleft'});
    lg.onAdd=()=>{const d=L.DomUtil.create('div','legend');d.innerHTML=
      '<i style="background:#1f4e8c;border:2px solid #fff;outline:1px solid #1f4e8c"></i>Likely successful<br>'+
      '<i style="background:#fff;border:3px solid #1f4e8c"></i>Potentially successful';return d;};
    lg.addTo(m);
  }
}
const mo=new IntersectionObserver(es=>es.forEach(e=>{
  if(!e.isIntersecting)return;mo.unobserve(e.target);const el=e.target;
  if(el.dataset.lake)makeMap(el,[byId[el.dataset.lake]],{single:true});
  else makeMap(el,MD.lakes,{single:false});
}),{rootMargin:'300px'});
document.querySelectorAll('.map').forEach(el=>mo.observe(el));

// Highlight the lake currently in view in the side menu.
const links=[...document.querySelectorAll('#side a[data-sec]')];
const secObs=new IntersectionObserver(es=>es.forEach(e=>{
  if(e.isIntersecting){links.forEach(a=>a.classList.toggle('active',a.dataset.sec===e.target.id));}
}),{rootMargin:'-30% 0px -60% 0px'});
document.querySelectorAll('main > section[id]').forEach(s=>secObs.observe(s));
})();
"""


def tabs_html(lid: str, tabs: list[tuple[str, str, str]]) -> str:
    """tabs: (key, label, panel html). First tab is open."""
    btns, panels = [], []
    for i, (key, label, body) in enumerate(tabs):
        pid = f"{lid}-{key}"
        btns.append(f'<button class="tab{" active" if i == 0 else ""}" role="tab" '
                    f'aria-selected="{"true" if i == 0 else "false"}" data-tab="{pid}">'
                    f'{esc(label)}</button>')
        panels.append(f'<div class="tabpanel" id="{pid}" role="tabpanel"'
                      f'{"" if i == 0 else " hidden"}>{body}</div>')
    return f'<div class="tabs" role="tablist">{"".join(btns)}</div>{"".join(panels)}'


def build_html(cfg, lakes_out, meta, plotly_js_mode) -> str:
    s = cfg["settings"]
    rep = cfg["report"]
    tz = meta["tz"]
    H = meta["html"]

    # ---- side menu, grouped by assessment class
    groups: dict[str, list] = {}
    for lo in lakes_out:
        groups.setdefault((lo["qa"] or {}).get("classification") or "Other lakes", []).append(lo)
    side = []
    for cls in sorted(groups, key=lambda c: QA_RANK.get(c, 9)):
        items = []
        for lo in groups[cls]:
            sm = lo["summary"]
            sub = f'{sm["latest"]["time"]:%b %d}' if "latest" in sm else "no data"
            items.append(f'<a href="#{lo["id"]}" data-sec="{lo["id"]}">'
                         f'<span>{esc(lo["lake"]["display_name"])}</span><small>{sub}</small></a>')
        side.append(f'<div class="box"><h4>{esc(cls)}</h4>{"".join(items)}</div>')
    side.append(
        '<div class="box"><h4>Related links</h4>'
        f'<a href="{GDL_MAIN_URL}" target="_blank" rel="noopener">NWS APRFC glacier-dammed lakes</a>'
        f'<a href="{DATA_RELEASE_URL}" target="_blank" rel="noopener">USGS data release</a>'
        '<a href="https://podaac.jpl.nasa.gov/SWOT" target="_blank" rel="noopener">SWOT at PO.DAAC</a></div>')
    side.append(f'<div class="box"><h4>Report</h4><div class="txt">Generated '
                f'{esc(meta["gen_txt"])}<br>Data through '
                f'{meta["data_through"]:%Y-%m-%d %H:%M} UTC</div></div>')

    # ---- lake sections
    sections = []
    for lo in lakes_out:
        L, lid, qa = lo["lake"], lo["id"], lo["qa"]
        warn = ""
        if lo["warnings"]:
            warn = ('<div class="warn"><ul>' +
                    "".join(f"<li>{esc(w)}</li>" for w in lo["warnings"]) + "</ul></div>")
        gdl = L.get("gdl_url") or rep.get("gdl_url", "")
        note_bits = [esc(L["notes"])] if L["notes"] else []
        if gdl:
            note_bits.append(f'<a href="{esc(gdl)}" target="_blank" rel="noopener">'
                             f'{esc(rep.get("gdl_link_text") or "Glacial Dammed Lake Data")}</a>')
        notes_html = (f'<p class="lake-notes">{" · ".join(note_bits)}</p>' if note_bits else "")

        ts = H.plot(lid + "-ts", lo["fig_ts"], autoy=True) if lo["fig_ts"] else '<p class="help">No valid observations.</p>'
        se = H.plot(lid + "-se", lo["fig_se"]) if lo["fig_se"] else ""
        wse_tab = f"""
  <p class="help">Filtered median WSE; bars and shading are ±1 SD. Diamonds are Civil Air
  Patrol (CAP) estimates shifted by that year's mean CAP–SWOT offset. Use the buttons above
  the chart to zoom (30 d, 90 d, YTD, All) and the dropdown to require more pixels per pass.
  Dashed lines are documented releases and dash-dot lines are outburst crests; hover for dates.
  A dotted line joins passes more than {s['max_gap_days']} days apart.
  Read the filling trajectory across several passes rather than any single value.</p>
  {ts}
  <h3>Seasonal comparison</h3>
  <p class="help">Same observations by day of year; current year is drawn heavier.
  Dashed lines mark documented releases; dash-dot lines mark crests.</p>
  {se}"""
        tabs = [("wse", "Water-surface elevation", wse_tab)]
        if lo["fig_pass"] is not None:
            pd_txt = ""
            if qa and qa.get("pass_diff_m") is not None:
                pd_txt = (f" In the quality assessment, the mean difference between passes "
                          f"at this lake was {qa_value(qa, 'pass_diff_m')}.")
            tabs.append(("pass", "By SWOT pass", f"""
  <p class="help">The same valid observations, coloured by SWOT pass number. Different
  passes view the lake from different geometries and can disagree by metres (an
  offset, not a real change in lake level).{esc(pd_txt)} Comparing points from one pass is
  the most reliable way to judge filling. Opens on the latest season.</p>
  {H.plot(lid + "-bp", lo["fig_pass"], autoy=True)}"""))
        tabs.append(("op", "Overpasses", lo["table_passes"]))
        tabs.append(("qa", "Quality assessment", qa_lake_table(qa, rep)))
        if lo["fig_pixc"] is not None:
            tabs.append(("3d", "Pixel cloud (3-D)", f"""
  <p class="help">{meta['pixc_label']} (the same pixels that make up the WSE) for the latest
  {lo["n_pixc_passes"]} SWOT passes; pick
  the date in the dropdown. The black line is the lake outline drawn at that pass's
  median WSE. At most {s["max_3d_points"]:,} points are shown per pass (random sample).</p>
  {H.plot(lid + "-3d", lo["fig_pixc"])}"""))

        sections.append(f"""
<section id="{lid}">
  <h2>{esc(L['display_name'])} {qa_badge(qa)}<a class="totop" href="#top">&uarr; top</a></h2>
  <div class="lake-top">
    <div>
      {notes_html}
      {qa_summary(qa)}
      {kpi_block(lo["summary"], s, tz)}
      {warn}
    </div>
    <div class="map map-sm" data-lake="{lid}" aria-label="Map of {esc(L['display_name'])}"></div>
  </div>
  {tabs_html(lid, tabs)}
</section>""")

    # ---- quality assessment section
    included = ", ".join(s["include_classes"]) if s["include_classes"] else "all classes"
    crit_rows = "".join(f"<tr><td>{esc(lab)}</td><td>{esc(thr)}</td></tr>"
                        for _, lab, _, _, thr in QA_CRITERIA)
    qa_sec = f"""
<section id="qa"><h2>SWOT quality assessment</h2>
<p>SWOT was evaluated at 22 glacier-dammed lakes monitored by the National Weather Service.
Each lake's whole time series was classed as <i>likely successful</i>, <i>potentially
successful</i>, <i>failure</i>, or <i>no data</i> ({esc(rep['assessment_source'])}).
This report shows only glacier-dammed lakes classed as: <b>{esc(included)}</b>.</p>
{qa_overview_table(lakes_out)}
<div style="display:flex;gap:28px;flex-wrap:wrap;align-items:flex-start">
<div><h3>Criteria</h3>
<table class="grid"><thead><tr><th>Criterion</th><th>Flag if</th></tr></thead>
<tbody>{crit_rows}</tbody></table></div>
<div style="flex:1;min-width:300px"><h3>Classes</h3><ul class="notes">
<li>{qa_badge({"classification": QA_LIKELY})} No criteria flagged.</li>
<li>{qa_badge({"classification": QA_POTENTIAL})} One or more criteria flagged. Every
potentially successful lake had mean differences between SWOT passes of 15–62 m, and
several had mean pixel standard deviations above 10 m.</li>
<li><b>Failure</b>: mean gap between valid observations over 15 days. <b>No data</b>: no
valid filtered observations.{" Neither is shown in this report." if s["include_classes"] else ""}</li>
</ul>
<h3>Using SWOT for glacier-dammed lakes</h3><ul class="notes">
<li>Combining ascending and descending passes, SWOT sees most lakes every ~7–10 days.
That is too slow for real-time outburst warning, but it can track pre-release filling.
It complements Civil Air Patrol flights rather than replacing them.</li>
<li>The assessment checks that a time series is internally consistent. It does not check
that SWOT samples often enough to resolve a lake's fill–drain cycle.</li>
<li>At five gaged validation lakes, filtered SWOT WSE agreed well with gages
(R² 0.83–1; unbiased RMSE 0.15–2.63 m). Absolute biases of −1.2 to +3.3 m remained, so
changes in level are more reliable than absolute elevations.</li>
<li>Glacier-dammed lakes return many dark-water and low-coherence pixels. Only open-water
pixels are used, so pixel counts are often low.</li>
</ul></div></div>
</section>"""

    years_txt = ", ".join(map(str, s["years"])) if s["years"] else "all years in the data"
    mon = lambda m: dt.date(2000, m, 1).strftime("%B")
    notes = f"""
<section id="notes"><h2>Processing notes</h2><ul class="notes">
<li><b>Data source:</b> {meta['source_html']}. Data release:
 <a href="{DATA_RELEASE_URL}">{DATA_RELEASE_URL}</a>.</li>
<li><b>Lakes:</b> {("glacier-dammed lakes in the SWOT quality assessment whose class is "
 + esc(included) + ". Gaged validation lakes and other lakes are not shown.") if s["include_classes"]
 else "all lakes in the data (no quality-assessment filter)."}</li>
<li><b>Seasonal window:</b> {mon(s['season_start_month'])}–{mon(s['season_end_month'])};
 years: {years_txt}. Passes outside the window are not shown anywhere.</li>
<li><b>Valid pass:</b> at least {s['min_pixels']} pixels and a non-missing
 WSE. WSE plotted and summarized is {meta['wse_label']};
 error bars and shading are {meta['sd_label']}.</li>
<li><b>Outliers:</b> Hampel filter on that WSE, per lake,
 window ±{s['hampel_k']} passes, threshold {s['hampel_t0']} × MAD (applied only when a lake
 has ≥ {s['min_obs_for_hampel']} valid passes). Flagged passes appear in the overpass table.</li>
<li><b>Lines:</b> passes more than {s['max_gap_days']} days apart are joined by a faint
 dotted line, to show that the level between them is unknown.</li>
<li><b>Crests:</b> outburst crest dates set in the report configuration
 (<code>crests</code> for each lake), drawn as dash-dot lines.</li>
<li><b>Trend:</b> least-squares slope of WSE over the {s['trend_window_days']} days ending at
 the latest valid pass (needs ≥ 3 passes). <b>vs prior year</b> interpolates last year's series
 to the same day of year. <b>vs pre-release peak</b> compares with the highest valid WSE in the
 {s['pre_release_window_days']} days before the most recent documented release.</li>
<li><b>Age:</b> days since the latest valid pass. It shows how current the data are, not
 how close a lake is to releasing.</li>
<li><b>Predicted overpasses:</b> each pass number's last observed time plus multiples of
 {SWOT_REPEAT_DAYS} days. Approximate; data latency after an overpass is not included.</li>
<li><b>CAP:</b> CAP elevation (ft) × 0.3048, plus the mean offset between each CAP estimate and
 the nearest valid SWOT pass in the same year. Years without SWOT data are not plotted.</li>
<li><b>Pixel cloud:</b> {meta['pixc_label']} from year
 {s['pixc_min_year'] or meta['now'].year} onward, latest {s['max_3d_passes']} pass dates per lake.</li>
<li><b>Maps:</b> Esri World Imagery / World Topographic basemaps; lake outlines from the
 report's GeoJSON files where available, otherwise the assessment's lake coordinates.</li>
<li><b>Generated by</b> <code>gdl_swot_report.py</code> v{__version__} on {meta['gen_txt']}.</li>
</ul></section>"""

    if plotly_js_mode == "cdn":
        pjs = (f'<script src="https://cdn.plot.ly/plotly-{get_plotlyjs_version()}.min.js" '
               'charset="utf-8"></script>')
    else:
        pjs = f'<script type="text/javascript">{get_plotlyjs()}</script>'
    leaflet = (f'<link rel="stylesheet" href="https://unpkg.com/leaflet@{LEAFLET_VERSION}/dist/leaflet.css">'
               f'<script src="https://unpkg.com/leaflet@{LEAFLET_VERSION}/dist/leaflet.js"></script>')

    specs = "\n".join(f'<script type="application/json" id="fig-{fid}">{js}</script>'
                      for fid, js in H.fig_specs)
    map_json = json.dumps(map_payload(lakes_out, rep)).replace("</", "<\\/")

    return f"""<!DOCTYPE html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{esc(rep['title'])}</title>
<style>{CSS}</style>
{leaflet}
{pjs}
</head><body><div id="page">
<header id="banner"><a id="top"></a>
  <div class="t1">{esc(rep['author'])}</div>
  <div class="t2">{esc(rep['title'])}</div>
  <div class="t3">SWOT satellite water-surface elevations for glacier-dammed lakes ·
  data through {meta['data_through']:%Y-%m-%d %H:%M} UTC</div>
</header>
<nav id="topnav"><a href="#overview">Overview</a><a href="#overview-map">Map</a>
<a href="#qa">Quality assessment</a><a href="#notes">Processing notes</a>
<a class="ext" href="{GDL_MAIN_URL}" target="_blank" rel="noopener">NWS APRFC GDL page &#8599;</a></nav>
<div id="crumb">Alaska &rsaquo; Glacier-dammed lakes &rsaquo; SWOT monitoring</div>
<div id="wrap">
<aside id="side">{''.join(side)}</aside>
<main>
<section id="overview">
  <h2>Overview</h2>
  <p class="help">{"Lakes shown are glacier-dammed lakes where SWOT was assessed as "
  + esc(" or ".join(c.lower() for c in s["include_classes"])) + ". " if s["include_classes"] else ""}Click a lake on the map or in the table to jump to it; click a
  column header to sort. <b>Age</b> is days since the latest valid pass.
  Source: {meta['source_html']}.</p>
  <div id="overview-map" class="map map-lg" aria-label="Map of lakes"></div>
  {overview_table(lakes_out, s)}
</section>
{''.join(sections)}
{qa_sec}
{notes}
</main>
</div>
<footer>
<p><b>Provisional data.</b> SWOT-derived elevations and the quality assessment are
provisional and subject to revision. This page is not a National Weather Service product;
for official glacier-dammed lake information see the
<a href="{GDL_MAIN_URL}">NWS Alaska-Pacific River Forecast Center</a>.</p>
<p>{esc(rep['author'])} · Basemaps &copy; Esri · Generated {esc(meta['gen_txt'])}</p>
</footer>
</div>
<script type="application/json" id="map-data">{map_json}</script>
{specs}
<script>{JS}</script>
</body></html>"""


# =============================================================================
# PUBLISH (GitHub Pages)
# =============================================================================

def _git(repo: Path | None, *args: str) -> str:
    cmd = ["git"] + (["-C", str(repo)] if repo else []) + list(args)
    log.debug("$ %s", " ".join(cmd))
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(f"`{' '.join(cmd)}` failed:\n{(r.stderr or r.stdout).strip()}")
    return r.stdout


def prepare_repo(pub: dict, base: Path) -> Path:
    """Make sure a clean, up-to-date clone of the Pages repo exists; return its path."""
    if shutil.which("git") is None:
        raise RuntimeError("git is not installed or not on PATH.")
    d = Path(pub["repo_dir"])
    d = d if d.is_absolute() else base / d
    br = pub["branch"]
    if not (d / ".git").exists():
        if d.exists() and any(d.iterdir()):
            raise RuntimeError(f"{d} exists but is not a git clone; set [publish].repo_dir.")
        log.info("Cloning %s into %s", pub["repo_url"], d)
        _git(None, "clone", "--branch", br, pub["repo_url"], str(d))
    else:
        log.info("Updating %s", d)
        _git(d, "checkout", br)
        _git(d, "pull", "--ff-only", "origin", br)
    return d


def publish_repo(d: Path, pub: dict, push: bool, when: pd.Timestamp) -> None:
    """Commit the report file and push it (the Pages workflow redeploys on push)."""
    f = pub["file"]
    _git(d, "add", "--", f)
    if not _git(d, "status", "--porcelain", "--", f).strip():
        log.info("%s unchanged; nothing to publish.", f)
        return
    _git(d, "commit", "-m", pub["commit_message"].format(date=f"{when:%Y-%m-%d}"), "--", f)
    if not push:
        log.info("Committed %s in %s (not pushed: --no-push).", f, d)
        return
    log.info("Pushing to %s (%s)…", pub["repo_url"], pub["branch"])
    _git(d, "push", "origin", pub["branch"])
    m = re.search(r"github\.com[/:]([^/]+)/([^/.]+)", pub["repo_url"])
    if m:
        log.info("Pushed. The site updates in a minute or two: https://%s.github.io/%s/",
                 m.group(1).lower(), m.group(2))


# =============================================================================
# MAIN
# =============================================================================

def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description="Build an HTML SWOT monitoring report for glacier-dammed lakes.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("-c", "--config", type=Path, default=Path("lakes.toml"),
                   help="TOML config file with settings and [[lakes]] entries")
    p.add_argument("-o", "--output", type=Path, default=None,
                   help="output HTML file (default: index.html; with --publish, "
                        "the file inside the GitHub Pages clone)")
    p.add_argument("--sciencebase-item", metavar="ID_OR_URL",
                   help="download inputs from this ScienceBase item "
                        "(overrides [sciencebase].item_id)")
    p.add_argument("--download-dir", type=Path,
                   help="where ScienceBase files are cached")
    p.add_argument("--refresh", action="store_true",
                   help="re-download ScienceBase files even if cached")
    p.add_argument("--sb-list", action="store_true",
                   help="list the files on the ScienceBase item and exit")
    p.add_argument("--offline", action="store_true",
                   help="don't contact ScienceBase; use cached/local files only")
    p.add_argument("--data-dir", type=Path, action="append", default=[],
                   help="extra folder to search for input files (repeatable)")
    p.add_argument("--lakes", nargs="+", metavar="NAME",
                   help="only these lakes (names as in the CSV lake_name column)")
    p.add_argument("--all-lakes", action="store_true",
                   help="also include lakes found in the CSV but not in the config")
    p.add_argument("--all-classes", action="store_true",
                   help="don't limit the report to glacier-dammed lakes assessed as "
                        "likely/potentially successful ([settings].include_classes)")
    p.add_argument("--no-3d", action="store_true",
                   help="skip the 3-D pixel-cloud plots (smaller, faster)")
    p.add_argument("--plotly-js", choices=["inline", "cdn"], default="inline",
                   help="embed plotly.js (works offline) or load it from the CDN")
    p.add_argument("--report-date", type=dt.date.fromisoformat,
                   help="treat this date (YYYY-MM-DD) as 'today'")
    p.add_argument("--list-lakes", action="store_true",
                   help="print the lakes found in the aggregated CSV and exit")
    p.add_argument("--publish", action="store_true",
                   help="write the report into the GitHub Pages repo clone "
                        "([publish] in the config), commit, and push")
    p.add_argument("--no-push", action="store_true",
                   help="with --publish: commit but don't push")
    p.add_argument("--timeseries", type=Path, metavar="PARQUET",
                   help="read data/aggregated/wse_timeseries.parquet from this repo's "
                        "pipeline instead of aggregated CSVs (ScienceBase is skipped)")
    p.add_argument("--pixc-root", type=Path, metavar="DIR",
                   help="pipeline pixel files (data/pixc) for the 3-D plots")
    p.add_argument("--wse-set", choices=list(WSE_SETS), default="polygon_qc",
                   help="with --timeseries: which pixel set's WSE to report")
    p.add_argument("-v", "--verbose", action="store_true")
    return p.parse_args(argv)


def main(argv=None) -> int:
    a = parse_args(argv)
    logging.basicConfig(level=logging.DEBUG if a.verbose else logging.INFO,
                        format="%(levelname)s: %(message)s")

    if not a.config.exists():
        log.error("Config file not found: %s", a.config)
        return 2
    cfg = load_config(a.config)
    s, base = cfg["settings"], cfg["base_dir"]
    tz = ZoneInfo(cfg["report"]["display_timezone"])

    pub_dir = None
    if a.publish:
        try:
            pub_dir = prepare_repo(cfg["publish"], base)
        except RuntimeError as e:
            log.error("Publish setup failed: %s", e)
            return 3

    now = (pd.Timestamp(a.report_date) + pd.Timedelta(hours=23, minutes=59)
           if a.report_date else pd.Timestamp.now(tz="UTC").tz_localize(None))

    # ---- locate inputs -----------------------------------------------------
    dirs = [d if d.is_absolute() else Path.cwd() / d for d in a.data_dir]
    dirs += [base / d for d in cfg["data"]["search_dirs"]]

    sb_id = None if a.timeseries else (a.sciencebase_item or cfg["sciencebase"]["item_id"])
    source_html = "local files"
    sb_msgs = []
    if sb_id:
        sb_id = parse_sb_id(sb_id)
        dl = (a.download_dir or base / cfg["sciencebase"]["download_dir"]) / sb_id
        url = f"https://www.sciencebase.gov/catalog/item/{sb_id}"
        source_html = f'ScienceBase item <a href="{url}">{sb_id}</a>'
        if not a.offline:
            try:
                item, entries = sciencebase_catalog(sb_id, cfg["sciencebase"]["include_children"])
                if a.sb_list:
                    print_sb_listing(item, entries)
                    return 0
                log.info("ScienceBase: %s (%d files)", item.get("title", sb_id), len(entries))
                upd = next((d.get("dateString") for d in item.get("dates") or []
                            if "revis" in str(d.get("type", "")).lower()
                            or "update" in str(d.get("type", "")).lower()), None)
                title = esc(item.get("title", ""))
                source_html = (f'ScienceBase item <a href="{url}">{title or sb_id}</a>'
                               + (f" (last updated {esc(upd)})" if upd else ""))
                cfg_lakes = cfg["lakes"]
                if a.lakes:
                    want = {name_key(n) for n in a.lakes}
                    cfg_lakes = [lk for lk in cfg_lakes if name_key(lk["name"]) in want] + \
                        [default_lake(n) for n in a.lakes
                         if name_key(n) not in {name_key(l["name"]) for l in cfg_lakes}]
                picked, sb_msgs = select_sb_files(
                    entries, cfg_lakes, a.all_lakes, not a.no_3d,
                    float(cfg["sciencebase"]["max_small_file_mb"]))
                for m in sb_msgs:
                    log.warning(m)
                download_sciencebase(picked, dl, a.refresh)
            except Exception as e:  # network errors: fall back to cache / local
                if a.sb_list:
                    log.error("Could not read ScienceBase item %s: %s", sb_id, e)
                    return 2
                log.warning("ScienceBase unavailable (%s); using cached/local files.", e)
        if dl.is_dir():
            dirs.insert(0, dl)
    elif a.sb_list:
        log.error("--sb-list needs --sciencebase-item or [sciencebase].item_id.")
        return 2
    finder = FileFinder(dirs)

    if a.timeseries:
        if not a.timeseries.exists():
            log.error("Time series not found: %s", a.timeseries)
            return 2
        log.info("Pipeline time series: %s (%s pixels)", a.timeseries, a.wse_set)
        agg = unify_lake_names(read_timeseries(a.timeseries, a.wse_set), cfg["lakes"])
        csv_label = a.timeseries.name
        source_html = (f"SWOT lake pipeline, <code>{esc(csv_label)}</code> "
                       f"(<code>{esc(a.wse_set)}</code> pixels)")
        wse_label = f"the median of the <code>{esc(a.wse_set)}</code> pixels"
        sd_label = "their standard deviation"
        pixc_label = {
            "refined_qc": "Segmentation-assigned, quality-filtered pixels drawn at their "
                          "refined positions (<code>refined_qc</code>)",
            "polygon_qc": "Quality-filtered pixels whose reported position falls inside the "
                          "lake polygon (<code>in_polygon_qc</code>; not the segmentation "
                          "assignment)",
            "polygon_all": "All pixels inside the lake polygon, unfiltered "
                           "(<code>in_polygon</code>)",
        }[a.wse_set]
    else:
        csv_paths = find_aggregated(finder, cfg["data"]["aggregated_csv"], cfg["columns"],
                                    dl if sb_id else None)
        if not csv_paths:
            log.error("No aggregated WSE CSV (%s, or any CSV with columns %s) in: %s",
                      ", ".join(cfg["data"]["aggregated_csv"]), ", ".join(AGG_REQUIRED),
                      ", ".join(str(d) for d in finder.dirs) or "(no existing folders)")
            return 2
        for pth in csv_paths:
            log.info("Aggregated WSE: %s", pth)
        agg = unify_lake_names(read_aggregated_many(csv_paths, cfg["columns"]), cfg["lakes"])
        csv_label = ", ".join(sorted({pth.name for pth in csv_paths}))
        if not sb_id:
            source_html = f"<code>{esc(csv_label)}</code> (local)"
        else:
            local = sorted({pth.name for pth in csv_paths if dl not in pth.parents})
            if local:
                source_html += f", combined with local <code>{esc(', '.join(local))}</code>"
        wse_label = f"the filtered median, <code>{WSE_COL}</code>"
        sd_label = f"<code>{SD_COL}</code>"
        pixc_label = "Filtered pixels (<code>filtered_boolean == TRUE</code>)"

    if a.list_lakes:
        c = agg.groupby("lake_name")["time"].agg(["count", "min", "max"])
        c.columns = ["passes", "first (UTC)", "last (UTC)"]
        print(c.to_string())
        return 0

    # ---- which lakes -------------------------------------------------------
    lakes = list(cfg["lakes"])
    known = {name_key(lk["name"]) for lk in lakes}
    if a.all_lakes or not lakes:
        lakes += [default_lake(n) for n in sorted(agg["lake_name"].unique())
                  if name_key(n) not in known]
    if a.lakes:
        want = {name_key(n) for n in a.lakes}
        lakes = [lk for lk in lakes if name_key(lk["name"]) in want]
        missing = want - {name_key(lk["name"]) for lk in lakes}
        csv_names = {name_key(n): n for n in agg["lake_name"].unique()}
        lakes += [default_lake(csv_names.get(name_key(n), n)) for n in a.lakes
                  if name_key(n) in missing]
    # ---- keep only assessed glacier-dammed lakes in the wanted classes -----
    for lk in lakes:
        lk["qa"] = lake_assessment(lk)
    if a.all_classes:
        log.info("--all-classes: not filtering lakes by quality assessment.")
        s["include_classes"] = []
    else:
        keep = {c.lower() for c in s["include_classes"]}
        kept = []
        for lk in lakes:
            cls = (lk["qa"] or {}).get("classification")
            if cls and cls.lower() in keep:
                kept.append(lk)
            else:
                log.info("  skip %s (%s)", lk["name"],
                         f"assessed '{cls}'" if cls else "not an assessed glacier-dammed lake")
        lakes = kept
        found = {lk["qa"]["name"] for lk in lakes if lk["qa"]}
        for nm, row in GDL_ASSESSMENT.items():
            if row[-1].lower() in keep and nm not in found:
                log.warning("Assessed lake '%s' (%s) not found in the data; set "
                            "aliases = [\"<name in data>\"] for it in the config.", nm, row[-1])
    lakes.sort(key=lambda lk: (QA_RANK.get((lk["qa"] or {}).get("classification"), 9),
                               lk["display_name"].lower()))
    if not lakes:
        log.error("No lakes selected.")
        return 2

    # ---- build each lake ---------------------------------------------------
    H = Html()
    lakes_out = []
    for lk in lakes:
        log.info("Lake %s", lk["name"])
        warnings = []
        raw = agg[agg["lake_name"] == lk["name"]]
        if raw.empty:
            loose = agg[agg["lake_name"].map(name_key) == name_key(lk["name"])]
            if len(loose):
                raw = loose
                log.info("  matched CSV lake_name %s", sorted(loose["lake_name"].unique()))
        if raw.empty:
            import difflib
            close = difflib.get_close_matches(lk["name"], sorted(agg["lake_name"].unique()), 3, 0.5)
            hint = f" Did you mean {', '.join(repr(c) for c in close)}?" if close else ""
            warnings.append(f"No rows with lake_name '{lk['name']}' in {csv_label}.{hint}")
            log.warning("  %s", warnings[-1])
        passes, notes = classify_lake(raw, s)
        warnings += notes
        used = passes[passes["status"] == STATUS_USED]
        cap = adjust_cap(lk["cap"], used)
        predicted = predict_overpasses(raw, now, s["predict_days_ahead"])
        summary = lake_summary(used, lk, s, now, predicted)

        if "latest" in summary and summary["latest"]["n_pix_filtered"] < 50:
            warnings.append(f"Latest pass has only {int(summary['latest']['n_pix_filtered'])} "
                            "pixels — treat with caution.")
        n_rej = int((passes["status"] == STATUS_OUTLIER).sum())
        if n_rej:
            warnings.append(f"{n_rej} pass(es) removed as outliers — see overpass table.")

        # Lake outline: used by the maps and the 3-D plot.
        try:
            ol_path, ol_sel = find_outline(finder, lk)
            rings = read_outline(ol_path, ol_sel)
        except Exception as e:
            ol_path, ol_sel, rings = None, None, []
            log.warning("  could not read outline: %s", e)

        fig_p, n3 = None, 0
        if not a.no_3d:
            pq_paths = (pipeline_pixc_files(a.pixc_root, lk["name"], s, now) if a.pixc_root
                        else find_pixc(finder, lk, cfg["columns"]))
            if not pq_paths:
                warnings.append(f"No pixel-cloud file found for {lk['name']}; no 3-D plot.")
            else:
                log.info("  pixel cloud: %s", ", ".join(pth.name for pth in pq_paths[:5])
                         + (f" (+{len(pq_paths) - 5} more)" if len(pq_paths) > 5 else ""))
                if ol_path is None:
                    warnings.append("Lake outline not found; 3-D plot has no outline.")
                try:
                    pixc = (read_pipeline_pixc(pq_paths, s, a.wse_set) if a.pixc_root
                            else read_pixc(pq_paths, s, now, cfg["columns"]))
                    fig_p = fig_pixc(pixc, rings, lk, s)
                    n3 = min(pixc["date"].nunique(), s["max_3d_passes"]) if len(pixc) else 0
                    if fig_p is None:
                        warnings.append("No filtered PIXC points in the window for the 3-D plot.")
                except Exception as e:
                    warnings.append(f"Could not build 3-D plot: {e}")
                    log.exception("3-D plot failed for %s", lk["name"])

        lid = "lake-" + slug(lk["name"])
        lakes_out.append({
            "lake": lk, "id": lid, "summary": summary, "warnings": warnings,
            "qa": lk.get("qa"), "rings": rings,
            "fig_ts": fig_timeseries(used, cap, lk, s, now) if len(used) else None,
            "fig_se": fig_seasonal(used, cap, lk, s, now) if len(used) else None,
            "fig_pass": fig_by_pass(used, s, now, lk) if len(used) else None,
            "fig_pixc": fig_p, "n_pixc_passes": n3,
            "table_passes": overpass_table(lid, passes.rename(columns={"pass": "pass_"}),
                                           predicted.rename(columns={"pass": "pass_"}),
                                           tz, lk["name"]),
        })

    gen_utc = pd.Timestamp.now(tz="UTC")
    gen_txt = (f"{gen_utc:%Y-%m-%d %H:%M} UTC "
               f"({gen_utc.tz_convert(tz):%b %d %H:%M %Z})")
    if a.report_date:
        gen_txt += f", as of {a.report_date:%Y-%m-%d}"
    meta = {"tz": tz, "html": H, "now": now, "gen_txt": gen_txt,
            "source_html": source_html, "data_through": agg["time"].max(),
            "wse_label": wse_label, "sd_label": sd_label, "pixc_label": pixc_label}

    out = a.output or Path("index.html")
    if pub_dir is not None and a.output is None:
        out = pub_dir / cfg["publish"]["file"]
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(build_html(cfg, lakes_out, meta, a.plotly_js), encoding="utf-8")
    log.info("Wrote %s (%.1f MB, %d lakes)", out, out.stat().st_size / 1e6, len(lakes_out))

    if pub_dir is not None:
        target = pub_dir / cfg["publish"]["file"]
        if out.resolve() != target.resolve():
            shutil.copyfile(out, target)
        try:
            publish_repo(pub_dir, cfg["publish"], not a.no_push, now)
        except RuntimeError as e:
            log.error("Publish failed: %s", e)
            return 3
    return 0


if __name__ == "__main__":
    sys.exit(main())