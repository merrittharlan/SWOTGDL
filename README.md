# SWOT Glacial Dammed Lakes Data Pipeline

## Project status
Project is under development - for internal use only.

## Description
Glacial-dammed lakes (GDLs) in Alaska are important features of the cryosphere, posing significant risks due to potential outburst floods that can cause downstream damage and influence hydrology and ecosystems. This study explores how new data from the Surface Water Ocean Topography (SWOT) satellite, launched in 2022, can be used to monitor and track changes in water elevation of glacial-dammed lakes in Alaska. The SWOT satellite provides high-resolution, global measurements of water surface elevation for large waterbodies, offering a unique capability for monitoring these remote, inaccessible regions. The study, of which this codebase is a small part, demonstrates the potential use of SWOT pixel cloud data to detect changes in elevations in GDLs.

This repository runs a scheduled GitLab pipeline that downloads SWOT pixel cloud (PIXC) data over the monitored lakes, assigns pixels to each lake, and maintains a water surface elevation (WSE) time series per lake and overpass.

```
config.yaml              all settings (collection, season, QC thresholds, assignment parameters)
src/export_pixc.py       search + download new granules, save the pixels for each lake
src/aggregate_wse.py     per-overpass WSE statistics from the saved pixels
src/gdl_swot_report.py   HTML report of the time series, published with GitLab Pages
src/common.py            shared helpers
report.toml              report settings, and optional per-lake notes, releases and CAP estimates
lakes/SWOT_monitoring_lakes.geojson   lake polygons and properties (add this)
```

Run everything from the repo root, so the scripts find `config.yaml`.

## Related Publications
Monitoring Glacial Dammed Lakes with the SWOT Satellite (Paper -- in progress)

[Satellite Glacial Dammed Lake Elevations from Surface Water Ocean Topography Satellite from April 2023 to May 2025](https://www.sciencebase.gov/catalog/item/67d1cc1fd34ecfe34cc88723) (Data Release)

[Water Information from Space (WISP)](https://apps.usgs.gov/wisp) (Web Application)

## Installation
The scheduled pipeline installs its own dependencies. To run the scripts locally:

1. <strong>Clone or download the repository</strong>:
    ```bash
    $ git clone <repository URL>
    $ cd <repository folder>
    ```

2. <strong>Set up the environment</strong>:

    Using Windows CMD:
    ```cmd
    $ python -m venv ./env
    $ .\env\Scripts\activate
    $ pip install -r requirements.txt
    ```
    OR using Linux/MacOS:
    ```bash
    $ python -m venv ./env
    $ source ./env/bin/activate
    $ pip install -r requirements.txt
    ```

## Output layout

```
data/
  processed_granules.csv               one row per (lake, granule): ok / no_pixels / error
  pixc/<lake>/<granule>.parquet        pixels for one lake from one granule (written once, never changed)
  aggregated/wse_timeseries.parquet    one row per lake and overpass
public/index.html                      the report (built in CI and published to Pages, not committed)
```

Each run only downloads granules that aren't in `processed_granules.csv`. That includes
granules that had no pixels and granules that failed (failed ones are retried). To
reprocess something, delete its rows from the manifest, or delete `data/` to start over.
Because the pixel files never change once written, git history doesn't grow by a full
copy of the dataset on every run.

## GitLab setup

1. **Add the lakes file** at `lakes/SWOT_monitoring_lakes.geojson`. Each feature needs
   `name`, `sw_lon/sw_lat/ne_lon/ne_lat`, `elev_min`, `elev_max` (NAVD88),
   `NAVD88-EGM2008_m` and `passes` (an empty list means any pass).
2. **CI/CD variables** (Settings > CI/CD > Variables, *Masked*):
   - `EARTHDATA_USERNAME`, `EARTHDATA_PASSWORD`
   - `DATA_PUSH_TOKEN`: a project access token (Settings > Access tokens) with the
     `write_repository` scope and the Developer role or higher. If the branch is protected,
     allow that role to push to it.
3. **Schedule** (Build > Pipeline schedules), e.g. daily `0 14 * * *`. Outside the
   June–November season the job finds nothing new and finishes in a minute or two.
4. **First run:** the pipeline starts with an empty `data/` folder and downloads the full
   record since `start_date` in `config.yaml`. Each run stops starting new downloads after
   `TIME_LIMIT_MIN` (90) and commits what it has, so this first download finishes over
   several runs. You can also press *Run pipeline* a few times to speed it up. If your
   runners allow a longer job timeout, raise `timeout` and `TIME_LIMIT_MIN`.

Before turning on the schedule, it's worth a small local test (below) on one lake, to check
your Earthdata login and the lakes GeoJSON fields.

## GitLab Pages report

After each data update, the CI job builds `public/index.html` with `src/gdl_swot_report.py`
and the `pages` job publishes it. The report shows every lake in the data: a forecaster
overview table, each lake's WSE time series and seasonal comparison, an overpass table
with predicted next passes, and a 3-D pixel cloud for the current season. It plots the
`polygon_qc` WSE; add `--wse-set refined_qc` to the report command in `.gitlab-ci.yml` to
plot the refined assignment instead.

Setup:
1. Check that Pages is enabled: *Settings → General → Visibility, project features,
   permissions → Pages*. Because the project is for internal use, set Pages access to
   **Only Project Members** there.
2. After the first successful pipeline, *Deploy → Pages* shows the site URL.
3. To add display names, notes, documented releases or CAP estimates for a lake, add a
   `[[lakes]]` block to `report.toml` (there's an example in the file).

The page updates every time the scheduled pipeline runs.

## Running locally

```bash
python src/export_pixc.py --lakes "Suicide Basin" --max-granules 5   # small test on one lake
python src/export_pixc.py                    # everything; uses ~/.netrc or prompts for Earthdata login
python src/aggregate_wse.py
python src/gdl_swot_report.py --config report.toml --all-lakes \
    --timeseries data/aggregated/wse_timeseries.parquet --pixc-root data/pixc -o public/index.html
```

The report needs Python 3.11 or newer.

Delete `data/` to start over from scratch.

For the `uncert_*` columns, clone https://github.com/SWOTAlgorithms/RiverObs and
`export PYTHONPATH=RiverObs/src`. Only `SWOTWater/aggregate.py` is used, so the full
package doesn't need to be installed.

## Column names

Pixel files keep the PIXC product variable names (`height`, `geoid`, `sig0`,
`phase_noise_std`, ...), so they match the SWOT product description. Lake metadata
(area, region, TRI, elevation bounds) is no longer repeated on every pixel or overpass.
It lives in the lakes GeoJSON; join on `lake_name` (the GeoJSON's `name`) when you need it.

| Old (pixel) | New |
|---|---|
| `elevation_m_EGM2008` / `elevation_m_NAVD88` | `wse_m_egm2008` / `wse_m_navd88` |
| `latitude_WGS84`, `longitude_WGS84` (+ `_refined`) | `latitude`, `longitude`, `latitude_refined`, `longitude_refined` |
| `height_m`, `geoid_m`, `*_tide_*_m`, `pixel_area_sqm`, `cross_track_dist_m` | `height`, `geoid`, `*_tide_*`, `pixel_area`, `cross_track` |
| `naive_in_polygon` | `in_polygon` |
| `filtered_boolean` | `in_polygon_qc` |
| `fancy_assigned` / `fancy_filtered_boolean` / `fancy_status` | `refined_assigned` / `refined_qc` / `refined_status` |
| `datetime_UTC`, `pixc_url`, `pixc_crid` | `datetime_utc`, `granule`, `crid` (+ `cycle_number`, `pass_number`, `tile_name`) |

Aggregated columns are `wse_m_navd88_<set>_<stat>` and `n_pix_<set>`:

| set | pixels | old name |
|---|---|---|
| `polygon_all` | in the lake polygon, no QC | `all` |
| `polygon_qc` | in the lake polygon, passing QC | `filtered` |
| `refined_qc` | refined assignment, passing QC | `fancy` |

| stat | meaning | old name |
|---|---|---|
| `median`, `sd` | median and standard deviation of pixel WSE | `med`, `sd` |
| `weighted_mean` | inverse-variance weighted mean, weights 1/(phase_noise_std·dheight_dphase)² | `weighted_med` |
| `uncert_std`, `uncert_multilook` | RiverObs aggregate.py uncertainties | same |

`wse_m_navd88_refined_minus_polygon_qc_median` replaces `wse_m_NAVD88_fancy_minus_filtered_med`.

## Authors

[Alexander Hernes](https://www.usgs.gov/staff-profiles/alek-hernes) - Hydrologic Technician, USGS California Water Science Center | ahernes@usgs.gov

[Merritt Harlan](https://www.usgs.gov/staff-profiles/merritt-harlan) - Research Hydrologist, USGS Water Resources Mission Area | mharlan@usgs.gov


## Code of Conduct

All contributions to- and interactions surrounding- this project will
abide by the [USGS Code of Scientific
Conduct](https://www.usgs.gov/fsp/fsp_code_of_scientific_conduct.asp).

## Support

This package is under development and extended support cannot be guaranteed.

## Disclaimer

[Disclaimer](DISCLAIMER.md)  
This software is preliminary or provisional and is subject to revision.
It is being provided to meet the need for timely best science. The
software has not received final approval by the U.S. Geological Survey
(USGS). No warranty, expressed or implied, is made by the USGS or the
U.S. Government as to the functionality of the software and related
material nor shall the fact of release constitute any such warranty. The
software is provided on the condition that neither the USGS nor the U.S.
Government shall be held liable for any damages resulting from the
authorized or unauthorized use of the software.

This software is in the public domain because it contains materials that originally came from the U.S. Geological Survey, an agency of the United States Department of Interior. For more information, see the [official USGS copyright policy](https://www.usgs.gov/visual-id/credit_usgs.html#copyright)

Although this software program has been used by the U.S. Geological Survey (USGS), no warranty, expressed or implied, is made by the USGS or the U.S. Government as to the accuracy and functioning of the program and related program material nor shall the fact of distribution constitute any such warranty, and no responsibility is assumed by the USGS in connection therewith.

Any use of trade, firm, or product names is for descriptive purposes only and does not imply endorsement by the U.S. Government.

This software is provided "AS IS."

[![CC0](http://i.creativecommons.org/p/zero/1.0/88x31.png)](http://creativecommons.org/publicdomain/zero/1.0/)

## License Information


[License](LICENSE.md)
This software is a product of the U.S. Geological Survey (USGS), an
agency of the United States Department of Interior, which is part of the
U.S. Government.

#### Cost

This software is freely distributed. There is no fee to download and
(or) use this software.

#### License

Users do not need a license or permission from the USGS to use this
software. Users can download and install as many copies of the software
as they need. See Exceptions below.

#### Public domain

As a work of the United States Government, this USGS product is in the
[public domain](http://www.usa.gov/publicdomain/label/1.0/) within the
United States. You can copy, modify, distribute, and perform the work,
even for commercial purposes, all without asking permission.
Additionally, USGS waives copyright and related rights in the work
worldwide through [CC0 1.0 Universal Public Domain
Dedication](https://creativecommons.org/publicdomain/zero/1.0/).  
[![CC0](http://i.creativecommons.org/p/zero/1.0/88x31.png)](http://creativecommons.org/publicdomain/zero/1.0/)

#### Software User Rights Notice

This software is a product of the U.S. Geological Survey, which is part of the U.S. Government. This software is freely distributed. There is no fee to download and (or) use this software. 

As a work of the United States Government, this USGS product is in the public domain (https://www.usa.gov/government-works) within the United States. You can copy, modify, distribute, and perform the work, even for commercial purposes, all without asking permission. Additionally, USGS waives copyright and related rights in the work worldwide through CC0 1.0 Universal Public Domain Dedication (https://creativecommons.org/publicdomain/zero/1.0/ ). 

This software has been approved for release by the U.S. Geological Survey (USGS). Although the software has been subjected to rigorous review, the USGS reserves the right to update the software as needed pursuant to further analysis and review. No warranty, expressed or implied, is made by the USGS or the U.S. Government as to the functionality of the software and related material nor shall the fact of release constitute any such warranty. Furthermore, the software is released on condition that neither the USGS nor the U.S. Government shall be held liable for any damages resulting from its authorized or unauthorized use. 

Any use of trade, firm, or product names is for descriptive purposes only and does not imply endorsement by the U.S. Government. 

#### Exceptions

This project may have numerous dependencies from other open-source
software projects. The re-use and distribution of those software
packages may be subject to the licenses of each of those dependencies.


## Additional Information 

The Digital Object Identifier for this project is: <br>
The IPDS Record Number is: IP-180833
