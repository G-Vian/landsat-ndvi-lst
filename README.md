# landsat-c2-zonal

**NDVI and land surface temperature time series from Landsat Collection 2,
aggregated to any set of polygons.**

A modular R pipeline that turns raw **Landsat Collection 2 Level-2** scenes into
analysis-ready **NDVI** and **LST** series for *any* area of interest, at three
spatial scales (whole AOI, per zone, per pixel) and two temporal scales
(monthly, annual).

Handles **Landsat 4, 5, 7, 8 and 9 in a single run**: band numbering,
per-scene calibration factors and sensor-specific quality flags are resolved
automatically.

Adapting it to a new study area means editing **one file**.

---

## Table of contents

1. [Can I use this on my own area?](#1-can-i-use-this-on-my-own-area)
2. [Why this pipeline](#2-why-this-pipeline)
3. [Requirements](#3-requirements)
4. [Quick start](#4-quick-start)
5. [Verifying the installation](#5-verifying-the-installation)
6. [Input data layout](#6-input-data-layout)
7. [Configuration reference](#7-configuration-reference)
8. [Adapting to a new study area](#8-adapting-to-a-new-study-area)
9. [Module reference](#9-module-reference)
10. [Outputs](#10-outputs)
11. [Output column dictionary](#11-output-column-dictionary)
12. [Methods](#12-methods)
13. [Quality control cascade](#13-quality-control-cascade)
14. [Running on an HPC cluster](#14-running-on-an-hpc-cluster)
15. [Troubleshooting](#15-troubleshooting)
16. [Known limitations](#16-known-limitations)
17. [References](#17-references)
18. [Citation and licence](#18-citation-and-licence)

---

## 1. Can I use this on my own area?

**Yes.** Nothing in the code is specific to the area it was first written for.
The automated self-test builds a synthetic study area from scratch and runs the
whole pipeline over it, so "works anywhere" is verified rather than claimed.

### What you provide

| | Requirement |
|---|---|
| **A polygon file** | `.shp`, `.gpkg` or `.geojson`. Must have a **defined CRS**, **non-overlapping** polygons, and a column holding the zone names. |
| **Landsat scenes** | Collection 2 **Level-2** products from [EarthExplorer](https://earthexplorer.usgs.gov/), each unpacked into its own folder, with the red, NIR, thermal and `QA_PIXEL` bands plus the `MTL` file. |
| **Five edits** | In `R/00_config.R` only. Listed below. |

### The five settings you must set

```r
# 1. Paths and identity
AOI_NAME <- "MyCity";  BASE_DIR <- "..."
DIRS_LANDSAT <- list(all = file.path(BASE_DIR, "landsat"))
SHP_FILE <- file.path(BASE_DIR, "shp/zones.gpkg")
ZONE_NAME_COLUMN <- "district_name"       # set it explicitly

# 2. Period
YEAR_START <- 2015;  YEAR_END <- 2024

# 3. Projection - MUST be projected (metres), not degrees
CRS_TARGET <- "EPSG:32723"

# 4. Physical LST bounds - CLIMATE-DEPENDENT   <-- silent trap
# Shipped defaults are deliberately WIDE (-20, 75) so that a first run on an
# unknown climate discards nothing. Narrow them to your region afterwards.
LST_MIN <- -20.0;  LST_MAX <- 75.0

# 5. Scene/AOI overlap - SIZE-DEPENDENT        <-- silent trap
MIN_AOI_OVERLAP_FRAC <- 0.10
```

### ⚠ Two settings fail *silently* if you leave the defaults

These do not raise an error. They produce a wrong or empty result that looks
plausible. They are the reason this section exists.

**`LST_MIN` / `LST_MAX` — wrong bounds truncate your series without warning.**
The shipped defaults (`-20`, `75`) are wide on purpose: they let a first run
complete anywhere without silently discarding data. They are *not* the right
final values. Set them too narrow and the failure is invisible — a floor of
`5 °C` in a temperate region discards *every genuine winter observation*, and
you get a series that looks fine while systematically missing its cold half.
Set them too wide and artefacts survive into the statistics.

| Climate | `LST_MIN` | `LST_MAX` |
|---|---|---|
| Tropical / subtropical coastal | 5 | 70 |
| Tropical inland / semi-arid | 5 | 75 |
| Temperate (with frost) | −15 | 60 |
| Boreal / high latitude | −45 | 45 |
| Desert | 0 | 80 |

*How to check:* after the first run, open `tables/anomalies/`. Pixels piling up
*at* a bound mean the bound is cutting into real data. Module 09 also warns at
runtime if more than half the scenes hit a limit.

**`MIN_AOI_OVERLAP_FRAC` — too high for a large area rejects every scene.**
One Landsat scene covers ~185 × 180 km. For a city it covers 100% of your area
and the default works. For a state, each scene covers a slice, so a 10%
threshold rejects *all* of them and the run ends with nothing.

| Your area vs one scene | Value |
|---|---|
| Smaller (city, park, watershed) | `0.10` |
| A few scenes (metro region, small state) | `0.01` |
| Many scenes (large state, country) | `0.0001` or `0` |

*Symptom:* the log repeats "no sufficient overlap" for every scene.

### What this pipeline is and is not suited to

| Scale | Status |
|---|---|
| City, municipality, watershed, park | **Ideal.** Defaults work; per-pixel tables tractable. |
| Metro region, small state | **Fine** with the overlap threshold lowered. |
| Large state, province, country | **Works**, but set `COMPUTE_PIXEL_TABLES <- FALSE`, use an equal-area CRS, and expect seams (below). |

**It does not build a mosaic.** Scenes are processed independently and
aggregated statistically. That is *correct* for zonal statistics — each zone
gets a pixel-count weighted mean of whichever scenes cover it. But the
per-pixel maps of an area spanning several scenes can show scene-boundary
seams. If a seamless composite image is your deliverable, this is the wrong
tool.

### Before you trust any result

```bash
Rscript tests/run_tests.R     # 33 checks, no real data needed, ~1 minute
```

Then do a **two-year test run** before committing to the full period.
Configuration mistakes surface in minutes instead of hours.

Full walkthrough, with a checklist: [`docs/ADAPTING_TO_NEW_AREA.md`](docs/ADAPTING_TO_NEW_AREA.md).
Three ready-made profiles: [`config_examples/`](config_examples/) — coastal
city, temperate city with frost, and state-scale region.

---

## 2. Why this pipeline

Three things distinguish it from the many scripts that compute NDVI from
Landsat:

**Multi-sensor by design.** Most published scripts handle Landsat 8/9 only.
Red and NIR sit on different band numbers in TM/ETM+ (B3/B4) and OLI (B4/B5);
getting that wrong yields a plausible-looking but entirely wrong NDVI. The
mapping is centralised and resolved from the scene identifier.

**Calibration read per scene, not hard-coded.** The Collection 2 scaling
factors are read from each scene's `MTL` file, scoped to the *Level-2*
parameter groups — the same keys also appear in the Level-1 groups with
different values, and a naive parser silently takes the wrong ones. The
published constants are used only as a fallback, and any deviation is logged.

**Generalisation is a first-class feature, not an afterthought.** Every
study-area-specific decision (paths, period, projection, physical bounds, mask
configuration, output scaling) lives in `R/00_config.R`, with a validator that
fails in seconds on a misconfiguration rather than after hours of processing.

The pipeline ships with an automated self-test (33 checks) that runs without
any real data.

---

## 3. Requirements

**R ≥ 4.1** (the native pipe `|>` is used throughout). Tested on R 4.3.3.

| Role | Packages |
|---|---|
| Required | `terra`, `sf`, `dplyr`, `lubridate`, `readr`, `jsonlite`, `stringr`, `tibble`, `tidyr` |
| Figures (optional) | `ggplot2`, `scales`, `ggrepel` |
| Diagnostics (optional) | `knitr`, `pROC` |

If the plotting packages are missing, **tables are still produced** — only the
figures are skipped.

> **Note on dependencies.** This pipeline deliberately does *not* use
> `tidyterra` or `ggspatial`. Both require recent R versions and are often
> unavailable on cluster installations, where their absence would silently
> disable every figure. Module 08 draws rasters, the north arrow and the scale
> bar with base `ggplot2` primitives instead.

`terra` and `sf` need system GDAL/PROJ/GEOS. On Debian/Ubuntu:

```bash
sudo apt install libgdal-dev libproj-dev libgeos-dev libudunits2-dev
```

**Hardware.** A city-sized AOI over 15 years (~200 scenes) runs in 1–3 hours on
a laptop with 8 GB RAM. Requirements scale roughly linearly with scene count;
see [§8](#8-adapting-to-a-new-study-area) for large-region guidance.

---

## 4. Quick start

```bash
git clone <your-repo-url> landsat-c2-zonal
cd landsat-c2-zonal

# 1. Install packages (once)
Rscript R/00_install_packages.R

# 2. Confirm everything works — no real data needed
Rscript tests/run_tests.R

# 3. Edit the configuration — the only file you need to touch
$EDITOR R/00_config.R

# 4. Run
Rscript R/06_main.R
```

---

## 5. Verifying the installation

```bash
Rscript tests/run_tests.R
```

The self-test generates synthetic but **structurally faithful** Landsat scenes
— correct scene-ID grammar, correct band names per sensor, real Collection 2
`QA_PIXEL` bit packing, real `MTL` group structure with deliberately different
Level-1 and Level-2 values — then runs the entire pipeline over them and checks
the results against independently recomputed values.

It verifies, among other things, that:

- the MTL parser takes the **Level-2** scaling factors, not the Level-1 decoys
- the per-pixel weighted mean uses a **per-pixel denominator** (a pixel valid
  in 1 of 3 scenes keeps its value instead of being diluted to a third of it)
- the thermal sentinel-DN threshold **tracks non-standard calibration factors**
- zonal means match an independent recomputation to 1e-7
- the tier filter, the cloud filter and the physical bounds all fire correctly
- the monthly grid is complete, with explicit `NA` for months with no scene

Expected ending:

```
  RESULT: 33 passed, 0 failed
```

Run it again after upgrading R, `terra` or `sf`, and after editing any module.

---

## 6. Input data layout

### Landsat scenes

Download **Collection 2 Level-2 Science Products** from
[USGS EarthExplorer](https://earthexplorer.usgs.gov/) and unpack each scene
into its own folder, exactly as delivered:

```
landsat_data/
└── l08_l09/
    ├── LC08_L2SP_219076_20200115_20200823_02_T1/
    │   ├── ..._SR_B4.TIF          Red
    │   ├── ..._SR_B5.TIF          NIR
    │   ├── ..._ST_B10.TIF         Thermal (surface temperature)
    │   ├── ..._QA_PIXEL.TIF       Quality flags
    │   └── ..._MTL.txt            Metadata (calibration factors)
    └── LC09_L2SP_219076_20200220_20200901_02_T1/
```

Unpacking a folder of tarballs:

```bash
for f in *.tar; do d="${f%.tar}"; mkdir -p "$d" && tar -xf "$f" -C "$d"; done
```

Required bands: **Red**, **NIR**, **ST** (thermal) and **QA_PIXEL**. Scenes
missing any are reported as incomplete and skipped.

The `l04_l05` / `l07` / `l08_l09` grouping is optional — the sensor is detected
from the scene ID, so a single flat folder works:

```r
DIRS_LANDSAT <- list(all = file.path(BASE_DIR, "landsat_data"))
```

**Most common mistake:** an extra nesting level. The pipeline expects exactly
`<folder>/<SCENE_ID>/<band files>`.

### Area of interest

One vector file (`.shp`, `.gpkg`, `.geojson` — anything `sf::st_read` opens),
serving two roles:

| Role | Derived how | Used for |
|---|---|---|
| Outer boundary | all polygons dissolved | clipping rasters |
| Zones | each polygon separately | per-zone statistics |

Requirements: **a defined CRS** (the pipeline reprojects but cannot guess a
missing one), **a name column** for the zones, and **non-overlapping polygons**
(overlaps double-count pixels).

GeoPackage is preferable to Shapefile: one file, no 10-character field-name
limit, no encoding surprises with accented names.

---

## 7. Configuration reference

Everything lives in `R/00_config.R`. Blocks you may want to change are marked
`### EDIT ###`; everything below the "DO NOT EDIT" line is USGS specification
or internal machinery.

| Setting | Meaning | Typical |
|---|---|---|
| `AOI_NAME` | Label for titles and filenames | `"Santos-SP"` |
| `BASE_DIR` | Project root | absolute path |
| `DIRS_LANDSAT` | Folder(s) of scene subfolders | see §5 |
| `SHP_FILE` | AOI polygon file | `.shp` / `.gpkg` |
| `ZONE_NAME_COLUMN` | Attribute holding zone names | `"NM_BAIRRO"` |
| `ZONE_LABEL` | Word for one polygon | `"neighbourhood"` |
| `YEAR_START` / `YEAR_END` | Period, inclusive | `2010` / `2025` |
| `MONTHS_KEEP` | Seasonal filter | `NULL` or `c(12,1,2)` |
| `TIER_KEEP` | Collection tier | `"T1"` |
| `CRS_TARGET` | Projected CRS for all outputs | `"EPSG:31983"` |
| `LST_MIN` / `LST_MAX` | Physical LST bounds (°C) | climate-dependent |
| `QA_BITS_MASK` | QA_PIXEL bits to discard | `c(0,1,3,4,5,7)` |
| `MIN_AOI_OVERLAP_FRAC` | Min scene footprint over AOI | `0.10` (city) |
| `COMPUTE_PIXEL_TABLES` | Per-pixel CSVs | `FALSE` unless small |

Three ready-made profiles are in `config_examples/`.

---

## 8. Adapting to a new study area

Full walkthrough in [`docs/ADAPTING_TO_NEW_AREA.md`](docs/ADAPTING_TO_NEW_AREA.md).
The five settings that actually matter:

**1. `CRS_TARGET`** — must be **projected** (metres), or areas and scale bars
are meaningless.

| Region | CRS |
|---|---|
| Brazil / South America | `EPSG:319xx` (SIRGAS 2000 UTM) |
| Anywhere (global) | `EPSG:326xx` N / `EPSG:327xx` S (WGS 84 UTM) |
| CONUS | `EPSG:5070` |
| Europe | `EPSG:3035` |
| Multi-UTM-zone AOI | equal-area, e.g. `ESRI:102033` |

**2. `LST_MIN` / `LST_MAX`** — climate-dependent, and the most common source of
silently wrong results.

| Climate | `LST_MIN` | `LST_MAX` |
|---|---|---|
| Tropical/subtropical coastal | 5 | 70 |
| Tropical inland / semi-arid | 5 | 75 |
| Temperate (with frost) | −15 | 60 |
| Boreal / high latitude | −45 | 45 |
| Desert | 0 | 80 |

Validate empirically after the first run: inspect `tables/anomalies/`. Pixels
piling up *at* a bound means the bound is cutting into real data. Module 09
also warns at runtime if more than half the scenes hit a bound.

**3. `QA_BITS_MASK`, bit 7 (water)** — masked by default, which is right for
land-surface studies. Remove it if you study water bodies. Watch for zones
dominated by water: they can end up with almost no valid pixels.

**4. `MIN_AOI_OVERLAP_FRAC`** — **the setting that breaks on large AOIs.** One
Landsat scene is ~185 × 180 km. For a city it covers 100% of the AOI; for a
state each scene covers a slice, so the default `0.10` rejects *every* scene
and the run finishes empty.

| AOI size | Value |
|---|---|
| Smaller than one scene | `0.10` |
| A few scenes | `0.01` |
| Many scenes (state, country) | `0.0001` or `0` |

**5. `COMPUTE_PIXEL_TABLES`** — one CSV row per 30 m pixel per period. A 40 km²
city is ~44,000 rows/month; a 250,000 km² state is ~280 **million**. Keep
`FALSE` above medium-city scale. Zonal statistics scale with polygon count, not
area, and are always safe.

---

## 9. Module reference

| Module | Role |
|---|---|
| `00_config.R` | All configuration, validation, logging, package loading |
| `00_install_packages.R` | One-off package installation (run first) |
| `01_scenes.R` | Scene discovery, MTL parsing, inventory |
| `02_qa_mask.R` | QA_PIXEL bit decoding and mask application |
| `03_calc_indices.R` | Radiometric calibration, NDVI, LST |
| `04_clip_aoi.R` | AOI loading, reprojection, clipping, GeoTIFF export |
| `05_export_stats.R` | AOI-wide statistics, temporal aggregation, CSV/TXT |
| `06_main.R` | **Orchestrator — run this one** |
| `07_zonal_stats.R` | Per-zone and per-pixel statistics |
| `08_plots.R` | Maps and charts |
| `09_anomaly_tracker.R` | Quality-control report |
| `99_diagnose_scene.R` | Deep inspection of a single scene |
| `tests/run_tests.R` | End-to-end self-test |

---

## 10. Outputs

```
results_landsat/
├── rasters/{ndvi,lst}/          per-scene clipped GeoTIFFs
├── tables/
│   ├── scenes/                  scene_inventory.csv, ndvi_lst_per_scene.csv
│   ├── monthly/                 ndvi_lst_monthly.csv, ndvi_lst_annual.csv
│   ├── spatial/zones/           per-zone monthly and annual CSV + TXT
│   ├── spatial/pixels/          per-pixel CSV (optional)
│   ├── anomalies/               quality-control reports
│   ├── full_report.txt
│   └── annual_summary.txt
├── plots/
│   ├── pixel/{monthly,annual}/  30 m maps
│   ├── zones/{monthly,annual}/  choropleths, heatmaps, annual series
│   └── timeseries/              AOI-wide dual-panel series
└── logs/processing.log
```

All figures are 300 dpi with English labels, a north arrow, a graphic scale bar
and a provenance caption (`Landsat Collection 2 Level-2 (USGS) | CRS | date`).

---

## 11. Output column dictionary

Statistics tables keep Portuguese column names (`ano`, `mes`, `data_aq`) for
compatibility with existing downstream analyses.

| Column | Meaning |
|---|---|
| `scene_id` | USGS scene identifier |
| `sensor` | `Landsat4/5/7/8/9` |
| `data_aq`, `ano`, `mes` | Acquisition date, year, month |
| `ndvi_media`, `ndvi_mediana` | Mean / median NDVI |
| `ndvi_min`, `ndvi_max`, `ndvi_dp` | Min, max, standard deviation |
| `ndvi_n` | **Valid pixel count — also the aggregation weight** |
| `lst_c_*` | Same set for LST in °C |
| `pct_valido` | % of scene surviving QA masking |
| `cloud_cover` | Scene cloud cover from the MTL |
| `n_cenas` | Scenes contributing to a period |
| `sensores_usados` | e.g. `"L7x2, L8x1"` |
| `zone_id`, `zone_name` | Zone identifier and label |

**Always check `*_n` before interpreting a mean.** A mean over 3 pixels is not
comparable to one over 30,000. Monthly and annual `_dp` columns are weighted
means of per-scene standard deviations, *not* a correct pooled SD — treat them
as a rough heterogeneity indicator only.

---

## 12. Methods


### Radiometric calibration

Collection 2 Level-2 products are distributed as unsigned 16-bit scaled
integers with a fill value of zero. Physical values are recovered by
multiplying each digital number by a scale factor and adding an offset, in that
order [[1]](#ref1):

```
Surface reflectance  ρ     = DN × 0.0000275 + (−0.2)
Surface temperature  Ts(K) = DN × 0.00341802 + 149.0
                     LST(°C) = Ts(K) − 273.15
```

Factors are read per scene from the MTL, scoped to the Level-2 parameter
groups; the constants above are the fallback.

### NDVI

```
NDVI = (NIR − Red) / (NIR + Red)
```

Rouse et al. (1974) [[4]](#ref4); see Tucker (1979) [[5]](#ref5). Band mapping:
Red = B3, NIR = B4 for Landsat 4–7; Red = B4, NIR = B5 for Landsat 8–9.

### LST

The **USGS Level-2 surface temperature product** is used directly, rather than
deriving LST from brightness temperature. That product applies a single-channel
algorithm with ASTER emissivity and atmospheric reanalysis
[[6]](#ref6) [[7]](#ref7). **No additional emissivity or atmospheric correction
is applied**, since doing so on an already corrected product would
double-correct.

Reported accuracy is condition- and surface-dependent rather than a single
figure, and the USGS characterises the product as having reached a
*provisional* level of maturity. Analyses based on within-month spatial
contrasts are robust to a bias common to all pixels of an acquisition.

### Cloud and quality masking

`QA_PIXEL` is produced by **CFMask**, the C implementation of Fmask developed
at Boston University and translated at USGS EROS for operational use
[[8]](#ref8), evaluated against alternative algorithms by Foga et al. (2017)
[[11]](#ref11).

| Bit | Flag | Masked by default |
|---|---|---|
| 0 | Fill | yes |
| 1 | Dilated cloud | yes |
| 2 | Cirrus (L8/9 only) | **no** — see below |
| 3 | Cloud | yes |
| 4 | Cloud shadow | yes |
| 5 | Snow / ice | yes |
| 6 | Clear | never (informational) |
| 7 | Water | yes |

Cirrus is deliberately **not** masked: only L8/9 can report it, so masking it
would make those sensors stricter than L4/5/7 and introduce a sensor-dependent
discontinuity in a long series [[10]](#ref10).

### Temporal aggregation

Scenes within a period are combined by a **mean weighted by valid pixel
count**:

```
weighted_mean = Σ(value_i × n_i) / Σ(n_i)
```

A scene contributing 40,000 clean pixels counts more than one contributing 400.
Months with no usable scene appear as explicit `NA` rows.

For per-pixel composites the denominator accumulates **per pixel**, counting
only scenes with a valid value at that location — otherwise every cell observed
in a subset of scenes would be biased downward.

---

## 13. Quality control cascade

Seven filters, in order, all configurable in `00_config.R` §8:

1. **Scene completeness** — all four required bands present.
2. **QA_PIXEL mask** — per-pixel cloud/shadow/fill/snow/water.
3. **Scene-level cloud** — `MIN_SCENE_VALID_FRAC` of the scene must survive.
4. **AOI footprint overlap** — `MIN_AOI_OVERLAP_FRAC`.
5. **AOI-level validity** — `MIN_AOI_VALID_FRAC` *after* clipping; a scene can
   be clear overall and clouded over your AOI.
6. **NDVI/LST consistency** — `MAX_NDVI_LST_RATIO`; a large excess of valid
   NDVI over valid LST means a degraded thermal band.
7. **Physical bounds** — NDVI ∈ [−1, 1], LST ∈ [`LST_MIN`, `LST_MAX`].
   Out-of-range pixels become `NA`: **discarded, never saturated to the bound**
   (saturating would drag every downstream mean toward the limit).

Everything discarded is recorded by module 09.

---

## 14. Running on an HPC cluster

Set `R_LIBS_USER` in both `00_config.R` and `00_install_packages.R`, then
install once from a node with internet:

```bash
singularity exec --bind $HOME:$HOME container.sif \
  Rscript R/00_install_packages.R
```

SLURM template in [`docs/slurm_example.sh`](docs/slurm_example.sh).

Notes: the pipeline resolves its own script directory, so it can be invoked
from anywhere; lower `TERRA_MEMFRAC` (e.g. `0.4`) on shared nodes; disable
`GENERATE_PLOTS` if the node lacks graphics libraries — tables are unaffected.

---

## 15. Troubleshooting

| Symptom | Cause | Fix |
|---|---|---|
| `No scenes found!` | Wrong `DIRS_LANDSAT`, or folders nested one level deeper | Expected: `<folder>/<SCENE_ID>/<bands>` |
| Every scene "no overlap" | `MIN_AOI_OVERLAP_FRAC` too high for a large AOI | Lower to `0.0001` or `0` |
| Run ends with no results | Usually the above, or wrong `LST_MIN/MAX` | Check `logs/processing.log` |
| All LST discarded | Climate bounds wrong | Widen them; run `99_diagnose_scene.R` |
| Zones all `NA` | `ZONE_NAME_COLUMN` wrong, or CRS mismatch | Set it explicitly; verify AOI has a CRS |
| Coastal zones empty | Water bit masking nearly everything | Check `*_n`; consider removing bit 7 |
| Out of memory | Per-pixel tables on a large AOI | `COMPUTE_PIXEL_TABLES <- FALSE` |
| `MTL unreadable` warnings | Missing MTL files | Harmless — defaults used |
| Maps missing, tables fine | `ggplot2`/`scales` unavailable | Expected degradation |
| Maps look nearly uniform | Colour limits wrong for your climate | Set `PLOT_LIMITS_* <- NULL` for automatic limits |

To inspect one problematic scene end-to-end:

```bash
Rscript R/99_diagnose_scene.R /path/to/SCENE_FOLDER
```

It reports the raw DN range, what those DNs mean in °C, the per-bit QA
breakdown, and whether your configured bounds are binding.

---

## 16. Known limitations

- **No mosaicking.** Scenes are processed independently and aggregated
  statistically. Correct for zonal statistics, but per-pixel maps of a
  multi-scene AOI can show scene-boundary seams.
- **CFMask commission errors over bright targets.** The USGS documents reduced
  performance over building tops and beaches — precisely the land cover of
  interest in coastal urban studies.
- **Thin cirrus** is intrinsically hard to detect and may survive masking;
  the physical range tests are the second line of defence.
- **Snow/ice flag not rigorously validated** by USGS; irrelevant in the
  tropics, relevant at higher latitudes.
- **Landsat 7 SLC-off.** After 2003-05-31, ~22% data gaps; masked as fill, but
  reduces effective sample size.
- **Sensor transitions.** LEDAPS (L4–7) and LaSRC (L8–9) are different
  algorithms; no harmonisation (e.g. Roy et al. 2016 [[12]](#ref12)) is
  applied. The `sensores_usados` column lets you test for sensor effects.
- **Overpass time.** Landsat crosses at ~10:00–10:30 local solar time. LST is a
  mid-morning surface measurement, not a daily mean, and not air temperature
  [[9]](#ref9).
- **Clear-sky bias.** Only cloud-free pixels are retained.
- **Weighted SD.** Monthly/annual `_dp` columns are not a correct pooled SD.

---

## 17. References

<a name="ref1"></a>**[1]** U.S. Geological Survey. *How do I use a scale factor
with Landsat Level-2 science products?*
https://www.usgs.gov/faqs/how-do-i-use-a-scale-factor-landsat-level-2-science-products

<a name="ref2"></a>**[2]** U.S. Geological Survey. *Landsat Collection 2
Level-2 Science Products.*
https://www.usgs.gov/landsat-missions/landsat-collection-2-level-2-science-products

<a name="ref3"></a>**[3]** U.S. Geological Survey. *Landsat 8-9 Collection 2
Level-2 Science Product Guide* (LSDS-1619) and *Landsat 4-7 Collection 2
Level-2 Science Product Guide* (LSDS-1618). Note that the two sensor families
are documented separately.

<a name="ref4"></a>**[4]** Rouse, J.W., Haas, R.H., Schell, J.A. & Deering,
D.W. (1974). Monitoring vegetation systems in the Great Plains with ERTS.
*Third ERTS Symposium*, NASA SP-351, 309–317.

<a name="ref5"></a>**[5]** Tucker, C.J. (1979). Red and photographic infrared
linear combinations for monitoring vegetation. *Remote Sensing of Environment*,
8(2), 127–150. https://doi.org/10.1016/0034-4257(79)90013-0

<a name="ref6"></a>**[6]** Malakar, N.K. et al. (2018). An operational land
surface temperature product for Landsat thermal data. *IEEE TGRS*, 56(10),
5717–5735. https://doi.org/10.1109/TGRS.2018.2824828

<a name="ref7"></a>**[7]** Cook, M. et al. (2014). Development of an
operational calibration methodology for the Landsat thermal data archive.
*Remote Sensing*, 6(11), 11244–11266. https://doi.org/10.3390/rs61111244

<a name="ref8"></a>**[8]** U.S. Geological Survey. *CFMask Algorithm.*
https://www.usgs.gov/landsat-missions/cfmask-algorithm

<a name="ref9"></a>**[9]** Voogt, J.A. & Oke, T.R. (2003). Thermal remote
sensing of urban climates. *Remote Sensing of Environment*, 86(3), 370–384.
https://doi.org/10.1016/S0034-4257(03)00079-8

<a name="ref10"></a>**[10]** Zhu, Z., Wang, S. & Woodcock, C.E. (2015).
Improvement and expansion of the Fmask algorithm. *Remote Sensing of
Environment*, 159, 269–277. https://doi.org/10.1016/j.rse.2014.12.014
See also Zhu, Z. & Woodcock, C.E. (2012), *RSE* 118, 83–94,
https://doi.org/10.1016/j.rse.2011.10.028

<a name="ref11"></a>**[11]** Foga, S. et al. (2017). Cloud detection algorithm
comparison and validation for operational Landsat data products. *Remote
Sensing of Environment*, 194, 379–390.
https://doi.org/10.1016/j.rse.2017.03.026

<a name="ref12"></a>**[12]** Roy, D.P. et al. (2016). Characterization of
Landsat-7 to Landsat-8 reflective wavelength and NDVI continuity. *Remote
Sensing of Environment*, 185, 57–70.
https://doi.org/10.1016/j.rse.2015.12.024

<a name="ref13"></a>**[13]** Hijmans, R.J. *terra: Spatial Data Analysis.*
R package. https://CRAN.R-project.org/package=terra

<a name="ref14"></a>**[14]** Pebesma, E. (2018). Simple Features for R.
*The R Journal*, 10(1), 439–446. https://doi.org/10.32614/RJ-2018-009

---

## 18. Citation and licence

If this pipeline supports a publication, cite the Landsat datasets by DOI (see
`docs/references.bib`), the ST algorithm [[6]](#ref6), the cloud mask
[[11]](#ref11), and this repository as:

```
@software{landsat-nvdi-lst,
  author = {Vian, Gabriel},
  title = {landsat-ndvi-lst},
  year = {2026},
  publisher = {GitHub},
  journal = {GitHub repository},
  howpublished = {\url{https://github.com/G-Vian/landsat-ndvi-lst/tree/main}} 
}
```

Landsat data are courtesy of the U.S. Geological Survey and are in the public
domain. Acknowledge them per the
[USGS data citation policy](https://www.usgs.gov/centers/eros/data-citation).

  
---

## Repository contents

```
landsat-c2-zonal/
├── README.md · LICENSE · .gitignore
├── R/
│   ├── 00_config.R                 ** the only file you edit **
│   ├── 00_install_packages.R       run once, first
│   ├── 01_scenes.R … 09_anomaly_tracker.R
│   ├── 06_main.R                   ** run this **
│   └── 99_diagnose_scene.R         single-scene deep inspection
├── config_examples/                coastal city · temperate · large region
├── tests/
│   ├── run_tests.R                 33-check end-to-end self-test
│   └── make_test_data.R            synthetic scene generator
└── docs/
    ├── ADAPTING_TO_NEW_AREA.md     step-by-step adaptation guide
    ├── CODE_REVIEW.md              review findings, fixes and verification
    ├── methods_section.tex         ready-to-use Methods section
    ├── references.bib              BibTeX for the above
    └── slurm_example.sh            HPC submission template
```
