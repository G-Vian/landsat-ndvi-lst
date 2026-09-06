# =============================================================================
# MODULE 00 - GLOBAL CONFIGURATION AND PATHS
# Landsat Collection 2 Level-2 | NDVI & LST time-series pipeline
# =============================================================================
#
#  >>> THIS IS THE ONLY FILE YOU NEED TO EDIT TO RUN THE PIPELINE ON A NEW <<<
#  >>> STUDY AREA. Every block you may want to change is marked with:      <<<
#  >>>                                                                     <<<
#  >>>        # ### EDIT ###                                               <<<
#  >>>                                                                     <<<
#  Everything below the "DO NOT EDIT" line near the bottom is machinery.
#
#  QUICK START (new study area):
#    1. Set AOI_NAME, BASE_DIR, SHP_FILE            (Section 4)
#    2. Set YEAR_START / YEAR_END                   (Section 5)
#    3. Set CRS_TARGET to a projected CRS for your area  (Section 5)
#    4. Set LST_MIN / LST_MAX for your climate      (Section 6)
#    5. Review QA_BITS_MASK — especially the WATER bit  (Section 7)
#    6. Review SCALE_FLAGS if your AOI is large     (Section 8)
#    Then run:  Rscript R/06_main.R
#
#  See config_examples/ for ready-made configs (coastal city, inland city,
#  state-level / large AOI) and docs/ADAPTING_TO_NEW_AREA.md for the full guide.
# =============================================================================


# -----------------------------------------------------------------------------
# 1. SCRIPT LOCATION (needed for Singularity / plain Rscript execution)
# -----------------------------------------------------------------------------
# When invoked as:  singularity exec ... Rscript /full/path/06_main.R
# the working directory is wherever you were in the shell, NOT where the script
# lives. This block detects the real script path and moves the wd there, so that
# source("01_scenes.R") works without absolute paths.
#
# You normally do NOT need to touch this block.

# Source guard: 00_config.R is sourced by every module, but the setup below
# only needs to run once per session. Without this, a pipeline run prints the
# same startup banner a dozen times.
if (exists(".CONFIG_LOADED", envir = globalenv()) &&
    isTRUE(get(".CONFIG_LOADED", envir = globalenv()))) {
  # already configured in this session - skip the noisy setup
} else {

SCRIPT_DIR <- tryCatch({
  args     <- commandArgs(trailingOnly = FALSE)
  file_arg <- args[grep("--file=", args)]
  if (length(file_arg) > 0) {
    dirname(normalizePath(sub("--file=", "", file_arg[1])))
  } else {
    getwd()   # fallback: interactive session (RStudio)
  }
}, error = function(e) getwd())

if (dir.exists(SCRIPT_DIR) && SCRIPT_DIR != getwd()) setwd(SCRIPT_DIR)
message(sprintf("[00_config] Working directory: %s", getwd()))
}   # end of the once-per-session setup block


# -----------------------------------------------------------------------------
# 2. R PACKAGE LIBRARY
# -----------------------------------------------------------------------------
# ### EDIT ### — only if you run on an HPC cluster with a container.
#
# On a normal laptop, leave R_LIBS_USER = NULL: the pipeline uses R's default
# library and installs missing packages there.
#
# On a cluster where the container filesystem is read-only, point this at a
# writable directory in your home, and install packages ONCE before submitting
# the job (nodes often have no internet):
#
#   R_LIBS_USER <- "~/R_libs"
#   singularity exec --bind $HOME:$HOME container.sif \
#     Rscript R/00_install_packages.R

R_LIBS_USER <- NULL     # e.g. "~/R_libs" on a cluster; NULL = default library

if (!is.null(R_LIBS_USER)) {
  R_LIBS_USER <- path.expand(R_LIBS_USER)
  if (!dir.exists(R_LIBS_USER)) {
    dir.create(R_LIBS_USER, recursive = TRUE)
    message("[00_config] Created package library: ", R_LIBS_USER)
  }
  .libPaths(c(R_LIBS_USER, .libPaths()))
  message(sprintf("[00_config] Active R_libs: %s", R_LIBS_USER))
}

# Mark configuration as loaded (see the source guard at the top of Section 1)
.CONFIG_LOADED <- TRUE


# -----------------------------------------------------------------------------
# 3. REQUIRED PACKAGES
# -----------------------------------------------------------------------------
# Core pipeline only. Plotting packages (module 08) are loaded lazily there,
# so the pipeline still produces tables if the plotting stack is unavailable.

required_packages <- c("terra", "sf", "dplyr", "lubridate",
                       "readr", "jsonlite", "stringr", "tibble", "tidyr")


# -----------------------------------------------------------------------------
# 4. STUDY AREA AND PATHS
# -----------------------------------------------------------------------------
# ### EDIT ### — this is the main block to change for a new area.

# Short label for your area of interest. Used in plot titles, report headers
# and output file names. Keep it short and filename-safe.
AOI_NAME <- "MyStudyArea"

# Project root. Everything else is resolved relative to this.
# Use an absolute path on a cluster; "~/..." works locally.
BASE_DIR <- path.expand("~/landsat_project")

# --- Input: Landsat scene folders -------------------------------------------
# Each entry points to a folder that CONTAINS ONE SUBFOLDER PER SCENE, exactly
# as EarthExplorer / USGS delivers them, e.g.:
#
#   landsat_data/l08_l09/
#     LC08_L2SP_219076_20200115_20200823_02_T1/
#       LC08_L2SP_219076_20200115_20200823_02_T1_SR_B4.TIF
#       LC08_L2SP_219076_20200115_20200823_02_T1_SR_B5.TIF
#       LC08_L2SP_219076_20200115_20200823_02_T1_ST_B10.TIF
#       LC08_L2SP_219076_20200115_20200823_02_T1_QA_PIXEL.TIF
#       LC08_L2SP_219076_20200115_20200823_02_T1_MTL.txt
#     LC09_L2SP_.../
#
# The grouping into separate folders is only for YOUR convenience — the
# pipeline detects the sensor from the scene ID, not from the folder name.
# You may use a single folder for everything:
#   
DIRS_LANDSAT <- list(
  l04_l05 = file.path(BASE_DIR, "landsat_data/l04_l05"),  # Landsat 4 & 5 (TM)
  l07     = file.path(BASE_DIR, "landsat_data/l07"),      # Landsat 7  (ETM+)
  l08_l09 = file.path(BASE_DIR, "landsat_data/l08_l09")   # Landsat 8 & 9 (OLI/TIRS)
)

# --- Input: area-of-interest polygon(s) --------------------------------------
# A single shapefile (or GeoPackage / GeoJSON — anything sf::st_read can open)
# defining your study area.
#
# TWO ROLES, both derived from this one file:
#   (a) OUTER BOUNDARY: all polygons dissolved into one → used to clip rasters.
#   (b) ZONES: each polygon individually → used for per-zone statistics
#       (neighbourhoods, districts, municipalities, census tracts, ...).
#
# REQUIREMENTS:
#   - Must have a valid CRS defined (the pipeline reprojects, but cannot guess).
#   - Should have a column with the zone NAME (see ZONE_NAME_COLUMN below).
#   - Polygons should not overlap (overlapping zones double-count pixels).

SHP_FILE <- file.path(BASE_DIR, "shapefiles/my_area.shp")

# Name of the attribute column holding each zone's label.
# Set to NULL to auto-detect (tries NOME, NAME, BAIRRO, NM_MUN, MUNICIPIO, ...
# then falls back to the first character column). Setting it explicitly is
# strongly recommended — auto-detection can pick the wrong column.
ZONE_NAME_COLUMN <- NULL          # e.g. "NM_BAIRRO", "NAME_2", "district"

# Generic label for what one polygon represents, used in report headers and
# column names in the output tables ("zone", "neighbourhood", "district", ...).
ZONE_LABEL <- "zone"

# --- Output ------------------------------------------------------------------
# Created automatically. Safe to delete and regenerate.
OUTPUT_DIR <- file.path(BASE_DIR, "results_landsat")


# -----------------------------------------------------------------------------
# 5. TIME PERIOD AND PROJECTION
# -----------------------------------------------------------------------------
# ### EDIT ###

# Period of interest (inclusive). Scenes outside this range are skipped at
# inventory time, so restricting the period is the cheapest way to speed up a
# test run. Landsat availability: L4/L5 1982-2013, L7 1999-2024 (SLC-off after
# 2003-05-31), L8 2013-, L9 2021-.
YEAR_START <- 2010
YEAR_END   <- 2025

# OPTIONAL month filter — keep only these calendar months (1-12).
# Useful for seasonal studies ("dry season only") or to cut processing time.
# NULL = all months.
MONTHS_KEEP <- NULL               # e.g. c(12, 1, 2, 3) for austral summer

# OPTIONAL Landsat collection tier filter.
#   "T1"  = Tier 1 only. Highest geometric quality, radiometrically consistent.
#           RECOMMENDED for time-series analysis.
#   "T2"  = Tier 2 only (lower geometric accuracy).
#   NULL  = accept all tiers, including Real-Time ("RT").
# Mixing tiers in a time series can introduce geolocation jitter between dates.
TIER_KEEP <- "T1"

# Target CRS for ALL outputs (rasters, clipping, area computations).
# ### EDIT ### — MUST be changed for a new region.
#
# Use a PROJECTED CRS (metres), not geographic (degrees), so that pixel areas
# and the graphic scale bar are meaningful.
#   - Brazil / South America:  SIRGAS 2000 UTM  → EPSG:319xx (xx = 01..25 zone)
#       e.g. "EPSG:31983" = SIRGAS 2000 / UTM 23S (São Paulo coast)
#   - Global default:          WGS 84 UTM       → EPSG:326xx (N) / 327xx (S)
#       e.g. "EPSG:32723" = WGS 84 / UTM 23S
#   - Continental USA:         "EPSG:5070" (NAD83 Conus Albers)
#   - Europe:                  "EPSG:3035" (ETRS89 LAEA)
#
# Find your UTM zone: floor((longitude + 180) / 6) + 1
# For AOIs spanning several UTM zones, prefer an equal-area CRS (Albers/LAEA).
CRS_TARGET <- "EPSG:32723"


# -----------------------------------------------------------------------------
# 6. PHYSICAL VALIDITY LIMITS
# -----------------------------------------------------------------------------
# ### EDIT ### — LST limits are CLIMATE-DEPENDENT and must be revised.
#
# Pixels outside these ranges are set to NA (DISCARDED, not clipped to the
# boundary value). They are residual artefacts: unmasked cirrus, cloud-shadow
# edges, Landsat-7 SLC-off stripe borders, and fill-value bleed.

# NDVI is mathematically bounded to [-1, 1]; values outside are numeric noise.
# There is normally no reason to change these.
NDVI_MIN <- -1.0
NDVI_MAX <-  1.0

# LST bounds in degrees Celsius. THE DEFAULTS BELOW ARE DELIBERATELY WIDE.
#
# How to choose them for your area:
#   LST_MIN — a few degrees BELOW the coldest plausible *surface* temperature.
#             Surface is not air temperature: it can be several degrees colder
#             than air on clear winter nights (Landsat overpass is ~10:00-10:30
#             local solar time, so this is a daytime measurement).
#   LST_MAX — a few degrees ABOVE the hottest plausible surface: dark asphalt
#             and metal roofs routinely reach 60-70 °C in tropical summer.
#
# Suggested starting points:
#   Tropical / subtropical coastal ...  LST_MIN =   5, LST_MAX = 70
#   Tropical inland / semi-arid ......  LST_MIN =   5, LST_MAX = 75
#   Temperate (with frost) ...........  LST_MIN = -15, LST_MAX = 60
#   Boreal / high latitude ...........  LST_MIN = -45, LST_MAX = 45
#   Desert ...........................  LST_MIN =   0, LST_MAX = 80
#
# HOW TO VALIDATE YOUR CHOICE: run the pipeline once, then inspect
# tables/anomalies/. If many scenes report pixels piling up exactly AT a limit,
# your bound is cutting into real data — widen it. If nothing is ever
# discarded, your bounds may be too loose to catch artefacts — tighten them.
LST_MIN <- -20.0
LST_MAX <-  75.0


# -----------------------------------------------------------------------------
# 7. QA_PIXEL MASKING
# -----------------------------------------------------------------------------
# ### EDIT ### — review the WATER bit for your application.
#
# Collection 2 QA_PIXEL bit layout (16-bit, identical for all Landsat sensors):
#   Bit 0  Fill (no data: scene border, SLC-off gaps)
#   Bit 1  Dilated Cloud (buffer around detected clouds)
#   Bit 2  Cirrus (Landsat 8/9 only; always 0 on L4/5/7)
#   Bit 3  Cloud
#   Bit 4  Cloud Shadow
#   Bit 5  Snow / Ice
#   Bit 6  Clear  (informational — do NOT put this in the mask list!)
#   Bit 7  Water
#   Bits 8-15 confidence levels (see USGS product guide, Table 6-3)
#
# A pixel is discarded if ANY listed bit is set.
#
# Decisions you should make consciously:
#
#   Bit 5 (Snow/Ice) — keep it masked unless snow itself is your subject.
#
#   Bit 7 (WATER) — THE MOST CONSEQUENTIAL CHOICE HERE.
#     Masking water is appropriate when you study the LAND surface: water has
#     very high thermal inertia (its LST says nothing about urban heat) and
#     strongly negative NDVI that would drag zone means down. This is the
#     default and is what you want for urban/vegetation studies.
#     REMOVE bit 7 from this list if you are studying water bodies, reservoirs,
#     flooding, or coastal/estuarine surface temperature.
#     WARNING FOR ZONES DOMINATED BY WATER: a coastal or riverine zone can end
#     up with almost no valid pixels once water is masked. Check the *_n columns
#     in the output tables before trusting such a zone's mean.
#
#   Bit 2 (Cirrus) — NOT masked by default. Adding it makes L8/L9 stricter
#     than L4/5/7 (which cannot report cirrus), introducing a SENSOR-DEPENDENT
#     bias in a long time series. Add it only if you process L8/L9 exclusively.
QA_BITS_MASK <- c(0, 1, 3, 4, 5, 7)


# -----------------------------------------------------------------------------
# 8. QUALITY THRESHOLDS AND SCALING BEHAVIOUR
# -----------------------------------------------------------------------------
# ### EDIT ### — especially MIN_AOI_OVERLAP_FRAC if your AOI is large.

# (a) Minimum fraction of the FULL SCENE that must survive QA masking before
#     the scene is processed at all. Cheap early rejection of cloud-covered
#     scenes. 0.10 = at least 10% of scene pixels usable.
MIN_SCENE_VALID_FRAC <- 0.10

# (b) Minimum fraction of the AOI that must have valid pixels AFTER clipping.
#     Guards against a scene that is clear elsewhere but cloudy over your AOI.
MIN_AOI_VALID_FRAC <- 0.10

# (c) Minimum fraction of the AOI a scene footprint must cover to be used.
#
#     ############################ IMPORTANT ############################
#     THIS IS THE #1 SETTING THAT BREAKS WHEN MOVING TO A LARGE AOI.
#
#     One Landsat scene covers roughly 185 x 180 km. If your AOI is a city,
#     a single scene covers 100% of it and the default 0.10 is fine.
#     If your AOI is a STATE, a PROVINCE or a COUNTRY, each individual scene
#     covers only a small fraction of it, and a 0.10 threshold would reject
#     EVERY SCENE, producing an empty run.
#
#     Rule of thumb:
#       AOI smaller than one scene (city, park, watershed) ..... 0.10
#       AOI of a few scenes (metro region, small state) ........ 0.01
#       AOI of many scenes (large state, country) ............. 0.0001 or 0
#     Setting 0 disables the check (any intersection is accepted).
#     ###################################################################
MIN_AOI_OVERLAP_FRAC <- 0.10

# (d) Maximum tolerated ratio between valid NDVI pixels and valid LST pixels.
#     A much larger NDVI count means the thermal band is degraded for that
#     scene (common in some L7 SLC-off acquisitions). Scene is discarded.
MAX_NDVI_LST_RATIO <- 5

# (e) Maximum fraction of valid LST pixels allowed to sit exactly at LST_MIN or
#     LST_MAX. Above this the thermal band is systematically saturated and the
#     scene is discarded. 0.01 = 1%.
MAX_LST_AT_BOUND_FRAC <- 0.01

# --- Scaling flags: turn OFF expensive products for large AOIs ---------------
# ### EDIT ### — set the pixel-level ones to FALSE for anything bigger than a
# medium city. A per-pixel CSV has ONE ROW PER 30 m PIXEL PER PERIOD: a city of
# 40 km² is ~44,000 rows per month (fine); a 250,000 km² state is ~280 MILLION
# rows per month (will exhaust memory and disk).
SAVE_CLIPPED_RASTERS  <- TRUE   # per-scene GeoTIFFs. Disk-heavy but reusable.
COMPUTE_ZONE_STATS    <- TRUE   # per-zone tables. Scales well; keep TRUE.
COMPUTE_PIXEL_TABLES  <- FALSE  # per-pixel CSVs. ONLY for small AOIs.
GENERATE_PLOTS        <- TRUE   # maps and charts. One PNG per period.

# terra memory budget: fraction of available RAM terra may use before it
# switches to on-disk processing. Lower it (0.3) on a shared cluster node.
TERRA_MEMFRAC <- 0.6


# =============================================================================
# ############################ DO NOT EDIT BELOW ##############################
# The values below are defined by the USGS product specification, or are
# internal machinery. Changing them without consulting the product guide
# (see docs/REFERENCES.md) will silently corrupt your results.
# =============================================================================

# -----------------------------------------------------------------------------
# 9. USGS COLLECTION 2 LEVEL-2 SCALING FACTORS
# -----------------------------------------------------------------------------
# These are FALLBACKS ONLY. The pipeline reads the actual per-scene values from
# each scene's MTL file and uses those; these constants are used only when the
# MTL is missing or unreadable, and any per-scene deviation is logged as a
# warning by verify_scaling() in 01_scenes.R.
#
# Surface Reflectance:        SR  = DN * 0.0000275 + (-0.2)
# Surface Temperature (K):    ST  = DN * 0.00341802 + 149.0
# Source: USGS Landsat Collection 2 Level-2 Science Product Guide (LSDS-1619).

SR_SCALE          <- 0.0000275   # surface reflectance multiplicative factor
SR_OFFSET         <- -0.2        # surface reflectance additive offset
ST_SCALE          <- 0.00341802  # surface temperature multiplicative factor
ST_OFFSET         <- 149.0       # surface temperature additive offset (Kelvin)
KELVIN_TO_CELSIUS <- 273.15


# -----------------------------------------------------------------------------
# 10. SENSOR -> BAND MAPPING
# -----------------------------------------------------------------------------
# Band numbering differs between sensor generations. This table is what lets a
# single code path handle Landsat 4 through 9.
#
#   TM (L4/L5) and ETM+ (L7):  Red = B3, NIR = B4, Thermal = ST_B6
#   OLI/TIRS   (L8/L9):        Red = B4, NIR = B5, Thermal = ST_B10
#
# Getting this wrong silently produces a plausible-looking but WRONG NDVI,
# which is why the mapping is centralised here rather than inlined.

SENSOR_MAP <- list(
  LC08 = list(sensor = "Landsat8", red = "SR_B4", nir = "SR_B5", lst = "ST_B10"),
  LC09 = list(sensor = "Landsat9", red = "SR_B4", nir = "SR_B5", lst = "ST_B10"),
  LE07 = list(sensor = "Landsat7", red = "SR_B3", nir = "SR_B4", lst = "ST_B6"),
  LT05 = list(sensor = "Landsat5", red = "SR_B3", nir = "SR_B4", lst = "ST_B6"),
  LT04 = list(sensor = "Landsat4", red = "SR_B3", nir = "SR_B4", lst = "ST_B6")
)

# Regex matching a valid Collection 2 Level-2 scene directory name.
SCENE_ID_PATTERN <- paste0("^(", paste(names(SENSOR_MAP), collapse = "|"), ")_L2SP_")


# -----------------------------------------------------------------------------
# 11. BACKWARD-COMPATIBILITY ALIASES
# -----------------------------------------------------------------------------
# Older versions of this pipeline used Portuguese / Santos-specific names.
# These aliases let existing downstream scripts keep working unchanged.
ANO_INICIO <- YEAR_START
ANO_FIM    <- YEAR_END
CRS_ALVO   <- CRS_TARGET


# -----------------------------------------------------------------------------
# 12. UTILITY FUNCTIONS
# -----------------------------------------------------------------------------

#' Check, install and load all required packages.
#'
#' On a cluster without internet on the compute nodes, installation will fail;
#' the error message then tells the user to run 00_install_packages.R from a
#' node that has internet. Package versions are always logged, for
#' reproducibility.
load_packages <- function() {

  if (!is.null(R_LIBS_USER) && !R_LIBS_USER %in% .libPaths())
    .libPaths(c(R_LIBS_USER, .libPaths()))

  missing_pkgs <- required_packages[
    !sapply(required_packages, requireNamespace, quietly = TRUE)
  ]

  if (length(missing_pkgs) > 0) {
    target_lib <- if (is.null(R_LIBS_USER)) .libPaths()[1] else R_LIBS_USER
    message(sprintf("[00_config] Missing packages: %s\n  Installing into: %s",
                    paste(missing_pkgs, collapse = ", "), target_lib))

    tryCatch({
      install.packages(missing_pkgs, lib = target_lib,
                       repos = "https://cloud.r-project.org",
                       quiet = FALSE, dependencies = TRUE)
    }, error = function(e) {
      stop(paste(
        "[00_config] FAILED to install packages.",
        "If your cluster nodes have no internet, run this ONCE beforehand",
        "from a node that does:",
        "  Rscript R/00_install_packages.R",
        sprintf("Original error: %s", conditionMessage(e)),
        sep = "\n"
      ))
    })
  }

  invisible(lapply(required_packages, library,
                   character.only = TRUE, lib.loc = .libPaths()))

  # terra memory budget (see TERRA_MEMFRAC above)
  if (requireNamespace("terra", quietly = TRUE))
    try(terra::terraOptions(memfrac = TERRA_MEMFRAC), silent = TRUE)

  versions <- sapply(required_packages,
                     function(p) as.character(packageVersion(p)))
  message("[00_config] Packages loaded:")
  for (p in names(versions)) message(sprintf("  %-12s v%s", p, versions[p]))
}


#' Create the output directory tree (idempotent).
create_output_dirs <- function() {
  dirs <- c(
    OUTPUT_DIR,
    file.path(OUTPUT_DIR, "rasters", "ndvi"),
    file.path(OUTPUT_DIR, "rasters", "lst"),
    file.path(OUTPUT_DIR, "tables", "scenes"),
    file.path(OUTPUT_DIR, "tables", "monthly"),
    file.path(OUTPUT_DIR, "tables", "spatial"),
    file.path(OUTPUT_DIR, "tables", "anomalies"),
    file.path(OUTPUT_DIR, "plots"),
    file.path(OUTPUT_DIR, "logs")
  )
  lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE)
  message("[00_config] Output tree created under: ", OUTPUT_DIR)
}


#' Timestamped logger — writes to console and to logs/processing.log.
log_msg <- function(msg, level = "INFO") {
  ts  <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  txt <- sprintf("[%s][%s] %s", ts, level, msg)
  message(txt)
  log_file <- file.path(OUTPUT_DIR, "logs", "processing.log")
  if (dir.exists(dirname(log_file)))
    cat(txt, "\n", file = log_file, append = TRUE)
}


#' Validate the configuration before any heavy processing starts.
#'
#' Fails fast and loudly on the mistakes that are most common when adapting
#' this pipeline to a new study area, instead of letting the job run for hours
#' and produce an empty result.
validate_config <- function() {

  problems <- character(0)
  warns    <- character(0)

  # --- paths ---
  if (!file.exists(SHP_FILE))
    problems <- c(problems, sprintf("SHP_FILE does not exist: %s", SHP_FILE))

  existing_dirs <- vapply(DIRS_LANDSAT, dir.exists, logical(1))
  if (!any(existing_dirs))
    problems <- c(problems, paste0(
      "None of the DIRS_LANDSAT folders exist:\n    ",
      paste(unlist(DIRS_LANDSAT), collapse = "\n    ")))
  else if (any(!existing_dirs))
    warns <- c(warns, sprintf("Landsat folder(s) not found (skipped): %s",
                              paste(names(DIRS_LANDSAT)[!existing_dirs], collapse = ", ")))

  # --- period ---
  if (YEAR_START > YEAR_END)
    problems <- c(problems, sprintf("YEAR_START (%d) > YEAR_END (%d)", YEAR_START, YEAR_END))
  if (YEAR_END > as.integer(format(Sys.Date(), "%Y")))
    warns <- c(warns, "YEAR_END is in the future — no scenes will exist for those years.")

  # --- CRS ---
  if (!grepl("^EPSG:[0-9]+$", CRS_TARGET))
    warns <- c(warns, sprintf("CRS_TARGET '%s' is not in 'EPSG:nnnnn' form; make sure terra understands it.",
                              CRS_TARGET))

  # --- physical limits ---
  if (LST_MIN >= LST_MAX)
    problems <- c(problems, sprintf("LST_MIN (%.1f) >= LST_MAX (%.1f)", LST_MIN, LST_MAX))
  if (NDVI_MIN < -1 || NDVI_MAX > 1)
    warns <- c(warns, "NDVI limits outside [-1, 1] — NDVI is mathematically bounded there.")

  # --- QA bits ---
  if (6 %in% QA_BITS_MASK)
    problems <- c(problems, paste(
      "QA_BITS_MASK contains bit 6, which is the CLEAR flag.",
      "Masking it would discard every GOOD pixel and keep only bad ones."))
  if (!(0 %in% QA_BITS_MASK))
    warns <- c(warns, "Bit 0 (Fill) is not masked — scene-border no-data will enter the statistics.")
  if (!(7 %in% QA_BITS_MASK))
    warns <- c(warns, "Bit 7 (Water) is NOT masked — water pixels will be included (intentional for aquatic studies).")

  # --- scaling sanity for large AOIs ---
  if (isTRUE(COMPUTE_PIXEL_TABLES))
    warns <- c(warns, paste(
      "COMPUTE_PIXEL_TABLES = TRUE produces one CSV row per 30 m pixel per period.",
      "Verify this is tractable for your AOI size before a long run."))

  if (length(warns) > 0) {
    message("\n[00_config] CONFIGURATION WARNINGS:")
    for (w in warns) message("  ! ", w)
  }

  if (length(problems) > 0) {
    stop(paste0("\n[00_config] CONFIGURATION ERRORS:\n  - ",
                paste(problems, collapse = "\n  - "),
                "\n\nFix R/00_config.R and re-run. ",
                "See docs/ADAPTING_TO_NEW_AREA.md for guidance.\n"))
  }

  message("[00_config] Configuration validated.")
  invisible(TRUE)
}
