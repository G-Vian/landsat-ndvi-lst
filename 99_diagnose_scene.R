# =============================================================================
# 99_diagnose_scene.R - SINGLE-SCENE DIAGNOSTIC TOOL
# =============================================================================
# (This replaces the original diagnostico_lst.R, which was hard-coded to one
#  Santos scene and one absolute path.)
#
# Inspects ONE scene in detail, from raw DN through to final LST/NDVI, to
# explain WHY a scene produced suspicious values.
#
# Run it when:
#   - the anomaly report (module 09) flags a scene as critical
#   - a scene is silently discarded and you want to know which filter caught it
#   - you are calibrating LST_MIN / LST_MAX for a new region and want to see
#     the actual DN and temperature distribution before choosing bounds
#
# USAGE
#   Pass the scene folder on the command line (nothing to edit):
#       Rscript R/99_diagnose_scene.R /path/to/LC08_L2SP_219076_20200115_..._T1
#
#   Or ### EDIT ### the SCENE_DIR default below and run with no argument.
#   With neither, it diagnoses the first scene it finds in DIRS_LANDSAT.
# =============================================================================

source("00_config.R")
load_packages()

# --- Scene selection ---------------------------------------------------------
# ### EDIT ### optional default scene folder.
SCENE_DIR <- NULL

args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0) SCENE_DIR <- args[1]

if (is.null(SCENE_DIR)) {
  for (d in DIRS_LANDSAT) {
    if (!dir.exists(d)) next
    subs <- list.dirs(d, recursive = FALSE, full.names = TRUE)
    subs <- subs[grepl(SCENE_ID_PATTERN, basename(subs))]
    if (length(subs) > 0) { SCENE_DIR <- subs[1]; break }
  }
}

if (is.null(SCENE_DIR) || !dir.exists(SCENE_DIR))
  stop("No scene to diagnose. Pass a scene folder as an argument, or set ",
       "SCENE_DIR at the top of this file.")

SCENE_ID <- basename(SCENE_DIR)

cat("\n===========================================================\n")
cat("  SINGLE-SCENE DIAGNOSTIC\n")
cat("===========================================================\n")
cat(sprintf("Scene : %s\n", SCENE_ID))
cat(sprintf("Folder: %s\n\n", SCENE_DIR))

# --- 1. Files present --------------------------------------------------------
files <- list.files(SCENE_DIR)
cat(sprintf("--- FILES IN FOLDER (%d) ---\n", length(files)))
cat(paste(" ", head(files, 20), collapse = "\n"), "\n")
if (length(files) > 20) cat(sprintf("  ... and %d more\n", length(files) - 20))
cat("\n")

# --- 2. Identify sensor and resolve band names -------------------------------
prefix <- sub("_.*$", "", SCENE_ID)
info   <- SENSOR_MAP[[prefix]]
if (is.null(info)) stop("Unknown sensor prefix: ", prefix)

cat("--- SENSOR ---\n")
cat(sprintf("  Satellite : %s\n", info$sensor))
cat(sprintf("  Red       : %s\n", info$red))
cat(sprintf("  NIR       : %s\n", info$nir))
cat(sprintf("  Thermal   : %s\n\n", info$lst))

file_st <- file.path(SCENE_DIR, paste0(SCENE_ID, "_", info$lst, ".TIF"))
file_qa <- file.path(SCENE_DIR, paste0(SCENE_ID, "_QA_PIXEL.TIF"))
if (!file.exists(file_st)) stop("Thermal band not found: ", file_st)

# --- 3. Raw thermal DN inspection --------------------------------------------
cat("--- RAW THERMAL BAND (before any filtering) ---\n")
st_raw <- terra::rast(file_st)
cat(sprintf("  Resolution : %.0f x %.0f m\n",
            terra::res(st_raw)[1], terra::res(st_raw)[2]))
cat(sprintf("  Dimensions : %d rows x %d cols (%d cells)\n",
            nrow(st_raw), ncol(st_raw), terra::ncell(st_raw)))
cat(sprintf("  Data type  : %s\n", terra::datatype(st_raw)[1]))
cat(sprintf("  CRS        : %s\n\n", terra::crs(st_raw, describe = TRUE)$name))

dn_min <- terra::global(st_raw, "min", na.rm = TRUE)[[1]]
dn_max <- terra::global(st_raw, "max", na.rm = TRUE)[[1]]
n_zero <- terra::global(st_raw == 0, "sum", na.rm = TRUE)[[1]]

to_c <- function(dn) dn * ST_SCALE + ST_OFFSET - KELVIN_TO_CELSIUS
cat(sprintf("  DN range      : %.0f to %.0f\n", dn_min, dn_max))
cat(sprintf("  -> in Celsius : %.1f to %.1f C  (UNFILTERED - includes fill)\n",
            to_c(dn_min), to_c(dn_max)))
cat(sprintf("  DN == 0 (fill): %.0f cells (%.2f%% of scene)\n\n",
            n_zero, 100 * n_zero / terra::ncell(st_raw)))

# The sentinel-DN cut-offs implied by the configured limits.
dn_lo <- (LST_MIN + KELVIN_TO_CELSIUS - ST_OFFSET) / ST_SCALE
dn_hi <- (LST_MAX + KELVIN_TO_CELSIUS - ST_OFFSET) / ST_SCALE
cat("--- CONFIGURED BOUNDS, TRANSLATED TO DN ---\n")
cat(sprintf("  LST_MIN = %6.1f C  ->  DN >= %8.0f\n", LST_MIN, dn_lo))
cat(sprintf("  LST_MAX = %6.1f C  ->  DN <= %8.0f\n", LST_MAX, dn_hi))
n_below <- terra::global(st_raw < dn_lo & st_raw > 0, "sum", na.rm = TRUE)[[1]]
n_above <- terra::global(st_raw > dn_hi, "sum", na.rm = TRUE)[[1]]
cat(sprintf("  Non-fill cells below LST_MIN : %.0f\n", n_below))
cat(sprintf("  Cells above LST_MAX          : %.0f\n", n_above))
if (n_below > 0 || n_above > 0)
  cat("  NOTE: these cells will be DISCARDED (set to NA), not truncated.\n")
cat("\n")

# --- 4. QA mask, with a per-bit breakdown ------------------------------------
mask <- NULL
if (file.exists(file_qa)) {
  source("02_qa_mask.R")
  cat("--- QA_PIXEL BREAKDOWN ---\n")
  cat(sprintf("  Bits currently masked: %s\n\n", paste(QA_BITS_MASK, collapse = ", ")))

  qa <- terra::rast(file_qa)
  n  <- terra::ncell(qa)
  bit_names <- c("0 Fill", "1 DilatedCloud", "2 Cirrus", "3 Cloud",
                 "4 CloudShadow", "5 Snow/Ice", "6 Clear", "7 Water")
  cat(sprintf("  %-18s %12s %8s   %s\n", "Flag", "Cells", "%", "In mask?"))
  cat("  ", strrep("-", 52), "\n", sep = "")
  for (b in 0:7) {
    cnt <- terra::global((qa %/% (2^b)) %% 2, "sum", na.rm = TRUE)[[1]]
    cat(sprintf("  %-18s %12.0f %7.2f%%   %s\n",
                bit_names[b + 1], cnt, 100 * cnt / n,
                if (b %in% QA_BITS_MASK) "MASKED" else "-"))
  }
  cat("\n")

  mask <- make_qa_mask(file_qa)
  st   <- qa_stats(mask, SCENE_ID)
  cat(sprintf("  => %.2f%% of the scene survives masking.\n", st$pct_valid))
  if (st$pct_valid < MIN_SCENE_VALID_FRAC * 100)
    cat(sprintf("  => BELOW MIN_SCENE_VALID_FRAC (%.0f%%): this scene WOULD BE REJECTED.\n",
                MIN_SCENE_VALID_FRAC * 100))
  cat("\n")
} else {
  cat("QA_PIXEL not found - skipping mask diagnostics.\n\n")
}

# --- 5. Final products after the full pipeline --------------------------------
if (!is.null(mask)) {
  source("03_calc_indices.R")
  cat("--- FINAL PRODUCTS (after masking, calibration and range tests) ---\n")

  lst <- calc_lst(file_st, mask, ST_SCALE, ST_OFFSET)
  lmin <- terra::global(lst, "min",  na.rm = TRUE)[[1]]
  lmax <- terra::global(lst, "max",  na.rm = TRUE)[[1]]
  lmea <- terra::global(lst, "mean", na.rm = TRUE)[[1]]
  nval <- terra::global(!is.na(lst), "sum")[[1]]
  cat(sprintf("  LST  : %.2f to %.2f C | mean %.2f | %.0f valid px (%.1f%%)\n",
              lmin, lmax, lmea, nval, 100 * nval / terra::ncell(lst)))

  if (abs(lmin - LST_MIN) < 0.1)
    cat("  WARNING: minimum sits exactly at LST_MIN - the bound is binding.\n")
  if (abs(lmax - LST_MAX) < 0.1)
    cat("  WARNING: maximum sits exactly at LST_MAX - the bound is binding.\n")

  file_red <- file.path(SCENE_DIR, paste0(SCENE_ID, "_", info$red, ".TIF"))
  file_nir <- file.path(SCENE_DIR, paste0(SCENE_ID, "_", info$nir, ".TIF"))
  if (file.exists(file_red) && file.exists(file_nir)) {
    ndvi <- calc_ndvi(file_red, file_nir, mask, SR_SCALE, SR_OFFSET)
    cat(sprintf("  NDVI : %.4f to %.4f | mean %.4f | %.0f valid px\n",
                terra::global(ndvi, "min",  na.rm = TRUE)[[1]],
                terra::global(ndvi, "max",  na.rm = TRUE)[[1]],
                terra::global(ndvi, "mean", na.rm = TRUE)[[1]],
                terra::global(!is.na(ndvi), "sum")[[1]]))
  }
  cat("\n")
}

cat("===========================================================\n")
cat("  HOW TO READ THIS\n")
cat("===========================================================\n")
cat("  * Large 'DN == 0' count      -> scene border or SLC-off gaps (normal\n")
cat("                                  for Landsat 7 after May 2003).\n")
cat("  * Cells below LST_MIN        -> if MANY, your LST_MIN may be too high\n")
cat("                                  for this climate (00_config.R Sec. 6).\n")
cat("  * Minimum exactly at bound   -> the bound is cutting into real data.\n")
cat("  * High 'Cloud' or 'Water' %  -> explains a low valid-pixel count.\n")
cat("  * Low survival %             -> compare with MIN_SCENE_VALID_FRAC.\n")
cat("===========================================================\n\n")
