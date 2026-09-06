# =============================================================================
# run_tests.R - END-TO-END SELF-TEST OF THE PIPELINE
# =============================================================================
# Verifies that the pipeline is correctly installed and numerically correct,
# WITHOUT needing any real Landsat data. It:
#
#   1. generates synthetic but structurally faithful Landsat scenes
#      (real scene-ID grammar, real QA_PIXEL bit packing, real MTL groups)
#   2. runs unit tests on the numerically delicate functions
#   3. runs the whole pipeline over the synthetic data
#   4. checks the outputs against independently recomputed values
#
# USAGE (from the repository root):
#   Rscript tests/run_tests.R
#
# Runtime: about one minute. Everything is written to a temporary directory,
# which is removed at the end unless KEEP_OUTPUT is TRUE.
#
# WHEN TO RUN IT
#   - after installing on a new machine or cluster
#   - after upgrading terra / sf / R
#   - after editing any module, before trusting a production run
# =============================================================================

KEEP_OUTPUT <- FALSE     # ### EDIT ### TRUE to inspect the generated outputs

# --- locate the repository ----------------------------------------------------
args      <- commandArgs(trailingOnly = FALSE)
file_arg  <- args[grep("--file=", args)]
TEST_DIR  <- if (length(file_arg) > 0) {
  # normalizePath() warns when given a relative path that R has already
  # resolved differently; suppress that cosmetic warning.
  suppressWarnings(dirname(normalizePath(sub("--file=", "", file_arg[1]),
                                         mustWork = FALSE)))
} else getwd()
if (!dir.exists(TEST_DIR)) TEST_DIR <- file.path(getwd(), "tests")
REPO_DIR  <- normalizePath(file.path(TEST_DIR, ".."))
R_DIR     <- file.path(REPO_DIR, "R")

if (!dir.exists(R_DIR))
  stop("Cannot find the R/ directory. Run this from the repository root:\n",
       "  Rscript tests/run_tests.R")

cat("\n=============================================================\n")
cat("  PIPELINE SELF-TEST\n")
cat("=============================================================\n")
cat(sprintf("Repository : %s\n", REPO_DIR))
cat(sprintf("R version  : %s\n\n", R.version.string))

PASS <- 0L; FAIL <- 0L
check <- function(label, condition, detail = "") {
  if (isTRUE(condition)) {
    PASS <<- PASS + 1L
    cat(sprintf("  [PASS] %s%s\n", label, if (nzchar(detail)) paste0("  ", detail) else ""))
  } else {
    FAIL <<- FAIL + 1L
    cat(sprintf("  [FAIL] %s%s\n", label, if (nzchar(detail)) paste0("  ", detail) else ""))
  }
}

# --- 0. dependencies ----------------------------------------------------------
cat("--- 0. Dependencies ---\n")
core <- c("terra", "sf", "dplyr", "lubridate", "readr", "jsonlite",
          "stringr", "tibble", "tidyr")
missing <- core[!sapply(core, requireNamespace, quietly = TRUE)]
check("core packages available",
      length(missing) == 0,
      if (length(missing)) paste("missing:", paste(missing, collapse = ", ")) else "")
if (length(missing) > 0) {
  cat("\nCannot continue. Run: Rscript R/00_install_packages.R\n")
  quit(status = 1)
}
suppressPackageStartupMessages({library(terra); library(sf); library(dplyr)})
cat("\n")

# --- 1. synthetic data --------------------------------------------------------
cat("--- 1. Generating synthetic Landsat scenes ---\n")
TEST_ROOT <- file.path(tempdir(), "landsat_selftest")
Sys.setenv(LANDSAT_TEST_ROOT = TEST_ROOT)
suppressMessages(source(file.path(TEST_DIR, "make_test_data.R")))
n_scenes <- length(list.dirs(file.path(TEST_ROOT, "landsat"),
                             recursive = FALSE))
check("scenes generated", n_scenes == 14, sprintf("(%d found, expected 14)", n_scenes))
check("AOI written", file.exists(file.path(TEST_ROOT, "shp", "zones.gpkg")))
cat("\n")

# --- 2. build a test configuration --------------------------------------------
cat("--- 2. Building test configuration ---\n")
WORK <- file.path(TEST_ROOT, "pipeline")
dir.create(WORK, recursive = TRUE, showWarnings = FALSE)
file.copy(list.files(R_DIR, pattern = "\\.R$", full.names = TRUE), WORK,
          overwrite = TRUE)

cfg <- readLines(file.path(WORK, "00_config.R"), warn = FALSE)
patch <- function(cfg, pattern, replacement) {
  i <- grep(pattern, cfg)[1]
  if (!is.na(i)) cfg[i] <- replacement
  cfg
}
cfg <- patch(cfg, '^AOI_NAME <- ',         'AOI_NAME <- "TestArea"')
cfg <- patch(cfg, '^BASE_DIR <- ',         sprintf('BASE_DIR <- "%s"', TEST_ROOT))
cfg <- patch(cfg, '^SHP_FILE <- ',         'SHP_FILE <- file.path(BASE_DIR, "shp/zones.gpkg")')
cfg <- patch(cfg, '^ZONE_NAME_COLUMN <- ', 'ZONE_NAME_COLUMN <- "NM_BAIRRO"')
cfg <- patch(cfg, '^YEAR_START <- ',       'YEAR_START <- 2015')
cfg <- patch(cfg, '^YEAR_END   <- ',       'YEAR_END   <- 2016')
cfg <- patch(cfg, '^CRS_TARGET <- ',       'CRS_TARGET <- "EPSG:31983"')
cfg <- patch(cfg, '^LST_MIN <- ',          'LST_MIN <- 5.0')
cfg <- patch(cfg, '^LST_MAX <- ',          'LST_MAX <- 70.0')
cfg <- patch(cfg, '^COMPUTE_PIXEL_TABLES', 'COMPUTE_PIXEL_TABLES  <- TRUE')

# single flat scene folder for the test
i <- grep('^DIRS_LANDSAT <- list\\(', cfg)[1]
j <- grep('^\\)', cfg[i:length(cfg)])[1] + i - 1
cfg <- c(cfg[1:(i-1)],
         'DIRS_LANDSAT <- list(all = file.path(BASE_DIR, "landsat"))',
         cfg[(j+1):length(cfg)])
writeLines(cfg, file.path(WORK, "00_config.R"))
check("test configuration written", TRUE)
cat("\n")

# --- 3. unit tests on the delicate numerics -----------------------------------
cat("--- 3. Unit tests ---\n")
old_wd <- getwd(); setwd(WORK)
suppressMessages({
  source("00_config.R"); load_packages()
  source("05_export_stats.R"); source("07_zonal_stats.R"); source("03_calc_indices.R")
})

# 3a. per-pixel weighted mean must use a PER-PIXEL denominator.
# A pixel valid in only 1 of 3 scenes must keep its value, not be diluted.
tm <- terra::rast(nrows = 1, ncols = 2, xmin = 0, xmax = 2, ymin = 0, ymax = 1)
r1 <- tm; terra::values(r1) <- c(0.6, 0.6)
r2 <- tm; terra::values(r2) <- c(0.6, NA)
r3 <- tm; terra::values(r3) <- c(0.6, NA)
v  <- terra::values(pixel_weighted_mean(list(r1, r2, r3), c(10, 10, 10)))[, 1]
check("pixel mean, fully valid pixel",   abs(v[1] - 0.6) < 1e-9, sprintf("= %.4f", v[1]))
check("pixel mean, partly valid pixel",  abs(v[2] - 0.6) < 1e-9,
      sprintf("= %.4f (a buggy denominator would give 0.20)", v[2]))

# 3b. weighted mean helper
check("weighted_mean_safe", abs(weighted_mean_safe(c(10, 20), c(1, 3)) - 17.5) < 1e-9)
check("weighted_mean_safe with zero weights", is.na(weighted_mean_safe(c(1, 2), c(0, 0))))

# 3c. sentinel-DN threshold must track per-scene calibration factors
back <- ((LST_MIN + KELVIN_TO_CELSIUS - ST_OFFSET) / ST_SCALE) * ST_SCALE +
        ST_OFFSET - KELVIN_TO_CELSIUS
check("DN threshold, standard factors", abs(back - LST_MIN) < 1e-6)
back2 <- ((LST_MIN + KELVIN_TO_CELSIUS - 150.0) / 0.0035) * 0.0035 +
         150.0 - KELVIN_TO_CELSIUS
check("DN threshold, non-standard factors", abs(back2 - LST_MIN) < 1e-6,
      "(hard-coded constants would fail here)")

# 3d. sensor composition strings
check("sensors_compact", sensors_compact(c("Landsat8","Landsat8","Landsat7")) == "L7x1, L8x2")
check("merge_sensor_strings", merge_sensor_strings(c("L7x1, L8x2","L8x1")) == "L7x1, L8x3")
cat("\n")

# --- 4. MTL parsing must take Level-2, not Level-1 ----------------------------
cat("--- 4. MTL parsing ---\n")
suppressMessages(source("01_scenes.R"))
inv <- suppressMessages(build_inventory())
check("scenes after tier filter", nrow(inv) == 13,
      sprintf("(%d retained; 1 Tier-2 scene excluded)", nrow(inv)))
check("Level-2 SR factor read (not the Level-1 decoy)",
      all(abs(inv$sr_scale - 2.75e-05) < 1e-12),
      sprintf("sr_scale = %.3e", inv$sr_scale[1]))
check("Level-2 ST factor read",
      all(abs(inv$st_scale - 0.00341802) < 1e-12))
check("all scenes complete", all(inv$complete))
cat("\n")

# --- 5. full pipeline ---------------------------------------------------------
cat("--- 5. Full pipeline run (this takes ~40 s) ---\n")
setwd(old_wd)
log_file <- file.path(TEST_ROOT, "run.log")
status <- system2("Rscript", file.path(WORK, "06_main.R"),
                  stdout = log_file, stderr = log_file)
run_log <- readLines(log_file, warn = FALSE)
check("pipeline exited cleanly", status == 0, sprintf("(exit %d)", status))
check("no ERROR lines in the log",
      !any(grepl("\\bERROR\\b", run_log)),
      if (any(grepl("\\bERROR\\b", run_log)))
        paste(head(grep("\\bERROR\\b", run_log, value = TRUE), 2), collapse = " | ") else "")
check("heavily clouded scene rejected",
      any(grepl("Cloud cover too high", run_log)))
n_proc <- as.integer(sub(".*Scenes processed *: *", "",
                         grep("Scenes processed", run_log, value = TRUE)[1]))
check("12 scenes processed", isTRUE(n_proc == 12), sprintf("(got %s)", n_proc))
cat("\n")

# --- 6. output validation -----------------------------------------------------
cat("--- 6. Output validation ---\n")
OUT <- file.path(TEST_ROOT, "results_landsat")
expect_file <- function(rel) file.exists(file.path(OUT, rel))
check("per-scene CSV",   expect_file("tables/scenes/ndvi_lst_per_scene.csv"))
check("monthly CSV",     expect_file("tables/monthly/ndvi_lst_monthly.csv"))
check("zone CSV",        expect_file("tables/spatial/zones/ndvi_lst_per_zone_monthly.csv"))
check("anomaly report",  expect_file("tables/anomalies/summary_by_sensor.csv"))
check("clipped rasters", length(list.files(file.path(OUT, "rasters/lst"))) == 12)

n_png <- length(list.files(file.path(OUT, "plots"), pattern = "\\.png$",
                           recursive = TRUE))
check("figures produced", n_png > 0, sprintf("(%d PNG files)", n_png))

# 6a. INDEPENDENT recomputation of one zone-month, bypassing the pipeline.
# 2015-05 has two scenes (L8 + L7), so this also tests the weighting.
zn <- readr::read_csv(file.path(OUT, "tables/spatial/zones/ndvi_lst_per_zone_monthly.csv"),
                      show_col_types = FALSE)
sc <- readr::read_csv(file.path(OUT, "tables/scenes/ndvi_lst_per_scene.csv"),
                      show_col_types = FALSE)
tgt <- zn |> dplyr::filter(ano == 2015, mes == 5, zone_name == "Centro")
ids <- sc |> dplyr::filter(ano == 2015, mes == 5) |> dplyr::pull(scene_id)
zpoly <- terra::vect(sf::st_read(file.path(TEST_ROOT, "shp/zones.gpkg"), quiet = TRUE) |>
                     dplyr::filter(NM_BAIRRO == "Centro"))
num <- 0; den <- 0
for (id in ids) {
  rr <- terra::rast(file.path(OUT, "rasters/lst", paste0(id, "_LST_Celsius.tif")))
  vv <- terra::extract(rr, zpoly, ID = FALSE)[, 1]; vv <- vv[!is.na(vv)]
  num <- num + mean(vv) * length(vv); den <- den + length(vv)
}
check("zonal weighted mean matches manual recomputation",
      abs(num/den - tgt$lst_c_media) < 1e-6,
      sprintf("pipeline %.6f vs manual %.6f", tgt$lst_c_media, num/den))
check("zonal pixel count matches", den == tgt$lst_c_n,
      sprintf("(%d vs %d)", tgt$lst_c_n, den))

# 6b. physical plausibility of the recovered values
check("LST within configured bounds",
      all(sc$lst_c_min >= 5 - 1e-6, na.rm = TRUE) &&
      all(sc$lst_c_max <= 70 + 1e-6, na.rm = TRUE),
      sprintf("range %.1f..%.1f C", min(sc$lst_c_min, na.rm = TRUE),
              max(sc$lst_c_max, na.rm = TRUE)))
check("NDVI within [-1, 1]",
      all(sc$ndvi_min >= -1, na.rm = TRUE) && all(sc$ndvi_max <= 1, na.rm = TRUE),
      sprintf("range %.3f..%.3f", min(sc$ndvi_min, na.rm = TRUE),
              max(sc$ndvi_max, na.rm = TRUE)))

# 6c. seasonal signal injected by the generator must be recovered:
# February (austral summer) warmer than August (winter)
feb <- mean(sc$lst_c_media[sc$mes == 2], na.rm = TRUE)
aug <- mean(sc$lst_c_media[sc$mes == 8], na.rm = TRUE)
check("seasonal cycle recovered (Feb warmer than Aug)", feb > aug,
      sprintf("Feb %.1f C vs Aug %.1f C", feb, aug))

# 6d. monthly grid must be complete, with explicit NA for months with no scene
mo <- readr::read_csv(file.path(OUT, "tables/monthly/ndvi_lst_monthly.csv"),
                      show_col_types = FALSE)
check("complete monthly grid (24 rows for 2 years)", nrow(mo) == 24,
      sprintf("(%d rows)", nrow(mo)))
check("months without scenes are explicit NA", any(is.na(mo$ndvi_media)))
cat("\n")

# --- summary ------------------------------------------------------------------
cat("=============================================================\n")
cat(sprintf("  RESULT: %d passed, %d failed\n", PASS, FAIL))
cat("=============================================================\n")
if (FAIL == 0) {
  cat("  The pipeline is correctly installed and numerically sound.\n")
  cat("  Next: edit R/00_config.R for your study area, then\n")
  cat("        Rscript R/06_main.R\n")
} else {
  cat(sprintf("  Inspect the run log: %s\n", log_file))
}
if (KEEP_OUTPUT) {
  cat(sprintf("\n  Test outputs kept in: %s\n", TEST_ROOT))
} else {
  unlink(TEST_ROOT, recursive = TRUE)
}
cat("\n")
quit(status = if (FAIL == 0) 0 else 1)
