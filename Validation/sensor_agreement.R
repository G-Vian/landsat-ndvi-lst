# =============================================================================
# sensor_agreement.R - QUANTIFYING INTER-SENSOR DIFFERENCES
# =============================================================================
# PURPOSE
#   Some months are composited from scenes acquired by more than one Landsat
#   instrument. This script measures how much LST and NDVI differ between those
#   instruments, at the neighbourhood level, and whether that difference
#   survives into the covariate actually used in the model.
#
# -----------------------------------------------------------------------------
# THE CENTRAL CAVEAT - READ BEFORE INTERPRETING ANY OUTPUT
# -----------------------------------------------------------------------------
#   Two scenes from different satellites are NEVER simultaneous. Landsat 7 and
#   Landsat 8 are offset by about 8 days over the same path. A difference
#   measured between them therefore confounds TWO things:
#
#       (a) instrument bias   - different spectral response, calibration
#       (b) real change       - the surface genuinely changed in those days
#
#   A raw difference CANNOT separate them. Three quantities are reported, of
#   which only the third is interpretable on its own:
#
#   1. RAW DIFFERENCE (confounded). An UPPER BOUND on the sensor effect.
#
#   2. SAME-SENSOR REFERENCE (control). Where a month holds two scenes from the
#      SAME instrument, their difference is attributable to date alone. This
#      calibrates how much of (1) is simply the passage of time.
#
#   3. ANOMALY DIFFERENCE (decision-relevant). The model uses each
#      neighbourhood's deviation from the city-wide mean of the same scene. If
#      an instrument shifts a whole scene, that shift enters both the
#      neighbourhood value and the city mean and cancels. Comparing anomalies
#      TESTS whether that cancellation happens, instead of assuming it.
#
#   Interpretation:
#     raw large + anomaly small  -> common-mode shift; the covariate is robust
#     raw large + anomaly large  -> instruments disagree about SPATIAL PATTERN;
#                                   this would justify harmonisation
#     raw small                  -> no issue to begin with
#
# -----------------------------------------------------------------------------
# ADDITIONAL VALIDATION (sections 4B and 5B)
# -----------------------------------------------------------------------------
#   The mean absolute difference alone does not say whether a discrepancy
#   MATTERS. A discrepancy of 0.03 NDVI units is negligible if neighbourhoods
#   differ from one another by 0.30, and fatal if they differ by 0.05. Five
#   further diagnostics are therefore computed:
#
#   (i)   SIGNAL-TO-NOISE and RELIABILITY. The spatial spread of the anomaly
#         across neighbourhoods is the signal the model exploits; the paired
#         difference estimates the measurement error. Their ratio, and the
#         classical reliability coefficient lambda = var_true / var_observed,
#         quantify how much a regression coefficient fitted on this covariate
#         is attenuated towards zero by measurement error. Attenuation is
#         roughly (1 - lambda) expressed as a percentage.
#
#   (ii)  RANK AND CONCORDANCE AGREEMENT, computed WITHIN each scene pair and
#         then summarised across pairs. For a spatial model what matters is
#         whether the instruments agree on the ORDERING of neighbourhoods, not
#         on the absolute level. Spearman rho answers exactly that. Lin's
#         concordance correlation coefficient additionally penalises departures
#         from the identity line, so it detects offset and scale problems that
#         Pearson r would miss.
#
#   (iii) SCALE (GAIN) DIFFERENCE via the standardised major axis slope. A
#         slope of one means the instruments differ at most by an offset, which
#         the anomaly removes. A slope away from one means they differ in
#         GAIN - one compresses or stretches the spatial contrast - and no
#         additive correction can fix that.
#
#   (iv)  DIRECTIONAL BIAS with a canonical ordering. In section 3 the labels
#         a and b follow the order scenes happen to appear, so the sign of
#         mean_diff there is NOT "newer instrument minus older". Section 4B
#         recomputes the signed difference with the older instrument always
#         first, which is the only form that can be quoted as an instrument
#         bias. Significance is assessed at the SCENE-PAIR level, because
#         neighbourhoods within a scene pair are not independent observations.
#
#   (v)   LIMITS OF AGREEMENT (Bland-Altman) and the DEPENDENCE ON THE
#         ACQUISITION GAP. If a difference is driven by real surface change it
#         should grow with the number of days separating the two acquisitions;
#         if it is instrumental it should not.
#
#   None of this replaces section 4 - it is added alongside it.
#
# -----------------------------------------------------------------------------
# COMPATIBILITY
#   Works with BOTH pipeline versions: the original Portuguese one
#   (04_clip_santos.R, tabelas/, load_santos_bairros) and the generalised one
#   (04_clip_aoi.R, tables/, load_aoi_zones). Folder names, function names, the
#   per-scene table and the neighbourhood-name column are DETECTED, not assumed.
#
# USAGE
#   Rscript sensor_agreement.R
#   Rscript sensor_agreement.R /path/to/pipeline/modules
#
#   Inside the container:
#     singularity exec --bind /home/g.vian \
#       --env R_LIBS="/home/g.vian/R_libs:/usr/local/lib/R/site-library:/usr/local/lib/R/library" \
#       /home/public/R_inla/r_inla.sif \
#       Rscript sensor_agreement.R
#
#   Requires a COMPLETED main run with clipped rasters saved to disk.
# =============================================================================


# =============================================================================
# 0. LOCATE THE PIPELINE AND ADAPT TO ITS VERSION
# =============================================================================

# ### EDIT ### absolute path to the folder holding 00_config.R. Used only if
# the script cannot work it out from the command line or its own location.
PIPELINE_DIR_FALLBACK <-
  "/home/g.vian/Pesquisa_Epidemic/PROJETO_SANTOS/Dados_Satélite"

.args <- commandArgs(trailingOnly = TRUE)

.script_dir <- tryCatch({
  a <- commandArgs(trailingOnly = FALSE)
  f <- a[grep("--file=", a)]
  if (length(f) > 0)
    suppressWarnings(dirname(normalizePath(sub("--file=", "", f[1]),
                                           mustWork = FALSE)))
  else NA_character_
}, error = function(e) NA_character_)

PIPELINE_DIR <- NULL
for (cand in c(if (length(.args) > 0) .args[1] else NULL,
               .script_dir, getwd(), PIPELINE_DIR_FALLBACK)) {
  if (!is.null(cand) && !is.na(cand) && nzchar(cand) &&
      file.exists(file.path(cand, "00_config.R"))) {
    PIPELINE_DIR <- normalizePath(cand); break
  }
}
if (is.null(PIPELINE_DIR))
  stop("Could not locate 00_config.R.\n",
       "  Pass the module directory as an argument:\n",
       "    Rscript sensor_agreement.R /path/to/modules\n",
       "  or fix PIPELINE_DIR_FALLBACK at the top of this script.")

message("\n=== INTER-SENSOR AGREEMENT ANALYSIS ===\n")
message("Pipeline modules : ", PIPELINE_DIR)

# The modules use relative source() calls internally, so the working directory
# must be the module folder.
setwd(PIPELINE_DIR)
source(file.path(PIPELINE_DIR, "00_config.R"))

# 00_config.R performs its own setwd() based on the script it thinks is
# running, which is this file rather than the pipeline entry point. Force the
# working directory back to the module folder, otherwise the relative
# source("00_config.R") inside 04_clip_*.R cannot resolve.
setwd(PIPELINE_DIR)

load_packages()
suppressPackageStartupMessages({
  library(dplyr); library(readr); library(terra); library(sf); library(tibble)
})

# --- clipping module: file name differs between pipeline versions ------------
.clip <- c("04_clip_aoi.R", "04_clip_santos.R")
.clip <- .clip[file.exists(file.path(PIPELINE_DIR, .clip))]
if (length(.clip) == 0)
  stop("No clipping module (04_*.R) found in ", PIPELINE_DIR)
source(file.path(PIPELINE_DIR, .clip[1]))
message("Clipping module  : ", .clip[1])

# --- zone loader: function name differs between versions ---------------------
.load_zones <- {
  if      (exists("load_aoi_zones"))      load_aoi_zones
  else if (exists("load_santos_bairros")) load_santos_bairros
  else stop("No zone-loading function found ",
            "(expected load_aoi_zones or load_santos_bairros).")
}

# --- output directory: variable name differs ---------------------------------
OUT_DIR <- {
  if      (exists("OUTPUT_DIR")) OUTPUT_DIR
  else if (exists("DIR_OUTPUT")) DIR_OUTPUT
  else stop("Neither OUTPUT_DIR nor DIR_OUTPUT is defined in 00_config.R.")
}
message("Output directory : ", OUT_DIR)

AOI_LABEL <- {
  if      (exists("AOI_NAME"))  AOI_NAME
  else if (exists("NOME_AREA")) NOME_AREA
  else "study area"
}

# --- tables folder: "tables" (English) or "tabelas" (Portuguese) -------------
.tbl_root <- c(file.path(OUT_DIR, "tables"), file.path(OUT_DIR, "tabelas"))
.tbl_root <- .tbl_root[dir.exists(.tbl_root)]
if (length(.tbl_root) == 0)
  stop("No tables/ or tabelas/ folder under ", OUT_DIR,
       "\nHas the main pipeline been run?")
TBL_DIR <- .tbl_root[1]

# --- read a CSV with whichever separator it actually uses --------------------
# Some tables were written with write.csv (comma) and others with write.csv2
# (semicolon, comma decimal), so the separator is sniffed from the header.
.read_any <- function(p, n = Inf) {
  h <- readLines(p, n = 1, warn = FALSE)
  if (grepl(";", h)) readr::read_csv2(p, n_max = n, show_col_types = FALSE)
  else               readr::read_csv (p, n_max = n, show_col_types = FALSE)
}
# --- per-scene table: located by CONTENT, so the file name does not matter ---
# The required columns are exactly the ones this script consumes. Note that
# other per-scene tables exist (e.g. the anomaly audit) which also carry
# scene_id/sensor/ano/mes but NOT the statistics, so the test must be specific.
.need <- c("sensor", "ano", "mes", "data_aq", "ndvi_media", "lst_c_media")
.id_candidates <- c("scene_id", "cena_id", "id_cena", "ID_cena")

.cands <- list.files(TBL_DIR, pattern = "\\.csv$",
                     recursive = TRUE, full.names = TRUE)
# Prefer files whose path suggests a per-scene table, but do not require it.
.cands <- c(.cands[grepl("cena|scene", .cands, ignore.case = TRUE)],
            .cands[!grepl("cena|scene", .cands, ignore.case = TRUE)])

F_SCENES <- NULL; COL_ID <- NULL
for (p in .cands) {
  hdr <- tryCatch(names(.read_any(p, n = 1)), error = function(e) character(0))
  idc <- intersect(.id_candidates, hdr)
  if (length(idc) >= 1 && all(.need %in% hdr)) {
    F_SCENES <- p; COL_ID <- idc[1]; break
  }
}
if (is.null(F_SCENES))
  stop("Could not find a per-scene statistics table under ", TBL_DIR,
       "\n  Required columns: a scene identifier plus ",
       paste(.need, collapse = ", "),
       "\n  (the anomaly audit table has some of these but not the statistics)")
message("Per-scene table  : ", basename(F_SCENES), "  (id column: ", COL_ID, ")")

# --- raster folders and file-name suffixes: read off what is on disk --------
DIR_NDVI <- file.path(OUT_DIR, "rasters", "ndvi")
DIR_LST  <- file.path(OUT_DIR, "rasters", "lst")
if (!dir.exists(DIR_NDVI) || !dir.exists(DIR_LST))
  stop("Raster folders not found under ", file.path(OUT_DIR, "rasters"),
       "\nThe main run must have saved clipped rasters.")

.f1 <- list.files(DIR_NDVI, pattern = "\\.tif$")[1]
.f2 <- list.files(DIR_LST,  pattern = "\\.tif$")[1]
if (is.na(.f1) || is.na(.f2))
  stop("Raster folders exist but contain no .tif files.")

.sc_head <- .read_any(F_SCENES)
names(.sc_head)[names(.sc_head) == COL_ID] <- "scene_id"
.ids <- unique(.sc_head$scene_id)

# The scene id prefixes the file name; the suffix is whatever follows it.
.strip_suffix <- function(fname, ids) {
  hit <- ids[vapply(ids, function(i) startsWith(fname, i), logical(1))]
  if (length(hit) == 0) return(NA_character_)
  sub(paste0("^", hit[which.max(nchar(hit))]), "", fname)
}
SUF_NDVI <- .strip_suffix(.f1, .ids)
SUF_LST  <- .strip_suffix(.f2, .ids)
if (is.na(SUF_NDVI) || is.na(SUF_LST)) {
  SUF_NDVI <- "_NDVI.tif"; SUF_LST <- "_LST_Celsius.tif"
  message("Raster suffixes  : not derivable; using conventional defaults")
} else {
  message("Raster suffixes  : ", SUF_NDVI, " | ", SUF_LST)
}

# --- output folders ----------------------------------------------------------
DIR_OUT <- file.path(TBL_DIR, "sensor_agreement")
DIR_FIG <- file.path(OUT_DIR, "plots", "sensor_agreement")
dir.create(DIR_OUT, recursive = TRUE, showWarnings = FALSE)
dir.create(DIR_FIG, recursive = TRUE, showWarnings = FALSE)
message("")


# -----------------------------------------------------------------------------
# 1. IDENTIFY MONTHS WITH MORE THAN ONE SCENE
# -----------------------------------------------------------------------------

sc <- .sc_head    # already read, id column normalised to `scene_id`

# A scene is only usable here if it carries a value for at least one variable
sc <- sc |> dplyr::filter(!is.na(ndvi_media) | !is.na(lst_c_media))

multi <- sc |>
  dplyr::group_by(ano, mes) |>
  dplyr::filter(dplyr::n() >= 2) |>
  dplyr::ungroup()

if (nrow(multi) == 0) {
  message("No month contains more than one scene. Nothing to compare.")
  quit(status = 0)
}

inventory <- multi |>
  dplyr::group_by(ano, mes) |>
  dplyr::summarise(
    n_scenes  = dplyr::n(),
    n_sensors = dplyr::n_distinct(sensor),
    sensors   = paste(sort(unique(sensor)), collapse = " + "),
    dates     = paste(sort(as.character(data_aq)), collapse = ", "),
    span_days = as.integer(max(as.Date(data_aq)) - min(as.Date(data_aq))),
    .groups   = "drop"
  ) |>
  dplyr::mutate(month_type = ifelse(n_sensors >= 2,
                                    "Mixed instruments", "Single instrument"))

readr::write_csv(inventory, file.path(DIR_OUT, "sensor_pairs_inventory.csv"))

message(sprintf("Months with >= 2 scenes : %d", nrow(inventory)))
message(sprintf("  mixed instruments     : %d",
                sum(inventory$month_type == "Mixed instruments")))
message(sprintf("  single instrument     : %d   (control group)",
                sum(inventory$month_type == "Single instrument")))
message(sprintf("Total months in series  : %d",
                dplyr::n_distinct(paste(sc$ano, sc$mes))))
message("")

# -----------------------------------------------------------------------------
# 2. PER-ZONE VALUES FOR EVERY SCENE IN THOSE MONTHS
# -----------------------------------------------------------------------------
# Both the absolute zone mean and the zone's anomaly relative to the city-wide
# mean OF THAT SAME SCENE are computed. The anomaly is the quantity the model
# is built on (see the covariate construction in the Methods), so it is what
# ultimately matters.

zones <- .load_zones()

# The neighbourhood-name column differs between pipeline versions; detect it.
.zcols <- names(sf::st_drop_geometry(zones))
.name_candidates <- c("zone_name", "nome_bairro", "NOME_BAIRRO", "Bairro_Limpo",
                      "BAIRRO", "bairro", "NOME", "Nome", "nome", "NM_BAIRRO")
COL_ZONE <- intersect(.name_candidates, .zcols)
if (length(COL_ZONE) == 0) {
  chr <- .zcols[vapply(sf::st_drop_geometry(zones), is.character, logical(1))]
  if (length(chr) == 0)
    stop("Could not identify a neighbourhood-name column. Columns present: ",
         paste(.zcols, collapse = ", "))
  COL_ZONE <- chr[1]
  message("Zone-name column : '", COL_ZONE, "'  (guessed - please verify)")
} else {
  COL_ZONE <- COL_ZONE[1]
  message("Zone-name column : '", COL_ZONE, "'")
}
znames <- as.character(sf::st_drop_geometry(zones)[[COL_ZONE]])
zvect  <- terra::vect(zones)
message(sprintf("Zones loaded     : %d\n", length(znames)))

zone_values <- function(scene_id, kind) {
  path <- if (kind == "lst") file.path(DIR_LST,  paste0(scene_id, SUF_LST))
  else               file.path(DIR_NDVI, paste0(scene_id, SUF_NDVI))
  if (!file.exists(path)) return(NULL)
  
  r  <- terra::rast(path)
  ex <- terra::extract(r, zvect, fun = NULL, na.rm = TRUE, ID = TRUE)
  names(ex)[2] <- "value"
  ex <- ex[!is.na(ex$value), , drop = FALSE]
  if (nrow(ex) == 0) return(NULL)
  
  agg <- ex |>
    dplyr::group_by(ID) |>
    dplyr::summarise(value = mean(value), n_px = dplyr::n(), .groups = "drop")
  
  # City-wide mean of THIS scene, weighted by valid pixel count, then the
  # zone's deviation from it: the same construction used for the covariate.
  city <- sum(agg$value * agg$n_px) / sum(agg$n_px)
  
  tibble::tibble(
    scene_id  = scene_id,
    variable  = kind,
    zone_name = znames[agg$ID],
    value     = agg$value,
    anomaly   = agg$value - city,
    n_px      = agg$n_px
  )
}

message("Extracting neighbourhood values scene by scene...")
vals <- list(); n_missing <- 0L
for (i in seq_len(nrow(multi))) {
  s <- multi[i, ]
  for (k in c("lst", "ndvi")) {
    v <- zone_values(s$scene_id, k)
    if (is.null(v)) { n_missing <- n_missing + 1L; next }
    vals[[length(vals) + 1]] <- v |>
      dplyr::mutate(ano = s$ano, mes = s$mes,
                    sensor = s$sensor, data_aq = s$data_aq)
  }
  if (i %% 25 == 0) message(sprintf("   %d / %d scenes", i, nrow(multi)))
}
if (length(vals) == 0)
  stop("No clipped rasters could be read. Check ",
       file.path(OUT_DIR, "rasters"))
if (n_missing > 0)
  message(sprintf("   (%d scene-variable rasters absent, skipped)", n_missing))
zv <- dplyr::bind_rows(vals)
message("")

# -----------------------------------------------------------------------------
# 3. PAIR SCENES WITHIN EACH MONTH
# -----------------------------------------------------------------------------
# All unordered pairs of scenes in the same month, matched zone by zone. Pairs
# are labelled by whether they cross instruments or not - the same-instrument
# pairs are the control that isolates the effect of acquisition date.

make_pairs <- function(df) {
  out <- list()
  keys <- df |> dplyr::distinct(ano, mes, variable)
  for (r in seq_len(nrow(keys))) {
    k   <- keys[r, ]
    sub <- df |> dplyr::filter(ano == k$ano, mes == k$mes, variable == k$variable)
    ids <- unique(sub$scene_id)
    if (length(ids) < 2) next
    cmb <- utils::combn(ids, 2)
    for (c_ in seq_len(ncol(cmb))) {
      a <- sub |> dplyr::filter(scene_id == cmb[1, c_])
      b <- sub |> dplyr::filter(scene_id == cmb[2, c_])
      j <- dplyr::inner_join(a, b, by = "zone_name", suffix = c("_a", "_b"))
      if (nrow(j) == 0) next
      out[[length(out) + 1]] <- j |>
        dplyr::transmute(
          ano = k$ano, mes = k$mes, variable = k$variable, zone_name,
          scene_a = scene_id_a, scene_b = scene_id_b,
          sensor_a = sensor_a,  sensor_b = sensor_b,
          date_a = data_aq_a,   date_b = data_aq_b,
          gap_days = as.integer(abs(as.Date(data_aq_b) - as.Date(data_aq_a))),
          value_a, value_b, anomaly_a, anomaly_b,
          n_px = pmin(n_px_a, n_px_b)
        )
    }
  }
  dplyr::bind_rows(out)
}

pairs <- make_pairs(zv)
if (nrow(pairs) == 0) stop("No comparable scene pairs could be formed.")

short <- function(x) gsub("Landsat", "L", x)

pairs <- pairs |>
  dplyr::mutate(
    comparison = ifelse(sensor_a == sensor_b,
                        "Same instrument", "Between instruments"),
    pair_label = ifelse(
      sensor_a == sensor_b,
      paste0(short(sensor_a), " vs ", short(sensor_b)),
      paste0(pmin(short(sensor_a), short(sensor_b)), " vs ",
             pmax(short(sensor_a), short(sensor_b)))),
    
    # --- absolute difference ---
    diff_abs = value_b - value_a,
    
    # --- percentage difference ---
    # For LST this is relative to the pair mean in degrees Celsius. For NDVI a
    # percentage is unstable near zero, so it is computed only where the pair
    # mean exceeds 0.1 in absolute value; elsewhere it is left missing and the
    # absolute difference should be used instead.
    mean_ab  = (value_a + value_b) / 2,
    diff_pct = ifelse(abs(mean_ab) > ifelse(variable == "ndvi", 0.1, 1e-6),
                      100 * (value_b - value_a) / abs(mean_ab), NA_real_),
    
    # --- anomaly difference: the decision-relevant quantity ---
    diff_anom = anomaly_b - anomaly_a
  )

readr::write_csv(pairs, file.path(DIR_OUT, "sensor_differences_per_zone.csv"))
message(sprintf("Neighbourhood-level paired comparisons: %d", nrow(pairs)))
message(sprintf("  between instruments : %d",
                sum(pairs$comparison == "Between instruments")))
message(sprintf("  same instrument     : %d\n",
                sum(pairs$comparison == "Same instrument")))

# -----------------------------------------------------------------------------
# 4. SUMMARY
# -----------------------------------------------------------------------------

summ <- pairs |>
  dplyr::group_by(variable, comparison) |>
  dplyr::summarise(
    n_pairs      = dplyr::n(),
    n_months     = dplyr::n_distinct(paste(ano, mes)),
    gap_days_med = stats::median(gap_days),
    mean_diff    = mean(diff_abs, na.rm = TRUE),          # signed: bias
    mad_diff     = mean(abs(diff_abs), na.rm = TRUE),     # unsigned: magnitude
    sd_diff      = stats::sd(diff_abs, na.rm = TRUE),
    p95_absdiff  = stats::quantile(abs(diff_abs), 0.95, na.rm = TRUE),
    mad_pct      = mean(abs(diff_pct), na.rm = TRUE),
    mad_anom     = mean(abs(diff_anom), na.rm = TRUE),    # after city-mean removal
    reduction    = 100 * (1 - mean(abs(diff_anom), na.rm = TRUE) /
                            mean(abs(diff_abs),  na.rm = TRUE)),
    .groups = "drop"
  ) |>
  dplyr::arrange(variable, comparison)

readr::write_csv(summ, file.path(DIR_OUT, "sensor_agreement_summary.csv"))

# Breakdown by instrument pair, to spot a specific problematic transition
by_pair <- pairs |>
  dplyr::group_by(variable, pair_label, comparison) |>
  dplyr::summarise(
    n_pairs   = dplyr::n(),
    mean_diff = mean(diff_abs, na.rm = TRUE),
    mad_diff  = mean(abs(diff_abs), na.rm = TRUE),
    mad_anom  = mean(abs(diff_anom), na.rm = TRUE),
    reduction = 100 * (1 - mean(abs(diff_anom), na.rm = TRUE) /
                         mean(abs(diff_abs),  na.rm = TRUE)),
    .groups = "drop"
  ) |>
  dplyr::arrange(variable, dplyr::desc(n_pairs))
readr::write_csv(by_pair, file.path(DIR_OUT, "sensor_agreement_by_pair.csv"))

message("--- SUMMARY ---")
print(as.data.frame(summ |>
                      dplyr::mutate(dplyr::across(dplyr::where(is.numeric), ~round(., 4)))),
      row.names = FALSE)
message("\n--- BY INSTRUMENT PAIR ---")
print(as.data.frame(by_pair |>
                      dplyr::mutate(dplyr::across(dplyr::where(is.numeric), ~round(., 4)))),
      row.names = FALSE)
message("")


# =============================================================================
# 4B. ADDITIONAL VALIDATION
# =============================================================================
# Everything below is added on top of section 4 and changes none of its output.

# --- small helpers, base R only so no new dependency -------------------------

# Lin's concordance correlation coefficient. Unlike Pearson r it is not
# invariant to a shift or a rescaling, so it drops when the points sit on a
# line that is not the identity line - which is precisely the failure mode of
# interest here.
.ccc <- function(a, b) {
  ok <- is.finite(a) & is.finite(b)
  a <- a[ok]; b <- b[ok]
  if (length(a) < 3) return(NA_real_)
  va <- stats::var(a); vb <- stats::var(b)
  den <- va + vb + (mean(a) - mean(b))^2
  if (!is.finite(den) || den <= 0) return(NA_real_)
  2 * stats::cov(a, b) / den
}

# Standardised major axis (reduced major axis) slope. Symmetric in a and b and
# not attenuated by error in the predictor, unlike an ordinary least squares
# slope, which makes it the right tool when BOTH variables are measured with
# error - the situation here.
.sma_slope <- function(a, b) {
  ok <- is.finite(a) & is.finite(b)
  a <- a[ok]; b <- b[ok]
  if (length(a) < 3) return(NA_real_)
  sa <- stats::sd(a); sb <- stats::sd(b)
  if (!is.finite(sa) || sa <= 0 || !is.finite(sb)) return(NA_real_)
  cv <- stats::cov(a, b)
  if (!is.finite(cv) || cv == 0) return(NA_real_)
  sign(cv) * sb / sa
}

.safe_cor <- function(a, b, method) {
  ok <- is.finite(a) & is.finite(b)
  if (sum(ok) < 3) return(NA_real_)
  out <- suppressWarnings(stats::cor(a[ok], b[ok], method = method))
  if (is.finite(out)) out else NA_real_
}

# Slope of y on x, fitted in an ordinary function rather than inside a
# summarise block, so that lm() resolves its formula in a normal environment
# instead of a dplyr data mask.
.slope_of <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) <= 5) return(NA_real_)
  xx <- x[ok]; yy <- y[ok]
  if (stats::var(xx) <= 0) return(NA_real_)
  unname(stats::coef(stats::lm(yy ~ xx))[2])
}

# Landsat generation as an integer, used only to fix a canonical direction for
# the signed bias. Anything unparseable sorts last.
.sensor_rank <- function(x) {
  n <- suppressWarnings(as.integer(gsub("[^0-9]", "", as.character(x))))
  ifelse(is.na(n), 99L, n)
}

# A scene pair must contain at least this many neighbourhoods before a
# correlation computed across them is worth reporting.
MIN_ZONES_PER_PAIR <- 10

# --- (a) per scene-pair agreement --------------------------------------------
# Correlations are computed WITHIN a scene pair, across neighbourhoods, and
# only then summarised across pairs. Pooling every zone-pair from every month
# into one correlation would mix spatial variation with seasonal variation and
# inflate the result.

scene_pair_metrics <- pairs |>
  dplyr::group_by(variable, comparison, pair_label,
                  ano, mes, scene_a, scene_b) |>
  dplyr::summarise(
    n_zones   = dplyr::n(),
    gap_days  = dplyr::first(gap_days),
    r_abs     = .safe_cor(value_a,   value_b,   "pearson"),
    rho_abs   = .safe_cor(value_a,   value_b,   "spearman"),
    r_anom    = .safe_cor(anomaly_a, anomaly_b, "pearson"),
    rho_anom  = .safe_cor(anomaly_a, anomaly_b, "spearman"),
    ccc_abs   = .ccc(value_a,   value_b),
    ccc_anom  = .ccc(anomaly_a, anomaly_b),
    sma_abs   = .sma_slope(value_a,   value_b),
    sma_anom  = .sma_slope(anomaly_a, anomaly_b),
    mdiff_abs  = mean(diff_abs,  na.rm = TRUE),
    mdiff_anom = mean(diff_anom, na.rm = TRUE),
    mad_abs_p  = mean(abs(diff_abs),  na.rm = TRUE),
    .groups   = "drop"
  ) |>
  dplyr::filter(n_zones >= MIN_ZONES_PER_PAIR)

readr::write_csv(scene_pair_metrics,
                 file.path(DIR_OUT, "sensor_agreement_per_scene_pair.csv"))

concord <- scene_pair_metrics |>
  dplyr::group_by(variable, comparison) |>
  dplyr::summarise(
    n_scene_pairs = dplyr::n(),
    zones_med     = stats::median(n_zones),
    r_abs_med     = stats::median(r_abs,    na.rm = TRUE),
    rho_abs_med   = stats::median(rho_abs,  na.rm = TRUE),
    ccc_abs_med   = stats::median(ccc_abs,  na.rm = TRUE),
    sma_abs_med   = stats::median(sma_abs,  na.rm = TRUE),
    r_anom_med    = stats::median(r_anom,   na.rm = TRUE),
    rho_anom_med  = stats::median(rho_anom, na.rm = TRUE),
    rho_anom_q25  = unname(stats::quantile(rho_anom, 0.25, na.rm = TRUE)),
    rho_anom_q75  = unname(stats::quantile(rho_anom, 0.75, na.rm = TRUE)),
    ccc_anom_med  = stats::median(ccc_anom, na.rm = TRUE),
    sma_anom_med  = stats::median(sma_anom, na.rm = TRUE),
    .groups = "drop"
  ) |>
  dplyr::arrange(variable, comparison)

readr::write_csv(concord, file.path(DIR_OUT, "sensor_concordance.csv"))

# --- (b) signal-to-noise and reliability -------------------------------------
# The signal is the spread of the anomaly ACROSS neighbourhoods, taken over all
# scenes: that is the contrast the spatial model is fitted to. The measurement
# error is recovered from the paired difference, whose variance is twice the
# error variance when the two errors are independent, hence the sqrt(2).
# Reliability is then the classical lambda = var_true / var_observed, and a
# regression coefficient fitted on a covariate with reliability lambda is
# attenuated towards zero by approximately a factor lambda.

signal_ref <- zv |>
  dplyr::group_by(variable) |>
  dplyr::summarise(
    sd_signal_abs  = stats::sd(value,   na.rm = TRUE),
    sd_signal_anom = stats::sd(anomaly, na.rm = TRUE),
    .groups = "drop"
  )

reliab <- pairs |>
  dplyr::group_by(variable, comparison) |>
  dplyr::summarise(
    n_pairs      = dplyr::n(),
    sd_err_abs   = stats::sd(diff_abs,  na.rm = TRUE) / sqrt(2),
    sd_err_anom  = stats::sd(diff_anom, na.rm = TRUE) / sqrt(2),
    .groups = "drop"
  ) |>
  dplyr::left_join(signal_ref, by = "variable") |>
  dplyr::mutate(
    var_true_anom  = pmax(sd_signal_anom^2 - sd_err_anom^2, 0),
    sd_true_anom   = sqrt(var_true_anom),
    snr_anom       = ifelse(sd_err_anom > 0, sd_true_anom / sd_err_anom, NA_real_),
    lambda_anom    = ifelse(sd_signal_anom > 0,
                            var_true_anom / sd_signal_anom^2, NA_real_),
    attenuation_pc = 100 * (1 - lambda_anom),
    var_true_abs   = pmax(sd_signal_abs^2 - sd_err_abs^2, 0),
    snr_abs        = ifelse(sd_err_abs > 0,
                            sqrt(var_true_abs) / sd_err_abs, NA_real_),
    lambda_abs     = ifelse(sd_signal_abs > 0,
                            var_true_abs / sd_signal_abs^2, NA_real_)
  ) |>
  dplyr::arrange(variable, comparison)

readr::write_csv(reliab, file.path(DIR_OUT, "sensor_reliability.csv"))

# --- (c) directional bias with a canonical ordering --------------------------
# In section 3 the roles of a and b follow the order in which scenes happen to
# be enumerated, so the sign of mean_diff there is arbitrary. Here the older
# instrument is always placed first, so the reported difference reads as
# "newer minus older" and can legitimately be quoted as an instrument bias.
# The test is applied to SCENE-PAIR means rather than to individual
# neighbourhoods, because neighbourhoods inside one scene pair share the same
# atmosphere, the same acquisition dates and the same calibration, and are
# therefore very far from independent; treating them as independent would
# shrink the standard error by roughly the square root of the number of zones.

.older_first <- pairs |>
  dplyr::filter(comparison == "Between instruments") |>
  dplyr::mutate(
    a_is_older = .sensor_rank(sensor_a) <= .sensor_rank(sensor_b),
    s_old   = ifelse(a_is_older, sensor_a, sensor_b),
    s_new   = ifelse(a_is_older, sensor_b, sensor_a),
    d_abs   = ifelse(a_is_older, value_b   - value_a,   value_a   - value_b),
    d_anom  = ifelse(a_is_older, anomaly_b - anomaly_a, anomaly_a - anomaly_b),
    direction = paste0(short(s_new), " - ", short(s_old))
  )

if (nrow(.older_first) > 0) {
  
  dir_scene <- .older_first |>
    dplyr::group_by(variable, direction, ano, mes, scene_a, scene_b) |>
    dplyr::summarise(d_abs = mean(d_abs, na.rm = TRUE),
                     d_anom = mean(d_anom, na.rm = TRUE),
                     .groups = "drop")
  
  directional <- dir_scene |>
    dplyr::group_by(variable, direction) |>
    dplyr::summarise(
      n_scene_pairs = dplyr::n(),
      bias_abs      = mean(d_abs, na.rm = TRUE),
      se_abs        = stats::sd(d_abs, na.rm = TRUE) / sqrt(dplyr::n()),
      bias_anom     = mean(d_anom, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      tstat = ifelse(n_scene_pairs > 1 & se_abs > 0, bias_abs / se_abs, NA_real_),
      pval  = ifelse(is.na(tstat), NA_real_,
                     2 * stats::pt(-abs(tstat), df = n_scene_pairs - 1)),
      tcrit = ifelse(n_scene_pairs > 1,
                     stats::qt(0.975, df = n_scene_pairs - 1), NA_real_),
      ci_lo = bias_abs - tcrit * se_abs,
      ci_hi = bias_abs + tcrit * se_abs
    ) |>
    dplyr::arrange(variable, dplyr::desc(n_scene_pairs))
  
  readr::write_csv(directional, file.path(DIR_OUT, "sensor_directional_bias.csv"))
} else {
  directional <- NULL
}

# --- (d) limits of agreement and dependence on the acquisition gap -----------
# Bland-Altman limits give the interval containing about 95 per cent of the
# discrepancies, which is a more honest statement of worst-case disagreement
# than a mean. The slope against the gap tests the competing explanation: a
# difference caused by real surface change should widen as the two acquisitions
# move apart in time, whereas an instrumental difference should not.

loa <- summ |>
  dplyr::transmute(
    variable, comparison, n_pairs,
    bias      = mean_diff,
    loa_lower = mean_diff - 1.96 * sd_diff,
    loa_upper = mean_diff + 1.96 * sd_diff,
    p95_absdiff
  )

gap_slope <- pairs |>
  dplyr::group_by(variable, comparison) |>
  dplyr::summarise(
    n_pairs      = dplyr::n(),
    gap_min      = min(gap_days, na.rm = TRUE),
    gap_max      = max(gap_days, na.rm = TRUE),
    gap_med      = stats::median(gap_days, na.rm = TRUE),
    slope_per_day = .slope_of(gap_days, abs(diff_abs)),
    cor_gap = .safe_cor(gap_days, abs(diff_abs), "spearman"),
    .groups = "drop"
  ) |>
  dplyr::arrange(variable, comparison)

validation <- loa |>
  dplyr::left_join(gap_slope |> dplyr::select(-n_pairs),
                   by = c("variable", "comparison")) |>
  dplyr::left_join(reliab |> dplyr::select(variable, comparison,
                                           sd_signal_anom, sd_err_anom,
                                           snr_anom, lambda_anom,
                                           attenuation_pc),
                   by = c("variable", "comparison")) |>
  dplyr::left_join(concord |> dplyr::select(variable, comparison,
                                            n_scene_pairs, rho_anom_med,
                                            ccc_anom_med, sma_anom_med,
                                            rho_abs_med, ccc_abs_med),
                   by = c("variable", "comparison"))

readr::write_csv(validation, file.path(DIR_OUT, "sensor_validation_summary.csv"))

message("--- RELIABILITY OF THE MODELLED COVARIATE (anomaly scale) ---")
print(as.data.frame(reliab |>
                      dplyr::select(variable, comparison, sd_signal_anom,
                                    sd_err_anom, snr_anom, lambda_anom,
                                    attenuation_pc) |>
                      dplyr::mutate(dplyr::across(dplyr::where(is.numeric), ~round(., 4)))),
      row.names = FALSE)

message("\n--- SPATIAL CONCORDANCE WITHIN SCENE PAIRS (medians) ---")
print(as.data.frame(concord |>
                      dplyr::select(variable, comparison, n_scene_pairs,
                                    rho_abs_med, ccc_abs_med, sma_abs_med,
                                    rho_anom_med, ccc_anom_med, sma_anom_med) |>
                      dplyr::mutate(dplyr::across(dplyr::where(is.numeric), ~round(., 3)))),
      row.names = FALSE)

if (!is.null(directional)) {
  message("\n--- DIRECTIONAL BIAS, NEWER MINUS OLDER INSTRUMENT ---")
  print(as.data.frame(directional |>
                        dplyr::select(variable, direction, n_scene_pairs,
                                      bias_abs, ci_lo, ci_hi, pval, bias_anom) |>
                        dplyr::mutate(dplyr::across(dplyr::where(is.numeric), ~round(., 4)))),
        row.names = FALSE)
}

message("\n--- LIMITS OF AGREEMENT AND GAP DEPENDENCE ---")
print(as.data.frame(validation |>
                      dplyr::select(variable, comparison, bias, loa_lower,
                                    loa_upper, gap_med, slope_per_day, cor_gap) |>
                      dplyr::mutate(dplyr::across(dplyr::where(is.numeric), ~round(., 4)))),
      row.names = FALSE)
message("")


# -----------------------------------------------------------------------------
# 5. LATEX TABLE
# -----------------------------------------------------------------------------

fmt <- function(x, d = 2) ifelse(is.na(x), "---", formatC(x, format = "f", digits = d))
vlab <- function(v) ifelse(v == "lst", "LST (\\si{\\celsius})", "NDVI")

tex <- c(
  "\\begin{table}[htbp]",
  "\\centering",
  "\\footnotesize",
  "\\caption{Agreement between Landsat scenes acquired within the same calendar",
  "month, measured at the zone level. \\emph{Between instruments} compares scenes",
  "from different satellites; \\emph{same instrument} compares scenes from the same",
  "satellite and therefore isolates the effect of acquisition date, serving as a",
  "reference for how much of the between-instrument difference is attributable to",
  "real surface change rather than to the sensor. The final columns report the",
  "same comparison after each zone is expressed as a deviation from the city-wide",
  "mean of its own scene, which is the form in which the variable enters the",
  "model; the reduction column gives the percentage decrease in mean absolute",
  "difference achieved by that transformation.}",
  "\\label{tab:sensor-agreement}",
  "\\begin{tabular}{llrrrrrr}",
  "\\toprule",
  " & & & \\multicolumn{3}{c}{Absolute values} & \\multicolumn{2}{c}{Anomalies} \\\\",
  "\\cmidrule(lr){4-6}\\cmidrule(lr){7-8}",
  "Variable & Comparison & $n$ & Mean & Mean abs. & Mean abs. & Mean abs. & Reduction \\\\",
  " & & pairs & diff. & diff. & diff. (\\%) & diff. & (\\%) \\\\",
  "\\midrule")

for (v in unique(summ$variable)) {
  s <- summ |> dplyr::filter(variable == v)
  for (i in seq_len(nrow(s))) {
    r <- s[i, ]
    d <- if (v == "ndvi") 4 else 2
    tex <- c(tex, sprintf("%s & %s & %d & %s & %s & %s & %s & %s \\\\",
                          if (i == 1) vlab(v) else "",
                          r$comparison, r$n_pairs,
                          fmt(r$mean_diff, d), fmt(r$mad_diff, d),
                          fmt(r$mad_pct, 1), fmt(r$mad_anom, d),
                          fmt(r$reduction, 1)))
  }
  if (v != utils::tail(unique(summ$variable), 1)) tex <- c(tex, "\\midrule")
}

tex <- c(tex, "\\bottomrule", "\\end{tabular}", "\\end{table}")
writeLines(tex, file.path(DIR_OUT, "TABLE_sensor_agreement.tex"))
message("LaTeX table written: TABLE_sensor_agreement.tex\n")


# =============================================================================
# 5B. LATEX TABLES FOR THE ADDITIONAL VALIDATION
# =============================================================================
# Three further tables. They are written individually and also concatenated
# into TABLE_sensor_validation_ALL.tex for convenience.

vars_present <- unique(summ$variable)

# --- Table 2: reliability and attenuation ------------------------------------
tex2 <- c(
  "\\begin{table}[htbp]",
  "\\centering",
  "\\footnotesize",
  "\\caption{Reliability of the covariate as it enters the model. The signal is",
  "the standard deviation of the zone anomaly across neighbourhoods and scenes;",
  "the noise is the standard deviation of the paired difference divided by",
  "$\\sqrt{2}$, which estimates the error of a single measurement. Reliability is",
  "$\\lambda = \\sigma^2_{\\text{true}} / \\sigma^2_{\\text{observed}}$; a regression",
  "coefficient fitted on a covariate with reliability $\\lambda$ is attenuated",
  "towards the null by approximately the factor $\\lambda$, and the final column",
  "reports that attenuation as a percentage. The same-instrument rows give the",
  "attenuation that would remain even if every scene came from one satellite,",
  "and therefore isolate the part of the degradation that harmonisation could",
  "not remove.}",
  "\\label{tab:sensor-reliability}",
  "\\begin{tabular}{llrrrrr}",
  "\\toprule",
  "Variable & Comparison & Signal & Noise & SNR & $\\lambda$ & Attenuation \\\\",
  " & & SD & SD & & & (\\%) \\\\",
  "\\midrule")

for (v in vars_present) {
  s <- reliab |> dplyr::filter(variable == v) |> dplyr::arrange(comparison)
  for (i in seq_len(nrow(s))) {
    r <- s[i, ]
    d <- if (v == "ndvi") 4 else 2
    tex2 <- c(tex2, sprintf("%s & %s & %s & %s & %s & %s & %s \\\\",
                            if (i == 1) vlab(v) else "",
                            r$comparison,
                            fmt(r$sd_signal_anom, d), fmt(r$sd_err_anom, d),
                            fmt(r$snr_anom, 2), fmt(r$lambda_anom, 3),
                            fmt(r$attenuation_pc, 1)))
  }
  if (v != utils::tail(vars_present, 1)) tex2 <- c(tex2, "\\midrule")
}
tex2 <- c(tex2, "\\bottomrule", "\\end{tabular}", "\\end{table}")
writeLines(tex2, file.path(DIR_OUT, "TABLE_sensor_reliability.tex"))

# --- Table 3: spatial concordance --------------------------------------------
tex3 <- c(
  "\\begin{table}[htbp]",
  "\\centering",
  "\\footnotesize",
  "\\caption{Spatial concordance between paired scenes, computed across",
  "neighbourhoods within each scene pair and then summarised as the median over",
  "pairs. Spearman $\\rho$ asks whether the two acquisitions rank the",
  "neighbourhoods in the same order, which is what a spatial model depends on.",
  "Lin's concordance correlation coefficient additionally penalises departures",
  "from the identity line, so it falls when the instruments agree on the pattern",
  "but not on its offset or scale. The standardised major axis slope isolates",
  "the scale component: a value of one means the instruments differ at most by",
  "an additive offset, which the anomaly transformation removes, whereas a value",
  "away from one means they differ in gain, which no additive correction can",
  "repair. Statistics are reported for the absolute values and for the anomalies.}",
  "\\label{tab:sensor-concordance}",
  "\\begin{tabular}{llrrrrrrr}",
  "\\toprule",
  " & & & \\multicolumn{3}{c}{Absolute values} & \\multicolumn{3}{c}{Anomalies} \\\\",
  "\\cmidrule(lr){4-6}\\cmidrule(lr){7-9}",
  "Variable & Comparison & Scene & $\\rho$ & CCC & SMA & $\\rho$ & CCC & SMA \\\\",
  " & & pairs & & & slope & & & slope \\\\",
  "\\midrule")

for (v in vars_present) {
  s <- concord |> dplyr::filter(variable == v) |> dplyr::arrange(comparison)
  for (i in seq_len(nrow(s))) {
    r <- s[i, ]
    tex3 <- c(tex3, sprintf("%s & %s & %d & %s & %s & %s & %s & %s & %s \\\\",
                            if (i == 1) vlab(v) else "",
                            r$comparison, as.integer(r$n_scene_pairs),
                            fmt(r$rho_abs_med, 3), fmt(r$ccc_abs_med, 3),
                            fmt(r$sma_abs_med, 3),
                            fmt(r$rho_anom_med, 3), fmt(r$ccc_anom_med, 3),
                            fmt(r$sma_anom_med, 3)))
  }
  if (v != utils::tail(vars_present, 1)) tex3 <- c(tex3, "\\midrule")
}
tex3 <- c(tex3, "\\bottomrule", "\\end{tabular}", "\\end{table}")
writeLines(tex3, file.path(DIR_OUT, "TABLE_sensor_concordance.tex"))

# --- Table 4: directional bias, limits of agreement, gap dependence ----------
tex4 <- c(
  "\\begin{table}[htbp]",
  "\\centering",
  "\\footnotesize",
  "\\caption{Directional bias and worst-case disagreement. The upper block",
  "reports the signed difference with the older instrument always taken as the",
  "reference, so the quantity reads as newer minus older; confidence intervals",
  "and $p$-values are computed over scene pairs rather than over neighbourhoods,",
  "since neighbourhoods within a scene pair share the same atmosphere and",
  "calibration and are not independent. The lower block gives Bland--Altman",
  "limits of agreement, which bracket about 95 per cent of the observed",
  "discrepancies, together with the slope of the absolute difference on the",
  "number of days separating the two acquisitions: a difference produced by real",
  "surface change should widen as that gap grows, whereas an instrumental",
  "difference should not.}",
  "\\label{tab:sensor-bias}",
  "\\begin{tabular}{llrrrr}",
  "\\toprule",
  "\\multicolumn{6}{l}{\\emph{Directional bias between instruments (newer $-$ older)}} \\\\",
  "\\midrule",
  "Variable & Direction & Scene & Bias & 95\\% CI & $p$ \\\\",
  " & & pairs & & & \\\\",
  "\\midrule")

if (!is.null(directional) && nrow(directional) > 0) {
  for (v in vars_present) {
    s <- directional |> dplyr::filter(variable == v)
    if (nrow(s) == 0) next
    for (i in seq_len(nrow(s))) {
      r <- s[i, ]
      d <- if (v == "ndvi") 4 else 2
      pv <- if (is.na(r$pval)) "---" else
        if (r$pval < 0.001) "$<$0.001" else formatC(r$pval, format = "f", digits = 3)
      tex4 <- c(tex4, sprintf("%s & %s & %d & %s & [%s, %s] & %s \\\\",
                              if (i == 1) vlab(v) else "",
                              r$direction, as.integer(r$n_scene_pairs),
                              fmt(r$bias_abs, d),
                              fmt(r$ci_lo, d), fmt(r$ci_hi, d), pv))
    }
  }
} else {
  tex4 <- c(tex4, "\\multicolumn{6}{l}{No between-instrument pairs available.} \\\\")
}

tex4 <- c(tex4,
          "\\midrule",
          "\\multicolumn{6}{l}{\\emph{Limits of agreement and dependence on the acquisition gap}} \\\\",
          "\\midrule",
          "Variable & Comparison & Median & \\multicolumn{2}{c}{95\\% limits of agreement} & Slope \\\\",
          " & & gap (d) & lower & upper & per day \\\\",
          "\\midrule")

for (v in vars_present) {
  s <- validation |> dplyr::filter(variable == v) |> dplyr::arrange(comparison)
  for (i in seq_len(nrow(s))) {
    r <- s[i, ]
    d <- if (v == "ndvi") 4 else 2
    tex4 <- c(tex4, sprintf("%s & %s & %s & %s & %s & %s \\\\",
                            if (i == 1) vlab(v) else "",
                            r$comparison, fmt(r$gap_med, 0),
                            fmt(r$loa_lower, d), fmt(r$loa_upper, d),
                            fmt(r$slope_per_day, if (v == "ndvi") 5 else 3)))
  }
  if (v != utils::tail(vars_present, 1)) tex4 <- c(tex4, "\\midrule")
}
tex4 <- c(tex4, "\\bottomrule", "\\end{tabular}", "\\end{table}")
writeLines(tex4, file.path(DIR_OUT, "TABLE_sensor_bias_loa.tex"))

writeLines(c(tex, "", tex2, "", tex3, "", tex4),
           file.path(DIR_OUT, "TABLE_sensor_validation_ALL.tex"))

message("LaTeX tables written:")
message("  TABLE_sensor_reliability.tex")
message("  TABLE_sensor_concordance.tex")
message("  TABLE_sensor_bias_loa.tex")
message("  TABLE_sensor_validation_ALL.tex   <- all four tables in one file\n")


# -----------------------------------------------------------------------------
# 6. FIGURE
# -----------------------------------------------------------------------------

if (requireNamespace("ggplot2", quietly = TRUE)) {
  suppressPackageStartupMessages(library(ggplot2))
  
  pl <- pairs |>
    dplyr::mutate(var_lab = ifelse(variable == "lst",
                                   "LST (\u00b0C)", "NDVI (dimensionless)"))
  
  COL <- c("Between instruments" = "#b2182b", "Same instrument" = "#2166ac")
  th  <- theme_bw(base_size = 10) +
    theme(plot.title    = element_text(face = "bold", size = 11),
          plot.subtitle = element_text(colour = "grey40", size = 8.5),
          plot.caption  = element_text(colour = "grey55", size = 7, hjust = 0),
          strip.text    = element_text(face = "bold"),
          panel.grid.minor = element_blank(),
          legend.position  = "bottom",
          legend.title     = element_blank(),
          plot.background  = element_rect(fill = "white", colour = NA))
  
  # (A) distribution of absolute-value differences
  pA <- ggplot(pl, aes(x = diff_abs, fill = comparison, colour = comparison)) +
    geom_density(alpha = 0.30, linewidth = 0.5) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
    facet_wrap(~var_lab, scales = "free") +
    scale_fill_manual(values = COL) + scale_colour_manual(values = COL) +
    labs(title = "A. Difference in absolute values between scenes of the same month",
         subtitle = "Between-instrument differences confound sensor bias with real change over the acquisition gap",
         x = "Difference between paired scenes", y = "Density") + th
  
  # (B) same, after removing the city-wide mean of each scene
  pB <- ggplot(pl, aes(x = diff_anom, fill = comparison, colour = comparison)) +
    geom_density(alpha = 0.30, linewidth = 0.5) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
    facet_wrap(~var_lab, scales = "free") +
    scale_fill_manual(values = COL) + scale_colour_manual(values = COL) +
    labs(title = "B. Difference in anomalies (deviation from the city-wide mean of each scene)",
         subtitle = "This is the form in which the variable enters the model; narrowing relative to panel A indicates a common-mode shift that cancels",
         x = "Difference in anomaly between paired scenes", y = "Density") + th
  
  # (C) scene-to-scene agreement, zone by zone
  pC <- ggplot(pl, aes(x = value_a, y = value_b, colour = comparison)) +
    geom_abline(slope = 1, intercept = 0, colour = "grey50", linetype = "dashed") +
    geom_point(alpha = 0.55, size = 1.1) +
    facet_wrap(~var_lab, scales = "free") +
    scale_colour_manual(values = COL) +
    labs(title = "C. Zone-level agreement between paired scenes",
         subtitle = "Each point is one zone in one month; the dashed line is exact agreement",
         x = "First scene of the pair", y = "Second scene of the pair") + th
  
  # (D) magnitude vs acquisition gap: does the difference grow with elapsed time?
  pD <- ggplot(pl, aes(x = gap_days, y = abs(diff_abs), colour = comparison)) +
    geom_point(alpha = 0.45, size = 1.1) +
    facet_wrap(~var_lab, scales = "free_y") +
    scale_colour_manual(values = COL) +
    labs(title = "D. Difference magnitude against the gap between acquisitions",
         subtitle = "A difference driven by real surface change should grow with the number of days separating the scenes",
         x = "Days between the two acquisitions", y = "Absolute difference") + th
  
  cap <- paste(
    "Landsat Collection 2 Level-2 (USGS). Between-instrument differences cannot be",
    "attributed to the sensor alone: paired scenes are acquired days apart, so the",
    "difference also contains genuine surface change. Same-instrument pairs, where",
    "available, provide the reference for that component.")
  
  save1 <- function(p, f, w, h) ggsave(file.path(DIR_FIG, f), p, width = w,
                                       height = h, dpi = 300, bg = "white")
  save1(pA, "FIG_sensor_diff_absolute.png", 9, 4)
  save1(pB, "FIG_sensor_diff_anomaly.png",  9, 4)
  save1(pC, "FIG_sensor_scatter.png",       9, 4)
  save1(pD, "FIG_sensor_gap.png",           9, 4)
  
  if (requireNamespace("patchwork", quietly = TRUE)) {
    library(patchwork)
    fig <- (pA / pB / pC / pD) +
      plot_annotation(
        title = sprintf("Inter-sensor agreement within calendar months \u2014 %s", AOI_LABEL),
        caption = cap,
        theme = theme(plot.title   = element_text(face = "bold", size = 13),
                      plot.caption = element_text(size = 7, colour = "grey45",
                                                  hjust = 0)))
    ggsave(file.path(DIR_FIG, "FIG_sensor_agreement.png"), fig,
           width = 9.5, height = 15, dpi = 300, bg = "white")
  }
  
  # --- (E) and (F): figures for the additional validation --------------------
  # E shows the per-scene-pair rank agreement, which is the quantity a spatial
  # model actually relies on. F shows the standardised major axis slope, whose
  # distance from one measures a gain difference that the anomaly cannot fix.
  
  spm <- scene_pair_metrics |>
    dplyr::mutate(var_lab = ifelse(variable == "lst",
                                   "LST (\u00b0C)", "NDVI (dimensionless)"))
  
  if (nrow(spm) > 0) {
    pE <- ggplot(spm, aes(x = comparison, y = rho_anom, fill = comparison)) +
      geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
      geom_boxplot(alpha = 0.35, outlier.size = 0.7, width = 0.55) +
      facet_wrap(~var_lab) +
      scale_fill_manual(values = COL) +
      labs(title = "E. Rank agreement between paired scenes, anomaly scale",
           subtitle = "One point per scene pair: Spearman correlation of neighbourhood anomalies. A spatial model relies on this ordering being stable",
           x = NULL, y = expression(Spearman~rho)) + th
    
    pF <- ggplot(spm, aes(x = comparison, y = sma_anom, fill = comparison)) +
      geom_hline(yintercept = 1, linetype = "dashed", colour = "grey40") +
      geom_boxplot(alpha = 0.35, outlier.size = 0.7, width = 0.55) +
      facet_wrap(~var_lab) +
      scale_fill_manual(values = COL) +
      labs(title = "F. Scale agreement: standardised major axis slope, anomaly scale",
           subtitle = "A slope of one (dashed) means the instruments differ only by an offset, which the anomaly removes; a departure indicates a gain difference, which it does not",
           x = NULL, y = "SMA slope") + th
    
    save1(pE, "FIG_sensor_rank_agreement.png", 9, 4)
    save1(pF, "FIG_sensor_scale_agreement.png", 9, 4)
    
    if (requireNamespace("patchwork", quietly = TRUE)) {
      fig2 <- (pE / pF) +
        plot_annotation(
          title = sprintf("Additional validation \u2014 %s", AOI_LABEL),
          caption = cap,
          theme = theme(plot.title   = element_text(face = "bold", size = 13),
                        plot.caption = element_text(size = 7, colour = "grey45",
                                                    hjust = 0)))
      ggsave(file.path(DIR_FIG, "FIG_sensor_validation.png"), fig2,
             width = 9.5, height = 8, dpi = 300, bg = "white")
    }
  }
  
  message("Figures written to: ", DIR_FIG, "\n")
} else {
  message("ggplot2 unavailable - tables written, figures skipped.")
}

# -----------------------------------------------------------------------------
# 7. PLAIN-LANGUAGE READOUT
# -----------------------------------------------------------------------------

message("\n=== HOW TO READ THIS ===")
for (v in unique(summ$variable)) {
  s   <- summ |> dplyr::filter(variable == v)
  bet <- s |> dplyr::filter(comparison == "Between instruments")
  sam <- s |> dplyr::filter(comparison == "Same instrument")
  vn  <- ifelse(v == "lst", "LST", "NDVI")
  
  if (nrow(bet) == 1) {
    message(sprintf("\n%s, between instruments:", vn))
    message(sprintf("  mean absolute difference %.4f (%.1f%%), over %d zone-pairs",
                    bet$mad_diff, bet$mad_pct, bet$n_pairs))
    message(sprintf("  after removing each scene's city-wide mean: %.4f  (%.1f%% smaller)",
                    bet$mad_anom, bet$reduction))
    if (is.finite(bet$reduction) && bet$reduction > 50)
      message("  -> most of the discrepancy is a shift common to the whole scene,")
    else
      message("  -> the discrepancy is NOT merely a common shift; the instruments")
    if (is.finite(bet$reduction) && bet$reduction > 50)
      message("     which cancels in the anomaly used by the model.")
    else
      message("     differ in spatial pattern, which does not cancel.")
  }
  if (nrow(sam) == 1 && nrow(bet) == 1) {
    message(sprintf("  same-instrument reference: %.4f", sam$mad_diff))
    if (sam$mad_diff >= bet$mad_diff * 0.8)
      message("  -> comparable to the between-instrument value: the difference is")
    else
      message("  -> clearly smaller: part of the between-instrument difference is")
    if (sam$mad_diff >= bet$mad_diff * 0.8)
      message("     largely explained by acquisition date, not by the sensor.")
    else
      message("     plausibly attributable to the instrument.")
  } else if (nrow(bet) == 1) {
    message("  no same-instrument pairs available: the date component cannot be")
    message("  separated, so the reported value is an UPPER BOUND on sensor effect.")
  }
}
message("")

# --- readout for the additional validation -----------------------------------

message("\n=== WHAT THE ADDITIONAL TESTS SAY ===")
for (v in unique(summ$variable)) {
  vn <- ifelse(v == "lst", "LST", "NDVI")
  rb <- reliab  |> dplyr::filter(variable == v,
                                 comparison == "Between instruments")
  rs <- reliab  |> dplyr::filter(variable == v,
                                 comparison == "Same instrument")
  cb <- concord |> dplyr::filter(variable == v,
                                 comparison == "Between instruments")
  
  if (nrow(rb) != 1) next
  message(sprintf("\n%s:", vn))
  message(sprintf("  spatial signal SD %.4f against measurement noise SD %.4f",
                  rb$sd_signal_anom, rb$sd_err_anom))
  message(sprintf("  signal-to-noise %.2f, reliability %.3f",
                  rb$snr_anom, rb$lambda_anom))
  
  if (is.finite(rb$attenuation_pc)) {
    message(sprintf("  -> a coefficient fitted on this covariate is attenuated"))
    message(sprintf("     towards the null by roughly %.0f%%.", rb$attenuation_pc))
    if (rb$attenuation_pc < 10)
      message("     That is negligible; no correction is called for.")
    else if (rb$attenuation_pc < 25)
      message("     That is modest but worth stating in the limitations.")
    else
      message("     That is substantial; effect sizes should be read as lower bounds.")
  }
  
  if (nrow(rs) == 1 && is.finite(rs$attenuation_pc))
    message(sprintf("  same-instrument reference attenuation: %.0f%% (the floor that",
                    rs$attenuation_pc))
  if (nrow(rs) == 1 && is.finite(rs$attenuation_pc))
    message("     would persist even with a single satellite)")
  
  if (nrow(cb) == 1 && is.finite(cb$rho_anom_med)) {
    message(sprintf("  median rank agreement across scene pairs: rho = %.2f",
                    cb$rho_anom_med))
    if (cb$rho_anom_med > 0.8)
      message("  -> the instruments order the neighbourhoods almost identically.")
    else if (cb$rho_anom_med > 0.5)
      message("  -> the ordering is broadly preserved but not tightly.")
    else
      message("  -> the ordering is NOT preserved; this is the serious case.")
  }
  
  if (nrow(cb) == 1 && is.finite(cb$sma_anom_med)) {
    message(sprintf("  median SMA slope: %.2f", cb$sma_anom_med))
    if (abs(cb$sma_anom_med - 1) < 0.15)
      message("  -> close to one: the disagreement is an offset, which the anomaly removes.")
    else
      message("  -> away from one: a gain difference the anomaly cannot remove.")
  }
}

message("")
message("Suggested next step: refit the model on months with a single instrument")
message("only, and compare the posteriors against the full series. If they agree,")
message("the inter-sensor issue does not drive the conclusions.")
message("")
message("\nOutputs: ", DIR_OUT, "\n")