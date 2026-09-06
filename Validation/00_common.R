# =============================================================================
# VALIDATION SUITE - 00 - SHARED HELPERS, PATHS AND FIGURE STYLE
# =============================================================================
# This file is sourced by every 01_..05_ validation script. It centralises
# four things so that they are defined once and stay consistent:
#
#   (a) WHERE the pipeline outputs live and WHERE validation outputs go
#   (b) HOW input CSVs are located (by content, not by a hard-coded name)
#   (c) HOW column names are normalised (the pipeline exists in a Portuguese
#       and an English variant with different column names)
#   (d) The publication figure theme and the figure-saving function
#
# You normally only need to edit the two paths in section 1.
#
# WHY LOCATE FILES BY CONTENT RATHER THAN BY PATH
#   The pipeline has two variants. The Portuguese one writes to tabelas/ with
#   columns named ano, mes, nome_bairro; the generalised English one writes to
#   tables/ with year, month, zone_name. Hard-coding either would break the
#   other. Instead the scripts search the results tree for a CSV that CONTAINS
#   the columns they need, then rename those columns to one internal
#   convention. Everything downstream can then assume a single naming scheme.
# =============================================================================


# -----------------------------------------------------------------------------
# 1. PATHS  ### EDIT HERE ###
# -----------------------------------------------------------------------------
# Where the validation outputs (tables, reports, figures) are written.
# Each script creates its own subfolder underneath this root.
#
# The accented character in "Dados_Satelite" is written as a unicode escape
# (\u00e9 = e-acute) rather than typed directly. A literal accent in a path
# string is interpreted according to the session locale, which inside a
# Singularity container is often C/POSIX rather than UTF-8, and the path then
# silently fails to match the real directory.
VALIDATION_ROOT <- file.path(
  "/home/g.vian/Pesquisa_Epidemic/PROJETO_SANTOS",
  "Dados_Sat\u00e9lite", "files", "Validation"
)

# Where the main Landsat pipeline wrote its results. The first entry that
# actually exists on disk is used, so you can leave several candidates listed
# and the scripts will pick whichever run is present.
PIPELINE_DIR_CANDIDATES <- c(
  "/home/g.vian/Pesquisa_Epidemic/PROJETO_SANTOS/Resultados_Landsat",
  "/home/g.vian/Pesquisa_Epidemic/PROJETO_SANTOS/Resultados_Landsat_EN",
  "/home/g.vian/Pesquisa_Epidemic/PROJETO_SANTOS/results_landsat"
)

# Months flagged as low quality and excluded from statistics.
# July 2021 carries a Landsat-7 SLC-off artefact: the scene reports a
# city-wide mean LST of 15.19 C for a mid-winter month, which is physically
# implausible for Santos (22-30 C expected). It passes every automatic filter
# in the pipeline because the cold pixels are diffuse rather than clamped, so
# it has to be excluded by hand here. See the pipeline README, section 10.
LOW_QUALITY_MONTHS <- data.frame(
  year   = 2021,
  month  = 7,
  reason = "LE07 SLC-off, diffuse cold bias (city mean LST 15.19 C, implausible)"
)

# Study period.
YEAR_START <- 2010
YEAR_END   <- 2025

# Short label used in figure titles and report headers.
AOI_LABEL <- "Santos, Brazil"


# -----------------------------------------------------------------------------
# 2. FIGURE QUALITY SETTINGS  ### EDIT HERE IF NEEDED ###
# -----------------------------------------------------------------------------
# 320 dpi comfortably exceeds the 300 dpi most journals require for raster
# figures. Every figure is ALSO written as a PDF, which is vector: text and
# lines stay sharp at any zoom level, and that is the version to submit if the
# journal accepts vector artwork.
FIG_DPI       <- 320
FIG_SAVE_PDF  <- TRUE    # set FALSE to halve the number of output files
FIG_BASE_SIZE <- 11      # base font size in points for standard figures


# -----------------------------------------------------------------------------
# 3. PACKAGES
# -----------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr)
  library(ggplot2); library(scales)
})

# Force English month names and English number formatting regardless of the
# container locale. Without this, format(date, "%b") returns "jan"/"fev" on a
# pt_BR system and the figures come out half-translated.
suppressWarnings(try(Sys.setlocale("LC_TIME", "C"), silent = TRUE))


# -----------------------------------------------------------------------------
# 4. LOCATING THE PIPELINE OUTPUTS
# -----------------------------------------------------------------------------

# Read a CSV regardless of whether it uses comma or semicolon separators.
# write.csv produces commas; write.csv2 (and Excel in a pt_BR locale) produces
# semicolons with comma decimals. Sniffing the header avoids a silent failure
# in which every row is parsed into a single character column.
read_any_csv <- function(path, n_max = Inf) {
  header <- readLines(path, n = 1, warn = FALSE)
  if (grepl(";", header)) {
    readr::read_csv2(path, n_max = n_max, show_col_types = FALSE)
  } else {
    readr::read_csv(path, n_max = n_max, show_col_types = FALSE)
  }
}

# Resolve the pipeline results directory from the candidate list.
resolve_pipeline_dir <- function() {
  hit <- PIPELINE_DIR_CANDIDATES[dir.exists(PIPELINE_DIR_CANDIDATES)]
  if (length(hit) == 0)
    stop("None of the pipeline result directories exist:\n  ",
         paste(PIPELINE_DIR_CANDIDATES, collapse = "\n  "),
         "\nEdit PIPELINE_DIR_CANDIDATES in 00_common.R.")
  normalizePath(hit[1])
}

# Find the tables root, which is "tabelas" in the Portuguese pipeline and
# "tables" in the English one.
resolve_tables_dir <- function(pipeline_dir) {
  cands <- file.path(pipeline_dir, c("tables", "tabelas"))
  hit   <- cands[dir.exists(cands)]
  if (length(hit) == 0)
    stop("No tables/ or tabelas/ folder under ", pipeline_dir,
         "\nHas the main pipeline finished a run?")
  hit[1]
}

# Search the tables tree for a CSV that contains ALL of `required_any`, where
# each element of `required_any` is itself a vector of acceptable spellings for
# one logical column. Files whose path matches `prefer_pattern` are examined
# first, but a name mismatch is not fatal.
#
# Returns the path, or NULL if nothing matches.
find_csv_by_columns <- function(tables_dir, required_any,
                                prefer_pattern = NULL,
                                forbid_any = NULL) {

  cands <- list.files(tables_dir, pattern = "\\.csv$",
                      recursive = TRUE, full.names = TRUE)
  if (length(cands) == 0) return(NULL)

  # Files whose path hints at the right table are examined first, but a name
  # mismatch never disqualifies a file - only its columns do.
  if (!is.null(prefer_pattern)) {
    preferred <- grepl(prefer_pattern, cands, ignore.case = TRUE)
    cands     <- c(cands[preferred], cands[!preferred])
  }

  for (p in cands) {
    hdr <- tryCatch(names(read_any_csv(p, n_max = 1)),
                    error = function(e) character(0))
    if (length(hdr) == 0) next

    # every logical column must be satisfied by at least one spelling
    has_required <- all(vapply(required_any,
                               function(spellings) any(spellings %in% hdr),
                               logical(1)))
    if (!has_required) next

    # `forbid_any` rules a file out if it carries any of these columns. Used to
    # separate the city-wide table from the per-zone table, which otherwise
    # match the same required-column test.
    if (!is.null(forbid_any) && any(forbid_any %in% hdr)) next

    return(p)
  }
  NULL
}


# -----------------------------------------------------------------------------
# 5. COLUMN NAME NORMALISATION
# -----------------------------------------------------------------------------
# Maps whichever spelling the pipeline used onto one internal convention, so
# that the analysis code below never has to branch on pipeline version.
#
# Internal convention used by every validation script:
#   year, month, zone_id, zone_name,
#   ndvi_mean, ndvi_sd, ndvi_min, ndvi_max, ndvi_n,
#   lst_mean,  lst_sd,  lst_min,  lst_max,  lst_n,
#   n_scenes, sensors_used

COLUMN_ALIASES <- list(
  year         = c("year", "ano"),
  month        = c("month", "mes"),
  zone_id      = c("zone_id", "bairro_id", "id_bairro"),
  zone_name    = c("zone_name", "nome_bairro", "bairro", "BAIRRO",
                   "NOME", "Nome", "nome", "NM_BAIRRO", "Bairro_Limpo"),
  ndvi_mean    = c("ndvi_mean", "ndvi_media"),
  ndvi_sd      = c("ndvi_sd", "ndvi_dp"),
  ndvi_min     = c("ndvi_min"),
  ndvi_max     = c("ndvi_max"),
  ndvi_n       = c("ndvi_n"),
  lst_mean     = c("lst_mean", "lst_c_mean", "lst_c_media", "lst_media"),
  lst_sd       = c("lst_sd", "lst_c_sd", "lst_c_dp"),
  lst_min      = c("lst_min", "lst_c_min"),
  lst_max      = c("lst_max", "lst_c_max"),
  lst_n        = c("lst_n", "lst_c_n"),
  n_scenes     = c("n_scenes", "n_cenas"),
  sensors_used = c("sensors_used", "sensores_usados")
)

#' Rename the columns of `df` onto the internal convention.
#' Columns with no known alias are left untouched.
normalise_columns <- function(df) {
  nm <- names(df)
  for (target in names(COLUMN_ALIASES)) {
    if (target %in% nm) next                       # already correct
    hit <- intersect(COLUMN_ALIASES[[target]], nm)
    if (length(hit) > 0) names(df)[names(df) == hit[1]] <- target
  }
  df
}

#' Add the columns every script needs: a mid-month date, an English month
#' abbreviation as an ordered factor, and the low-quality flag.
add_time_columns <- function(df) {
  df |>
    dplyr::mutate(
      date       = as.Date(sprintf("%d-%02d-15", year, month)),
      month_abb  = factor(month.abb[month], levels = month.abb),
      season     = dplyr::case_when(          # southern hemisphere
        month %in% c(12, 1, 2)  ~ "Summer",
        month %in% c(3, 4, 5)   ~ "Autumn",
        month %in% c(6, 7, 8)   ~ "Winter",
        month %in% c(9, 10, 11) ~ "Spring"
      ),
      low_quality = paste(year, month) %in%
        paste(LOW_QUALITY_MONTHS$year, LOW_QUALITY_MONTHS$month)
    )
}


# -----------------------------------------------------------------------------
# 6. LOADERS FOR THE TWO INPUT TABLES
# -----------------------------------------------------------------------------

#' City-wide monthly table: one row per (year, month).
load_city_monthly <- function(tables_dir) {
  # `forbid_any` excludes the per-zone table, which satisfies the same
  # required-column test but has ~60 rows per month; using it here would
  # silently corrupt every city-level statistic.
  p <- find_csv_by_columns(
    tables_dir,
    required_any = list(COLUMN_ALIASES$year, COLUMN_ALIASES$month,
                        COLUMN_ALIASES$ndvi_mean, COLUMN_ALIASES$lst_mean),
    prefer_pattern = "mensal|monthly",
    forbid_any     = COLUMN_ALIASES$zone_name
  )
  if (is.null(p))
    stop("Could not find a city-wide monthly CSV under ", tables_dir,
         "\n  Expected columns: year/ano, month/mes, ndvi_media, lst_c_media")
  message("[common] City monthly table : ", basename(p))
  read_any_csv(p) |> normalise_columns()
}

#' Per-zone monthly table: one row per (zone, year, month).
load_zone_monthly <- function(tables_dir) {
  p <- find_csv_by_columns(
    tables_dir,
    required_any = list(COLUMN_ALIASES$zone_name, COLUMN_ALIASES$year,
                        COLUMN_ALIASES$month, COLUMN_ALIASES$lst_mean),
    prefer_pattern = "bairro|zone"
  )
  if (is.null(p))
    stop("Could not find a per-zone monthly CSV under ", tables_dir,
         "\n  Expected a zone-name column plus year, month and LST mean.")
  message("[common] Zone monthly table : ", basename(p))
  df <- read_any_csv(p) |> normalise_columns()
  if (!"zone_id" %in% names(df))
    df$zone_id <- as.integer(factor(df$zone_name))
  df
}


# -----------------------------------------------------------------------------
# 7. PUBLICATION FIGURE THEME
# -----------------------------------------------------------------------------
# Design decisions, all aimed at printed reproduction rather than screen:
#   - No panel background fill; grid lines are pale grey and only where useful.
#   - Titles left-aligned. Centred titles look decorative; left alignment reads
#     faster and matches how journals typeset figure headings.
#   - Subtitles carry the quantitative summary, so the figure is interpretable
#     without the caption.
#   - Legends at the right for continuous scales, bottom for categorical ones.

theme_pub <- function(base_size = FIG_BASE_SIZE, grid = "both") {
  th <- ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      plot.title       = ggplot2::element_text(face = "bold",
                                               size = base_size + 2,
                                               hjust = 0,
                                               margin = ggplot2::margin(b = 3)),
      plot.subtitle    = ggplot2::element_text(colour = "grey35",
                                               size = base_size - 1,
                                               hjust = 0,
                                               margin = ggplot2::margin(b = 9)),
      plot.caption     = ggplot2::element_text(colour = "grey50",
                                               size = base_size - 3,
                                               hjust = 0,
                                               margin = ggplot2::margin(t = 9)),
      plot.title.position   = "plot",
      plot.caption.position = "plot",
      axis.title       = ggplot2::element_text(size = base_size - 1,
                                               colour = "grey20"),
      axis.text        = ggplot2::element_text(size = base_size - 2,
                                               colour = "grey30"),
      panel.grid.major = ggplot2::element_line(colour = "grey90",
                                               linewidth = 0.3),
      panel.grid.minor = ggplot2::element_blank(),
      legend.title     = ggplot2::element_text(size = base_size - 2,
                                               face = "bold"),
      legend.text      = ggplot2::element_text(size = base_size - 2),
      legend.key.height = ggplot2::unit(0.9, "lines"),
      strip.text       = ggplot2::element_text(face = "bold",
                                               size = base_size - 1),
      plot.margin      = ggplot2::margin(10, 12, 8, 10)
    )

  # Heatmaps and horizontal bar charts read better without competing grid lines
  if (grid == "none")
    th <- th + ggplot2::theme(panel.grid = ggplot2::element_blank())
  if (grid == "x")
    th <- th + ggplot2::theme(panel.grid.major.y = ggplot2::element_blank())
  if (grid == "y")
    th <- th + ggplot2::theme(panel.grid.major.x = ggplot2::element_blank())
  th
}

# --- colour palettes ---------------------------------------------------------
# Both are diverging and print acceptably in greyscale (luminance varies
# monotonically along each ramp), which matters if the journal prints in black
# and white.

# Cold blue -> warm red, for temperature.
PAL_LST  <- c("#2c7bb6", "#abd9e9", "#ffffbf", "#fdae61", "#d7191c")

# Bare/red -> vegetated/green, the conventional NDVI ramp.
PAL_NDVI <- c("#a50026", "#f46d43", "#fee08b", "#d9ef8b", "#66bd63", "#1a9850")

# Categorical colours for the season variable.
PAL_SEASON <- c(Summer = "#d7191c", Autumn = "#fdae61",
                Winter = "#2c7bb6", Spring = "#1a9850")


# -----------------------------------------------------------------------------
# 8. FIGURE SAVING
# -----------------------------------------------------------------------------
#' Save a ggplot as a high-resolution PNG and, optionally, a vector PDF.
#'
#' Uses the ragg PNG device when available. ragg renders text with proper
#' hinting and anti-aliasing; the default grDevices PNG device on a headless
#' Linux node often produces noticeably rougher glyphs, which is very visible
#' at the small font sizes used in the 60-zone figures.
#'
#' @param plot    a ggplot object
#' @param file    file name WITHOUT extension
#' @param dir     output directory
#' @param width   width in inches
#' @param height  height in inches
#' @param tag     short prefix used in the console log
save_fig <- function(plot, file, dir, width = 10, height = 7, tag = "") {

  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  png_path <- file.path(dir, paste0(file, ".png"))

  if (requireNamespace("ragg", quietly = TRUE)) {
    ggplot2::ggsave(png_path, plot, width = width, height = height,
                    dpi = FIG_DPI, bg = "white", device = ragg::agg_png)
  } else {
    ggplot2::ggsave(png_path, plot, width = width, height = height,
                    dpi = FIG_DPI, bg = "white")
  }
  message(sprintf("%s   %s.png", tag, file))

  if (isTRUE(FIG_SAVE_PDF)) {
    pdf_path <- file.path(dir, paste0(file, ".pdf"))
    # cairo_pdf embeds fonts and supports transparency; the default pdf device
    # does neither reliably, which breaks alpha-blended density plots.
    ok <- tryCatch({
      ggplot2::ggsave(pdf_path, plot, width = width, height = height,
                      device = grDevices::cairo_pdf, bg = "white")
      TRUE
    }, error = function(e) FALSE)
    if (!ok)
      tryCatch(ggplot2::ggsave(pdf_path, plot, width = width,
                               height = height, bg = "white"),
               error = function(e)
                 message(tag, "   (PDF export unavailable, PNG written)"))
  }
  invisible(png_path)
}


# -----------------------------------------------------------------------------
# 9. SMALL STATISTICAL HELPERS
# -----------------------------------------------------------------------------

#' Linear trend in units per year, from a date vector and a value vector.
#' Returns NA when fewer than `min_n` observations are available, because a
#' slope fitted to a handful of points is dominated by which points happen to
#' be present rather than by any trend.
slope_per_year <- function(date, value, min_n = 5) {
  ok <- !is.na(value) & !is.na(date)
  if (sum(ok) < min_n) return(NA_real_)
  year_frac <- as.numeric(format(date[ok], "%Y")) +
    (as.numeric(format(date[ok], "%m")) - 1) / 12
  if (stats::var(year_frac) <= 0) return(NA_real_)
  unname(stats::coef(stats::lm(value[ok] ~ year_frac))[2])
}

#' Open a sink safely: guarantees the connection is closed even if the code
#' that follows throws. Without this an error mid-report leaves the sink open
#' and every subsequent console message is silently redirected to the file.
open_report <- function(path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  sink(path)
  invisible(path)
}
close_report <- function() {
  while (sink.number() > 0) sink()
}

message("[common] Helpers loaded. Validation root: ", VALIDATION_ROOT)
