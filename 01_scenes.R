# =============================================================================
# MODULE 01 - SCENE DISCOVERY AND METADATA READING
# =============================================================================
# Responsibility:
#   - Walk the DIRS_LANDSAT folders and list every available scene
#   - Parse the scene ID to obtain sensor, acquisition date and path/row
#   - Read per-scene calibration factors from the MTL file
#   - Filter by period (YEAR_START..YEAR_END), month and collection tier
#
# NOTHING HERE IS AREA-SPECIFIC - this module works for any region as long as
# the folder layout follows the USGS EarthExplorer convention described in
# 00_config.R, Section 4.
# =============================================================================

source("00_config.R")

# -----------------------------------------------------------------------------
# 1. SCENE LISTING
# -----------------------------------------------------------------------------

#' Scan all Landsat folders and return one row per scene.
#'
#' Scene ID anatomy (Collection 2 Level-2):
#'   LC08_L2SP_219076_20200115_20200823_02_T1
#'   |    |    |      |        |        |  |
#'   |    |    |      |        |        |  +-- tier (T1 / T2 / RT)
#'   |    |    |      |        |        +----- collection number (02)
#'   |    |    |      |        +-------------- processing date
#'   |    |    |      +----------------------- ACQUISITION date (what we use)
#'   |    |    +------------------------------ path (219) + row (076)
#'   |    +----------------------------------- product level
#'   +---------------------------------------- sensor prefix
#'
#' @return tibble: scene_id, sensor_prefix, sensor, acq_date, year, month,
#'                 pathrow, tier, scene_dir, group
discover_scenes <- function() {

  all_scenes <- list()
  n_skipped_period <- 0L
  n_skipped_tier   <- 0L
  n_skipped_month  <- 0L

  for (group in names(DIRS_LANDSAT)) {
    base_dir <- DIRS_LANDSAT[[group]]

    if (!dir.exists(base_dir)) {
      log_msg(sprintf("Folder not found, skipping: %s", base_dir), "WARN")
      next
    }

    # Each immediate subdirectory is one scene
    subdirs <- list.dirs(base_dir, recursive = FALSE, full.names = TRUE)

    for (d in subdirs) {
      scene_id <- basename(d)

      # Reject anything that is not a Collection 2 L2SP scene folder
      if (!grepl(SCENE_ID_PATTERN, scene_id)) {
        log_msg(sprintf("Ignored (unexpected name pattern): %s", scene_id), "WARN")
        next
      }

      parts         <- strsplit(scene_id, "_")[[1]]
      sensor_prefix <- parts[1]
      pathrow       <- parts[3]
      date_str      <- parts[4]           # YYYYMMDD, acquisition
      tier          <- if (length(parts) >= 7) parts[7] else NA_character_

      acq_date <- tryCatch(as.Date(date_str, format = "%Y%m%d"),
                           error = function(e) NA)
      if (is.na(acq_date)) {
        log_msg(sprintf("Invalid date in scene ID: %s", scene_id), "WARN")
        next
      }

      yr <- lubridate::year(acq_date)
      mo <- lubridate::month(acq_date)

      # ### Period filter (YEAR_START / YEAR_END in 00_config.R) ###
      if (yr < YEAR_START || yr > YEAR_END) {
        n_skipped_period <- n_skipped_period + 1L
        next
      }

      # ### Optional month filter (MONTHS_KEEP in 00_config.R) ###
      if (!is.null(MONTHS_KEEP) && !(mo %in% MONTHS_KEEP)) {
        n_skipped_month <- n_skipped_month + 1L
        next
      }

      # ### Optional tier filter (TIER_KEEP in 00_config.R) ###
      # Mixing T1 and T2 in one time series introduces geolocation jitter.
      if (!is.null(TIER_KEEP) && !isTRUE(tier %in% TIER_KEEP)) {
        n_skipped_tier <- n_skipped_tier + 1L
        next
      }

      sensor_info <- SENSOR_MAP[[sensor_prefix]]
      if (is.null(sensor_info)) {
        log_msg(sprintf("Unknown sensor prefix: %s", sensor_prefix), "WARN")
        next
      }

      all_scenes[[length(all_scenes) + 1]] <- tibble::tibble(
        scene_id      = scene_id,
        sensor_prefix = sensor_prefix,
        sensor        = sensor_info$sensor,
        acq_date      = acq_date,
        year          = yr,
        month         = mo,
        pathrow       = pathrow,
        tier          = tier,
        scene_dir     = d,
        group         = group
      )
    }
  }

  if (length(all_scenes) == 0) {
    stop(paste0(
      "[01_scenes] No scenes found!\n",
      "  Checked folders:\n    ",
      paste(unlist(DIRS_LANDSAT), collapse = "\n    "), "\n",
      "  Common causes:\n",
      "    - DIRS_LANDSAT paths are wrong (00_config.R Section 4)\n",
      "    - Scene folders are nested one level deeper than expected\n",
      "      (this pipeline expects: <folder>/<SCENE_ID>/<band files>)\n",
      "    - The period YEAR_START..YEAR_END excludes everything\n",
      "    - TIER_KEEP excludes everything (try TIER_KEEP <- NULL)\n"))
  }

  result <- dplyr::bind_rows(all_scenes)
  result <- result[order(result$acq_date), ]

  log_msg(sprintf("[01_scenes] %d scenes retained (%d-%d).",
                  nrow(result), YEAR_START, YEAR_END))
  if (n_skipped_period > 0)
    log_msg(sprintf("[01_scenes] %d scenes outside the period.", n_skipped_period))
  if (n_skipped_month > 0)
    log_msg(sprintf("[01_scenes] %d scenes filtered out by MONTHS_KEEP.", n_skipped_month))
  if (n_skipped_tier > 0)
    log_msg(sprintf("[01_scenes] %d scenes filtered out by TIER_KEEP ('%s').",
                    n_skipped_tier, paste(TIER_KEEP, collapse = "/")))

  result
}

# -----------------------------------------------------------------------------
# 2. MTL METADATA READING
# -----------------------------------------------------------------------------
# Why read the MTL at all, when the scaling factors are constant in
# Collection 2? Because (a) it lets the pipeline detect a scene reprocessed
# with different factors, and (b) it supplies cloud cover and sun elevation
# used later as diagnostics. The constants in 00_config.R are the fallback.

#' Read MTL metadata for one scene. Tries _MTL.txt, then _MTL.json.
#'
#' @return list(sr_scale, sr_offset, st_scale, st_offset, sun_elevation,
#'              cloud_cover, station_id, mtl_source) or NULL if unreadable
read_mtl <- function(scene_dir, scene_id) {

  all_files <- list.files(scene_dir, full.names = TRUE)

  # --- Attempt 1: plain-text MTL (the usual case) ---
  mtl_txt <- all_files[grepl("_MTL\\.txt$", all_files, ignore.case = TRUE)]
  if (length(mtl_txt) > 0) {
    res <- tryCatch(read_mtl_txt(mtl_txt[1], scene_id), error = function(e) NULL)
    if (!is.null(res)) { res$mtl_source <- "txt"; return(res) }
  }

  # --- Attempt 2: JSON MTL (newer USGS deliveries) ---
  mtl_json <- all_files[grepl("_MTL\\.json$", all_files, ignore.case = TRUE)]
  if (length(mtl_json) > 0) {
    res <- tryCatch(read_mtl_json(mtl_json[1], scene_id), error = function(e) NULL)
    if (!is.null(res)) {
      res$mtl_source <- "json"
      log_msg(sprintf("[01_scenes] MTL read via JSON: %s", scene_id), "WARN")
      return(res)
    }
  }

  found <- basename(all_files[grepl("_MTL\\.", all_files, ignore.case = TRUE)])
  log_msg(sprintf("[01_scenes] MTL unreadable for %s | MTL files present: [%s]",
                  scene_id,
                  if (length(found) == 0) "none" else paste(found, collapse = ", ")),
          "WARN")
  NULL
}


#' Parse the key = value text MTL.
#'
#' CRITICAL DETAIL: keys such as REFLECTANCE_MULT_BAND_4 appear in MULTIPLE
#' GROUP sections of the MTL (Level-1 and Level-2) with DIFFERENT values.
#' We must always take the Level-2 values, so every lookup is scoped to a named
#' GROUP block rather than grepping the whole file.
read_mtl_txt <- function(filepath, scene_id) {
  lines <- readLines(filepath, warn = FALSE)

  # Numeric lookup restricted to one GROUP = <section> ... END_GROUP block
  get_val_section <- function(section_name, key) {
    ini <- which(grepl(paste0("GROUP\\s*=\\s*", section_name), lines))[1]
    fim <- which(grepl(paste0("END_GROUP\\s*=\\s*", section_name), lines))[1]
    if (is.na(ini) || is.na(fim)) return(NA_real_)
    block <- lines[ini:fim]
    ln    <- block[grepl(paste0("^\\s*", key, "\\s*="), block)]
    if (length(ln) == 0) return(NA_real_)
    suppressWarnings(as.numeric(trimws(sub(".*=\\s*", "", ln[1]))))
  }

  get_str_section <- function(section_name, key) {
    ini <- which(grepl(paste0("GROUP\\s*=\\s*", section_name), lines))[1]
    fim <- which(grepl(paste0("END_GROUP\\s*=\\s*", section_name), lines))[1]
    if (is.na(ini) || is.na(fim)) return(NA_character_)
    block <- lines[ini:fim]
    ln    <- block[grepl(paste0("^\\s*", key, "\\s*="), block)]
    if (length(ln) == 0) return(NA_character_)
    trimws(gsub('"', '', sub(".*=\\s*", "", ln[1])))
  }

  # --- Surface reflectance factors (Level-2 section) ---
  # L8/9 use BAND_4 as Red; L4-7 use BAND_3. Both carry the SAME Collection 2
  # scaling factor, so either serves: try BAND_3 first, then BAND_4.
  sr_mult <- get_val_section("LEVEL2_SURFACE_REFLECTANCE_PARAMETERS", "REFLECTANCE_MULT_BAND_3")
  sr_add  <- get_val_section("LEVEL2_SURFACE_REFLECTANCE_PARAMETERS", "REFLECTANCE_ADD_BAND_3")
  if (is.na(sr_mult)) sr_mult <- get_val_section("LEVEL2_SURFACE_REFLECTANCE_PARAMETERS", "REFLECTANCE_MULT_BAND_4")
  if (is.na(sr_add))  sr_add  <- get_val_section("LEVEL2_SURFACE_REFLECTANCE_PARAMETERS", "REFLECTANCE_ADD_BAND_4")

  # --- Surface temperature factors (B6 = L4-7, B10 = L8/9) ---
  st_mult <- get_val_section("LEVEL2_SURFACE_TEMPERATURE_PARAMETERS", "TEMPERATURE_MULT_BAND_ST_B6")
  st_add  <- get_val_section("LEVEL2_SURFACE_TEMPERATURE_PARAMETERS", "TEMPERATURE_ADD_BAND_ST_B6")
  if (is.na(st_mult)) st_mult <- get_val_section("LEVEL2_SURFACE_TEMPERATURE_PARAMETERS", "TEMPERATURE_MULT_BAND_ST_B10")
  if (is.na(st_add))  st_add  <- get_val_section("LEVEL2_SURFACE_TEMPERATURE_PARAMETERS", "TEMPERATURE_ADD_BAND_ST_B10")

  # --- Image attributes (diagnostics only) ---
  sun_elev    <- get_val_section("IMAGE_ATTRIBUTES", "SUN_ELEVATION")
  cloud_cover <- get_val_section("IMAGE_ATTRIBUTES", "CLOUD_COVER")
  station_id  <- get_str_section("IMAGE_ATTRIBUTES", "STATION_ID")

  verify_scaling(sr_mult, sr_add, st_mult, st_add, scene_id)

  list(
    sr_scale      = if (!is.na(sr_mult)) sr_mult else SR_SCALE,
    sr_offset     = if (!is.na(sr_add))  sr_add  else SR_OFFSET,
    st_scale      = if (!is.na(st_mult)) st_mult else ST_SCALE,
    st_offset     = if (!is.na(st_add))  st_add  else ST_OFFSET,
    sun_elevation = sun_elev,
    cloud_cover   = cloud_cover,
    station_id    = station_id
  )
}


#' Parse the nested JSON MTL delivered by newer USGS orders.
read_mtl_json <- function(filepath, scene_id) {
  j    <- jsonlite::read_json(filepath)
  root <- if (!is.null(j$LANDSAT_METADATA_FILE)) j$LANDSAT_METADATA_FILE else j

  jget <- function(...) {
    val <- root
    for (k in c(...)) { val <- val[[k]]; if (is.null(val)) return(NA_real_) }
    suppressWarnings(as.numeric(val))
  }
  jget_str <- function(...) {
    val <- root
    for (k in c(...)) { val <- val[[k]]; if (is.null(val)) return(NA_character_) }
    as.character(val)
  }

  sr_sections <- c("LEVEL2_SURFACE_REFLECTANCE_PARAMETERS", "LEVEL2_PROCESSING_RECORD")
  sr_mult <- NA_real_; sr_add <- NA_real_
  for (sec in sr_sections) {
    if (is.na(sr_mult)) sr_mult <- jget(sec, "REFLECTANCE_MULT_BAND_4")
    if (is.na(sr_mult)) sr_mult <- jget(sec, "REFLECTANCE_MULT_BAND_3")
    if (is.na(sr_add))  sr_add  <- jget(sec, "REFLECTANCE_ADD_BAND_4")
    if (is.na(sr_add))  sr_add  <- jget(sec, "REFLECTANCE_ADD_BAND_3")
  }

  st_mult <- jget("LEVEL2_SURFACE_TEMPERATURE_PARAMETERS", "TEMPERATURE_MULT_BAND_ST_B10")
  st_add  <- jget("LEVEL2_SURFACE_TEMPERATURE_PARAMETERS", "TEMPERATURE_ADD_BAND_ST_B10")
  if (is.na(st_mult)) st_mult <- jget("LEVEL2_SURFACE_TEMPERATURE_PARAMETERS", "TEMPERATURE_MULT_BAND_ST_B6")
  if (is.na(st_add))  st_add  <- jget("LEVEL2_SURFACE_TEMPERATURE_PARAMETERS", "TEMPERATURE_ADD_BAND_ST_B6")

  verify_scaling(sr_mult, sr_add, st_mult, st_add, scene_id)

  list(
    sr_scale      = if (!is.na(sr_mult)) sr_mult else SR_SCALE,
    sr_offset     = if (!is.na(sr_add))  sr_add  else SR_OFFSET,
    st_scale      = if (!is.na(st_mult)) st_mult else ST_SCALE,
    st_offset     = if (!is.na(st_add))  st_add  else ST_OFFSET,
    sun_elevation = jget("IMAGE_ATTRIBUTES", "SUN_ELEVATION"),
    cloud_cover   = jget("IMAGE_ATTRIBUTES", "CLOUD_COVER"),
    station_id    = jget_str("IMAGE_ATTRIBUTES", "STATION_ID")
  )
}


#' Warn when a scene's MTL reports factors different from the Collection 2
#' defaults. This is rare and usually means a reprocessed or non-standard
#' product - worth knowing before it silently shifts your time series.
verify_scaling <- function(sr_mult, sr_add, st_mult, st_add, scene_id) {
  chk <- function(val, expected, name) {
    if (!is.na(val) && abs(val - expected) > 1e-8)
      log_msg(sprintf("NON-STANDARD SCALING in %s: %s=%.10f (expected %.10f)",
                      scene_id, name, val, expected), "WARN")
  }
  chk(sr_mult, SR_SCALE,  "SR_MULT")
  chk(sr_add,  SR_OFFSET, "SR_ADD")
  chk(st_mult, ST_SCALE,  "ST_MULT")
  chk(st_add,  ST_OFFSET, "ST_ADD")
}

# -----------------------------------------------------------------------------
# 3. FULL INVENTORY
# -----------------------------------------------------------------------------

#' Safely pull a field out of an MTL list, with fallback.
#' (A plain ifelse() would fail here because it evaluates both branches.)
mtl_field <- function(mtl_list, field, fallback) {
  if (is.null(mtl_list)) return(fallback)
  val <- mtl_list[[field]]
  if (is.null(val) || length(val) == 0 || is.na(val)) return(fallback)
  val
}


#' Combine discover_scenes() + read_mtl() + band file lookup into one table.
#'
#' The `complete` column marks scenes that have all four required rasters
#' (Red, NIR, thermal, QA). Only complete scenes are processed downstream.
build_inventory <- function() {

  scenes <- discover_scenes()
  rows   <- vector("list", nrow(scenes))

  for (i in seq_len(nrow(scenes))) {
    s <- scenes[i, ]

    file_qa  <- find_band(s$scene_dir, s$scene_id, "QA_PIXEL")
    file_red <- find_band(s$scene_dir, s$scene_id, SENSOR_MAP[[s$sensor_prefix]]$red)
    file_nir <- find_band(s$scene_dir, s$scene_id, SENSOR_MAP[[s$sensor_prefix]]$nir)
    file_lst <- find_band(s$scene_dir, s$scene_id, SENSOR_MAP[[s$sensor_prefix]]$lst)

    mtl <- read_mtl(s$scene_dir, s$scene_id)

    rows[[i]] <- tibble::tibble(
      scene_id      = s$scene_id,
      sensor_prefix = s$sensor_prefix,
      sensor        = s$sensor,
      acq_date      = s$acq_date,
      year          = s$year,
      month         = s$month,
      pathrow       = s$pathrow,
      tier          = s$tier,
      scene_dir     = s$scene_dir,
      group         = s$group,
      file_qa       = file_qa,
      file_red      = file_red,
      file_nir      = file_nir,
      file_lst      = file_lst,
      # Per-scene MTL values take priority; config constants are the fallback
      sr_scale      = mtl_field(mtl, "sr_scale",      SR_SCALE),
      sr_offset     = mtl_field(mtl, "sr_offset",     SR_OFFSET),
      st_scale      = mtl_field(mtl, "st_scale",      ST_SCALE),
      st_offset     = mtl_field(mtl, "st_offset",     ST_OFFSET),
      cloud_cover   = mtl_field(mtl, "cloud_cover",   NA_real_),
      sun_elevation = mtl_field(mtl, "sun_elevation", NA_real_),
      mtl_ok        = !is.null(mtl),
      complete      = !is.na(file_qa) & !is.na(file_red) &
                      !is.na(file_nir) & !is.na(file_lst)
    )
  }

  out <- dplyr::bind_rows(rows)

  n_complete   <- sum(out$complete)
  n_incomplete <- sum(!out$complete)
  n_no_mtl     <- sum(!out$mtl_ok)

  log_msg(sprintf("[01_scenes] Inventory: %d complete, %d incomplete.",
                  n_complete, n_incomplete))
  if (n_no_mtl > 0)
    log_msg(sprintf("[01_scenes] %d scenes without MTL (using Collection 2 default factors).",
                    n_no_mtl), "WARN")
  if (n_incomplete > 0) {
    inc <- out$scene_id[!out$complete]
    log_msg(sprintf("Incomplete scenes (missing bands): %s",
                    paste(inc, collapse = ", ")), "WARN")
  }

  out
}

# -----------------------------------------------------------------------------
# INTERNAL HELPER
# -----------------------------------------------------------------------------

#' Locate one band file inside a scene folder.
#' @param band_tag e.g. "SR_B4", "ST_B10", "QA_PIXEL"
find_band <- function(scene_dir, scene_id, band_tag) {
  pattern <- paste0(scene_id, "_", band_tag, "\\.TIF$")
  f <- list.files(scene_dir, pattern = pattern,
                  full.names = TRUE, ignore.case = TRUE)
  if (length(f) == 0) return(NA_character_)
  f[1]
}
