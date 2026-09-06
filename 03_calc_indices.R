# =============================================================================
# MODULE 03 - NDVI AND LST COMPUTATION
# =============================================================================
# Responsibility:
#   - Convert optical band DNs to Surface Reflectance (SR)
#   - Compute NDVI = (NIR - Red) / (NIR + Red)
#   - Convert the thermal ST band from DN to Kelvin, then to Celsius
#
# THE SCIENCE HERE IS AREA-INDEPENDENT. The only area-dependent inputs are
# LST_MIN / LST_MAX (climate bounds), set in 00_config.R Section 6.
#
# Formulas (USGS Collection 2 Level-2, LSDS-1619):
#   SR      = DN * 0.0000275 + (-0.2)
#   ST(K)   = DN * 0.00341802 + 149.0
#   LST(C)  = ST(K) - 273.15
#   NDVI    = (NIR - Red) / (NIR + Red)          [Rouse et al. 1974]
# =============================================================================

source("00_config.R")
source("02_qa_mask.R")

# -----------------------------------------------------------------------------
# 1. SURFACE REFLECTANCE CALIBRATION
# -----------------------------------------------------------------------------

#' Convert a DN band to surface reflectance in [0, 1].
#'
#' @param file_band Path to the SR band GeoTIFF
#' @param mask      QA mask from make_qa_mask()
#' @param scale     Multiplicative factor (per-scene value from the MTL)
#' @param offset    Additive offset (per-scene value from the MTL)
#' @return SpatRaster of surface reflectance
calibrate_sr <- function(file_band, mask,
                         scale = SR_SCALE, offset = SR_OFFSET) {

  if (is.na(file_band) || !file.exists(file_band))
    stop("[03_calc] Band file not found: ", file_band)

  b <- terra::rast(file_band)

  # Remove the fill DN before scaling.
  # DN = 0 is the Collection 2 fill value. Scaled it would become exactly the
  # offset (-0.2), a physically impossible reflectance that would otherwise
  # only be caught by the clamp below. Removing it explicitly is clearer and
  # also protects the NDVI denominator from a degenerate (nir + red) value.
  b[b == 0] <- NA

  b <- apply_qa_mask(b, mask)

  sr <- b * scale + offset

  # Enforce the physical range. values = FALSE means OUT-OF-RANGE PIXELS
  # BECOME NA (discarded).
  # NEVER use values = TRUE / NA here: that would SATURATE pixels to exactly
  # 0 or 1 instead of dropping them, silently biasing every downstream mean.
  terra::clamp(sr, lower = 0, upper = 1, values = FALSE)
}

# -----------------------------------------------------------------------------
# 2. NDVI
# -----------------------------------------------------------------------------

#' Compute NDVI from the Red and NIR band files.
#'
#' NDVI = (NIR - Red) / (NIR + Red), Rouse et al. (1974).
#' Interpretation: < 0 water/snow; 0-0.2 bare soil, rock, built-up;
#' 0.2-0.5 sparse vegetation, grass; > 0.5 dense healthy vegetation.
#'
#' Band numbers differ by sensor (B3/B4 for L4-7, B4/B5 for L8/9) — that
#' mapping lives in SENSOR_MAP (00_config.R) and is resolved by the caller.
#'
#' @return SpatRaster named "NDVI"
calc_ndvi <- function(file_red, file_nir, mask,
                      sr_scale = SR_SCALE, sr_offset = SR_OFFSET) {

  red <- calibrate_sr(file_red, mask, sr_scale, sr_offset)
  nir <- calibrate_sr(file_nir, mask, sr_scale, sr_offset)

  # Bilinear is acceptable here because reflectance is a continuous quantity.
  # In Collection 2 both bands are already on the same grid, so this rarely runs.
  if (!terra::compareGeom(red, nir, stopOnError = FALSE))
    nir <- terra::resample(nir, red, method = "bilinear")

  ndvi <- (nir - red) / (nir + red)

  # See the warning in calibrate_sr(): values = FALSE discards, it does not saturate.
  ndvi <- terra::clamp(ndvi, lower = NDVI_MIN, upper = NDVI_MAX, values = FALSE)

  names(ndvi) <- "NDVI"
  ndvi
}

# -----------------------------------------------------------------------------
# 3. LAND SURFACE TEMPERATURE
# -----------------------------------------------------------------------------

#' Convert the thermal ST band from DN to degrees Celsius.
#'
#' Two-stage rejection of impossible values:
#'
#'   (a) SENTINEL-DN FILTER, applied in DN space BEFORE scaling. Fill and
#'       near-fill DNs scale to absurd temperatures (e.g. -120 C). We compute
#'       the DN corresponding to LST_MIN and drop everything below it.
#'
#'       The threshold is DERIVED FROM THE ACTUAL PER-SCENE st_scale/st_offset,
#'       inverting the calibration equation:
#'           LST_MIN = DN * st_scale + st_offset - 273.15
#'       =>  DN_min  = (LST_MIN + 273.15 - st_offset) / st_scale
#'
#'       (Earlier versions of this pipeline hard-coded the Collection 2
#'       constants in this formula. That silently produced a WRONG cut-off for
#'       any scene whose MTL reported different factors, and broke entirely if
#'       a user edited ST_OFFSET. Always derive it, never hard-code it.)
#'
#'   (b) PHYSICAL CLAMP in Celsius, enforcing the LST_MAX ceiling.
#'
#' @param file_lst  Path to the ST band (ST_B6 or ST_B10)
#' @param mask      QA mask
#' @param st_scale  Per-scene multiplicative factor (from MTL)
#' @param st_offset Per-scene additive offset in Kelvin (from MTL)
#' @return SpatRaster named "LST_Celsius"
calc_lst <- function(file_lst, mask,
                     st_scale = ST_SCALE, st_offset = ST_OFFSET) {

  if (is.na(file_lst) || !file.exists(file_lst))
    stop("[03_calc] Thermal band file not found: ", file_lst)

  st_raw <- terra::rast(file_lst)
  st_raw <- apply_qa_mask(st_raw, mask)

  # (a) Sentinel-DN filter, derived from THIS scene's calibration factors
  dn_min_lst <- (LST_MIN + KELVIN_TO_CELSIUS - st_offset) / st_scale
  st_raw[st_raw < dn_min_lst] <- NA

  # Conversion: DN -> Kelvin -> Celsius
  st_kelvin   <- st_raw * st_scale + st_offset
  lst_celsius <- st_kelvin - KELVIN_TO_CELSIUS

  # (b) Physical clamp. values = FALSE -> out-of-range pixels become NA.
  # Using values = TRUE/NA here would pin artefacts to exactly LST_MIN or
  # LST_MAX and pull every zone mean toward the bound.
  lst_celsius <- terra::clamp(lst_celsius,
                              lower  = LST_MIN,
                              upper  = LST_MAX,
                              values = FALSE)

  names(lst_celsius) <- "LST_Celsius"
  lst_celsius
}

# -----------------------------------------------------------------------------
# HELPER: RESOLUTION CHECK
# -----------------------------------------------------------------------------

#' Warn (but do not stop) if a band is not at the expected resolution.
#'
#' USGS delivers all Collection 2 Level-2 bands resampled to 30 m, including
#' the thermal bands (natively 100-120 m, oversampled to 30 m). A deviation
#' here means a non-standard product.
check_resolution <- function(r, name, scene_id, expected_res = 30) {

  res_now <- terra::res(r)
  res_x   <- round(res_now[1], 2)
  res_y   <- round(res_now[2], 2)

  if (abs(res_x - expected_res) > 0.5 || abs(res_y - expected_res) > 0.5) {
    log_msg(sprintf(
      "[03_calc] UNEXPECTED RESOLUTION in %s | band %s: %.1f x %.1f m (expected %dm)",
      scene_id, name, res_x, res_y, expected_res), "WARN")
  }
}

# -----------------------------------------------------------------------------
# 4. FULL SCENE PROCESSING
# -----------------------------------------------------------------------------

#' Process one scene end-to-end: QA mask, NDVI, LST.
#'
#' @param scene One row of the inventory tibble from 01_scenes.R
#' @return list(ndvi, lst, pct_valid) or NULL if the scene is rejected
process_scene <- function(scene) {

  sid <- scene$scene_id

  if (!scene$complete) {
    log_msg(sprintf("[03_calc] Incomplete scene skipped: %s", sid), "WARN")
    return(NULL)
  }

  tryCatch({
    log_msg(sprintf("[03_calc] Processing: %s", sid))

    # 0. Resolution sanity check (opens headers only, does not read pixels)
    check_resolution(terra::rast(scene$file_red), basename(scene$file_red), sid)
    check_resolution(terra::rast(scene$file_nir), basename(scene$file_nir), sid)
    check_resolution(terra::rast(scene$file_lst), basename(scene$file_lst), sid)
    check_resolution(terra::rast(scene$file_qa),  "QA_PIXEL",               sid)

    # 1. Quality mask
    mask  <- make_qa_mask(scene$file_qa)
    stats <- qa_stats(mask, sid)

    # Early rejection of heavily clouded scenes.
    # ### Threshold: MIN_SCENE_VALID_FRAC in 00_config.R Section 8 ###
    if (stats$pct_valid < MIN_SCENE_VALID_FRAC * 100) {
      log_msg(sprintf("[03_calc] Cloud cover too high (%.1f%% valid): %s — skipping.",
                      stats$pct_valid, sid), "WARN")
      return(NULL)
    }

    # 2. NDVI — uses this scene's own MTL scaling factors
    ndvi <- calc_ndvi(scene$file_red, scene$file_nir, mask,
                      scene$sr_scale, scene$sr_offset)

    # 3. LST — likewise
    lst <- calc_lst(scene$file_lst, mask,
                    scene$st_scale, scene$st_offset)

    log_msg(sprintf(
      "[03_calc] OK: %s | NDVI [%.3f, %.3f] | LST [%.1f, %.1f] C",
      sid,
      terra::global(ndvi, "min", na.rm = TRUE)[[1]],
      terra::global(ndvi, "max", na.rm = TRUE)[[1]],
      terra::global(lst,  "min", na.rm = TRUE)[[1]],
      terra::global(lst,  "max", na.rm = TRUE)[[1]]
    ))

    list(ndvi = ndvi, lst = lst, pct_valid = stats$pct_valid)

  }, error = function(e) {
    log_msg(sprintf("[03_calc] ERROR in %s: %s", sid, conditionMessage(e)), "ERROR")
    NULL
  })
}
