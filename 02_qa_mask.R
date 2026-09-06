# =============================================================================
# MODULE 02 - QUALITY MASK (QA_PIXEL)
# =============================================================================
# Responsibility:
#   - Decode the 16-bit QA_PIXEL band of Landsat Collection 2
#   - Build a binary mask: 1 = usable pixel, NA = discard
#   - Apply that mask to any band raster
#
# WHICH BITS ARE MASKED IS A CONFIGURATION DECISION, not a code decision:
# see QA_BITS_MASK in 00_config.R Section 7, where each bit is documented
# along with when you would want to change it (the WATER bit in particular).
#
# QA_PIXEL bit reference (Collection 2, all sensors):
#   Bit 0  Fill                Bit 4  Cloud Shadow
#   Bit 1  Dilated Cloud       Bit 5  Snow / Ice
#   Bit 2  Cirrus (L8/9 only)  Bit 6  Clear  <- informational, never mask
#   Bit 3  Cloud               Bit 7  Water
# Source: USGS Landsat Collection 2 Level-2 Science Product Guide.
# =============================================================================

source("00_config.R")

# -----------------------------------------------------------------------------
# 1. MASK CONSTRUCTION
# -----------------------------------------------------------------------------

#' Build a binary quality mask from a QA_PIXEL file.
#'
#' A pixel is discarded if ANY of the requested bits is set (logical OR).
#' This is deliberately conservative: for time-series work, a false negative
#' (keeping a cloudy pixel) corrupts the mean, whereas a false positive
#' (discarding a good pixel) only costs sample size.
#'
#' @param file_qa  Path to the *_QA_PIXEL.TIF file
#' @param bits     Integer vector of bits to mask (default: QA_BITS_MASK)
#' @return SpatRaster with values 1 (keep) or NA (discard)
make_qa_mask <- function(file_qa, bits = QA_BITS_MASK) {

  if (is.na(file_qa) || !file.exists(file_qa))
    stop("[02_qa] QA_PIXEL file not found: ", file_qa)

  qa <- terra::rast(file_qa)

  # Accumulate a "bad pixel" mask, starting from all-good
  bad_mask <- terra::init(qa, fun = 0)

  for (bit in bits) {
    # Bit extraction: (value %/% 2^bit) %% 2 == 1  means the bit is set
    bit_layer <- (qa %/% (2^bit)) %% 2
    bad_mask  <- bad_mask | (bit_layer == 1)
  }

  # 0 (good) -> 1 ; 1 (bad) -> NA
  terra::ifel(bad_mask == 0, 1, NA)
}

# -----------------------------------------------------------------------------
# 2. MASK STATISTICS (DIAGNOSTIC)
# -----------------------------------------------------------------------------

#' Report how much of the scene survived masking.
#'
#' NOTE ON INTERPRETATION: this is computed over the FULL SCENE, which is much
#' larger than most study areas. A scene can be 60% valid overall and still be
#' completely clouded over your AOI. That is why 06_main.R performs a SECOND,
#' AOI-specific validity check after clipping (MIN_AOI_VALID_FRAC).
#'
#' @return list(n_valid, n_masked, pct_valid)
qa_stats <- function(mask, scene_id = "") {

  freq_tab <- terra::freq(mask, value = 1)
  n_valid  <- if (nrow(freq_tab) > 0) freq_tab$count[1] else 0
  n_total  <- terra::ncell(mask)
  n_masked <- n_total - n_valid
  pct_valid <- round(100 * n_valid / n_total, 2)

  log_msg(sprintf("[02_qa] %s | Valid: %d (%.1f%%) | Masked: %d",
                  scene_id, n_valid, pct_valid, n_masked))

  list(n_valid = n_valid, n_masked = n_masked, pct_valid = pct_valid)
}

# -----------------------------------------------------------------------------
# 3. MASK APPLICATION
# -----------------------------------------------------------------------------

#' Apply a QA mask to a band raster.
#'
#' Nearest-neighbour resampling is used if geometries differ, because the mask
#' is CATEGORICAL - interpolating quality flags would invent meaningless
#' intermediate values. In Collection 2 all bands share the same 30 m grid, so
#' this branch is a safety net rather than the normal path.
#'
#' @return SpatRaster with invalid pixels set to NA
apply_qa_mask <- function(band_raster, mask) {

  if (!terra::compareGeom(band_raster, mask, stopOnError = FALSE)) {
    mask <- terra::resample(mask, band_raster, method = "near")
  }

  # Multiplying by a 1/NA mask propagates NA into the discarded pixels
  band_raster * mask
}
