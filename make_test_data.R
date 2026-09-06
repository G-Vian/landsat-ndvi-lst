# =============================================================================
# make_test_data.R - synthetic Landsat scenes + AOI, for pipeline testing
# =============================================================================
# Builds a miniature but STRUCTURALLY FAITHFUL replica of a real download:
# correct scene-ID grammar, correct band names per sensor, real QA_PIXEL bit
# packing, real MTL group structure, and DN values in the real valid ranges.
#
# It is deliberately tiny (200x200 px) so the whole pipeline runs in seconds.
# =============================================================================
suppressPackageStartupMessages({library(terra); library(sf)})
set.seed(42)   # fixed seed: the test data is reproducible

# Output root. run_tests.R sets LANDSAT_TEST_ROOT to a temporary directory;
# when run standalone, it falls back to a folder next to this script.
ROOT <- Sys.getenv("LANDSAT_TEST_ROOT",
                   unset = file.path(tempdir(), "landsat_testdata"))
unlink(ROOT, recursive = TRUE)
dir.create(file.path(ROOT, "landsat"), recursive = TRUE)
dir.create(file.path(ROOT, "shp"), recursive = TRUE)

# --- AOI: 4 contiguous zones, EPSG:31983 (SIRGAS 2000 / UTM 23S) -------------
x0 <- 360000; y0 <- 7350000; w <- 3000; h <- 3000
mk <- function(i, j) {
  st_polygon(list(cbind(
    c(x0+i*w, x0+(i+1)*w, x0+(i+1)*w, x0+i*w, x0+i*w),
    c(y0+j*h, y0+j*h, y0+(j+1)*h, y0+(j+1)*h, y0+j*h))))
}
zones <- st_sf(
  NM_BAIRRO = c("Centro", "Gonzaga", "Ponta da Praia", "Morro Sao Bento"),
  COD       = 1:4,
  geometry  = st_sfc(mk(0,0), mk(1,0), mk(0,1), mk(1,1), crs = 31983))
st_write(zones, file.path(ROOT, "shp", "zones.gpkg"), quiet = TRUE)
cat(sprintf("AOI written: %d zones, total %.1f km2\n",
            nrow(zones), as.numeric(sum(st_area(zones)))/1e6))

# --- scene template: 200x200 px at 30 m, covering the AOI with margin -------
tmpl <- rast(xmin = x0-1500, xmax = x0+2*w+1500,
             ymin = y0-1500, ymax = y0+2*h+1500,
             resolution = 30, crs = "EPSG:31983")
cat(sprintf("Scene grid: %d x %d px\n", nrow(tmpl), ncol(tmpl)))

# --- MTL writer: mirrors the real GROUP structure ---------------------------
# Level-1 values are deliberately DIFFERENT from Level-2 values, to verify the
# parser reads from the correct group (a real trap in genuine MTL files).
write_mtl <- function(path, sid, thermal_band, cloud, sun_elev) {
  writeLines(c(
    "GROUP = LANDSAT_METADATA_FILE",
    "  GROUP = IMAGE_ATTRIBUTES",
    sprintf("    CLOUD_COVER = %.2f", cloud),
    sprintf("    SUN_ELEVATION = %.5f", sun_elev),
    "    STATION_ID = \"LGN\"",
    "  END_GROUP = IMAGE_ATTRIBUTES",
    "  GROUP = LEVEL1_RADIOMETRIC_RESCALING",
    "    REFLECTANCE_MULT_BAND_3 = 9.9999E-05",   # <- WRONG on purpose (L1)
    "    REFLECTANCE_ADD_BAND_3 = -0.099999",     # <- WRONG on purpose (L1)
    "  END_GROUP = LEVEL1_RADIOMETRIC_RESCALING",
    "  GROUP = LEVEL2_SURFACE_REFLECTANCE_PARAMETERS",
    "    REFLECTANCE_MULT_BAND_3 = 2.75e-05",
    "    REFLECTANCE_ADD_BAND_3 = -0.2",
    "    REFLECTANCE_MULT_BAND_4 = 2.75e-05",
    "    REFLECTANCE_ADD_BAND_4 = -0.2",
    "  END_GROUP = LEVEL2_SURFACE_REFLECTANCE_PARAMETERS",
    "  GROUP = LEVEL2_SURFACE_TEMPERATURE_PARAMETERS",
    sprintf("    TEMPERATURE_MULT_BAND_%s = 0.00341802", thermal_band),
    sprintf("    TEMPERATURE_ADD_BAND_%s = 149.0", thermal_band),
    "  END_GROUP = LEVEL2_SURFACE_TEMPERATURE_PARAMETERS",
    "END_GROUP = LANDSAT_METADATA_FILE",
    "END"), path)
}

# --- QA_PIXEL builder: real Collection 2 bit packing -------------------------
# bit0 fill | bit1 dilated | bit3 cloud | bit4 shadow | bit6 clear | bit7 water
make_qa <- function(tmpl, cloud_frac, water_zone = TRUE) {
  qa <- tmpl; values(qa) <- 0L
  n  <- ncell(qa); v <- rep(0L, n)
  v <- bitwOr(v, bitwShiftL(1L, 6))                       # start: all clear
  # cloud blob in a corner
  nc <- round(n * cloud_frac)
  if (nc > 0) {
    idx <- seq_len(nc)
    v[idx] <- bitwOr(bitwAnd(v[idx], bitwNot(bitwShiftL(1L,6))), bitwShiftL(1L,3))
    dil <- (nc+1):min(n, nc + round(n*0.02))
    if (length(dil) > 0 && dil[1] <= n)
      v[dil] <- bitwOr(bitwAnd(v[dil], bitwNot(bitwShiftL(1L,6))), bitwShiftL(1L,1))
  }
  # a water strip (bottom rows) -> exercises the water-mask decision
  if (water_zone) {
    wr <- (n - round(n*0.08)):n
    v[wr] <- bitwOr(bitwAnd(v[wr], bitwNot(bitwShiftL(1L,6))), bitwShiftL(1L,7))
  }
  # fill border (scene edge) -> exercises bit 0 and the DN==0 rule
  v[1:round(n*0.01)] <- bitwShiftL(1L, 0)
  values(qa) <- v
  qa
}

# --- scene generator ---------------------------------------------------------
# DN ranges follow the USGS valid ranges: SR 7273-43636, ST 293-65535.
gen_scene <- function(sensor, date_str, pathrow = "219076", tier = "T1",
                      cloud_frac = 0.12, mean_lst_c = 26, base_ndvi = 0.35) {
  cfg <- switch(sensor,
    LC08 = list(red="SR_B4", nir="SR_B5", th="ST_B10"),
    LC09 = list(red="SR_B4", nir="SR_B5", th="ST_B10"),
    LE07 = list(red="SR_B3", nir="SR_B4", th="ST_B6"),
    LT05 = list(red="SR_B3", nir="SR_B4", th="ST_B6"))
  sid <- sprintf("%s_L2SP_%s_%s_%s_02_%s", sensor, pathrow, date_str,
                 format(as.Date(date_str,"%Y%m%d")+200, "%Y%m%d"), tier)
  d <- file.path(ROOT, "landsat", sid); dir.create(d, recursive = TRUE)

  qa <- make_qa(tmpl, cloud_frac)
  n  <- ncell(tmpl)

  # spatial gradient so zones differ from each other (testable signal)
  gx <- init(tmpl, "x"); gy <- init(tmpl, "y")
  gx <- (gx - min(values(gx))) / diff(range(values(gx)))
  gy <- (gy - min(values(gy))) / diff(range(values(gy)))

  # NDVI target -> invert to red/nir DN
  ndvi_t <- base_ndvi + 0.30*values(gy) - 0.15*values(gx) +
            rnorm(n, 0, 0.03)
  ndvi_t <- pmax(pmin(ndvi_t, 0.85), -0.1)
  red_sr <- 0.06 + 0.04*values(gx) + rnorm(n, 0, 0.005)
  red_sr <- pmax(red_sr, 0.01)
  nir_sr <- red_sr * (1 + ndvi_t) / (1 - ndvi_t)
  nir_sr <- pmin(nir_sr, 0.95)

  to_dn_sr <- function(sr) as.integer(round((sr + 0.2) / 0.0000275))
  red <- tmpl; values(red) <- to_dn_sr(red_sr)
  nir <- tmpl; values(nir) <- to_dn_sr(nir_sr)

  # LST: warmer where less vegetated (realistic inverse relation)
  lst_c <- mean_lst_c + 6*(1-values(gy)) + 3*values(gx) + rnorm(n, 0, 0.7)
  st_dn <- as.integer(round((lst_c + 273.15 - 149.0) / 0.00341802))
  st  <- tmpl; values(st) <- st_dn

  # fill pixels -> DN 0 in every band (as USGS delivers)
  fillpx <- which(bitwAnd(values(qa), 1L) == 1L)
  if (length(fillpx)) {
    red[fillpx] <- 0L; nir[fillpx] <- 0L; st[fillpx] <- 0L
  }

  wopt <- list(datatype = "INT2U", gdal = "COMPRESS=LZW")
  writeRaster(red, file.path(d, sprintf("%s_%s.TIF", sid, cfg$red)), overwrite=TRUE, wopt=wopt)
  writeRaster(nir, file.path(d, sprintf("%s_%s.TIF", sid, cfg$nir)), overwrite=TRUE, wopt=wopt)
  writeRaster(st,  file.path(d, sprintf("%s_%s.TIF", sid, cfg$th )), overwrite=TRUE, wopt=wopt)
  writeRaster(qa,  file.path(d, sprintf("%s_QA_PIXEL.TIF", sid)),   overwrite=TRUE, wopt=wopt)
  write_mtl(file.path(d, sprintf("%s_MTL.txt", sid)), sid, cfg$th,
            cloud_frac*100, 45 + rnorm(1,0,5))
  sid
}

# --- build a 2-year, multi-sensor series ------------------------------------
# 2015 and 2016, mixing L7 and L8 -> exercises the sensor-composition logic.
ids <- character(0)
for (yr in c(2015, 2016)) {
  for (mo in c(2, 5, 8, 11)) {
    seas <- 4 * cos((mo - 2) / 12 * 2 * pi)      # summer warm, winter cool
    ids <- c(ids, gen_scene("LC08", sprintf("%d%02d05", yr, mo),
                            cloud_frac = runif(1, .05, .25),
                            mean_lst_c = 26 + seas,
                            base_ndvi  = 0.35 + 0.05*cos((mo-2)/12*2*pi)))
    if (mo %in% c(5, 11))
      ids <- c(ids, gen_scene("LE07", sprintf("%d%02d13", yr, mo),
                              cloud_frac = runif(1, .10, .30),
                              mean_lst_c = 25 + seas,
                              base_ndvi  = 0.33 + 0.05*cos((mo-2)/12*2*pi)))
  }
}
# one heavily clouded scene -> must be REJECTED by MIN_SCENE_VALID_FRAC
ids <- c(ids, gen_scene("LC08", "20160710", cloud_frac = 0.95))
# one Tier 2 scene -> must be REJECTED by TIER_KEEP
ids <- c(ids, gen_scene("LC08", "20150902", tier = "T2"))

cat(sprintf("\nGenerated %d scenes in %s\n", length(ids), file.path(ROOT,"landsat")))
cat("Expected: 12 usable, 1 rejected (cloud), 1 rejected (Tier 2)\n")
