# =============================================================================
# 00_install_packages.R - ONE-OFF PACKAGE INSTALLATION
# =============================================================================
# Run this ONCE, from a machine/node that HAS INTERNET, before the first run.
# HPC compute nodes frequently have no outbound network access, so installing
# at job time fails.
#
#   Rscript R/00_install_packages.R
#
# With a container:
#   singularity exec --bind $HOME:$HOME container.sif \
#     Rscript R/00_install_packages.R
#
# ### EDIT ### Set this to the same value as R_LIBS_USER in 00_config.R.
# NULL = install into R's default library (normal laptop use).
# =============================================================================

R_LIBS_USER <- NULL      # e.g. "~/R_libs" on a cluster

if (!is.null(R_LIBS_USER)) {
  R_LIBS_USER <- path.expand(R_LIBS_USER)
  if (!dir.exists(R_LIBS_USER)) {
    dir.create(R_LIBS_USER, recursive = TRUE)
    cat(sprintf("Created: %s\n", R_LIBS_USER))
  }
  .libPaths(c(R_LIBS_USER, .libPaths()))
}

target_lib <- if (is.null(R_LIBS_USER)) .libPaths()[1] else R_LIBS_USER
cat(sprintf("R version : %s\n", R.version.string))
cat(sprintf("Installing into: %s\n\n", target_lib))

# --- REQUIRED: the pipeline cannot run without these -------------------------
core <- c(
  "terra",      # raster processing
  "sf",         # vector / shapefile handling
  "dplyr",      # data manipulation
  "lubridate",  # dates
  "readr",      # CSV I/O
  "jsonlite",   # JSON MTL parsing
  "stringr",    # string handling
  "tibble",     # modern data frames
  "tidyr"       # reshaping
)

# --- OPTIONAL: figures (module 08). Tables are produced without them. --------
# NOTE: this pipeline deliberately does NOT depend on tidyterra or ggspatial.
# Those packages require recent R versions and are often unavailable on
# cluster installations; module 08 draws its north arrow, scale bar and
# rasters with base ggplot2 instead.
plotting <- c("ggplot2", "scales", "ggrepel")

# --- OPTIONAL: extras used by a few diagnostics ------------------------------
extras <- c("knitr", "pROC")

install_set <- function(pkgs, label) {
  have <- pkgs[sapply(pkgs, requireNamespace, quietly = TRUE)]
  need <- setdiff(pkgs, have)
  cat(sprintf("--- %s ---\n", label))
  if (length(have)) cat(sprintf("  already present: %s\n", paste(have, collapse = ", ")))
  if (!length(need)) { cat("  nothing to install\n\n"); return(invisible(NULL)) }
  cat(sprintf("  installing: %s\n", paste(need, collapse = ", ")))
  try(install.packages(need, lib = target_lib,
                       repos = "https://cloud.r-project.org",
                       dependencies = TRUE), silent = TRUE)
  cat("\n")
}

install_set(core,     "CORE (required)")
install_set(plotting, "PLOTTING (optional)")
install_set(extras,   "EXTRAS (optional)")

# --- verification -------------------------------------------------------------
cat("=== VERIFICATION ===\n")
report <- function(pkgs, label) {
  cat(sprintf("\n%s\n", label))
  bad <- character(0)
  for (p in pkgs) {
    ok <- requireNamespace(p, quietly = TRUE)
    cat(sprintf("  %-12s %s\n", p,
                if (ok) sprintf("OK   v%s", as.character(packageVersion(p))) else "MISSING"))
    if (!ok) bad <- c(bad, p)
  }
  bad
}
bad_core <- report(core,     "CORE:")
bad_plot <- report(plotting, "PLOTTING:")
report(extras, "EXTRAS:")

cat("\n")
if (length(bad_core) > 0) {
  cat(sprintf("FAILED (required): %s\n", paste(bad_core, collapse = ", ")))
  cat("The pipeline will NOT run. Missing system libraries are the usual cause.\n")
  cat("On Debian/Ubuntu:\n")
  cat("  sudo apt install libgdal-dev libproj-dev libgeos-dev libudunits2-dev\n")
} else if (length(bad_plot) > 0) {
  cat(sprintf("Figures unavailable (%s missing), but TABLES WILL STILL BE PRODUCED.\n",
              paste(bad_plot, collapse = ", ")))
  cat("Set GENERATE_PLOTS <- FALSE in R/00_config.R to silence the warnings.\n")
} else {
  cat("All good. Next: edit R/00_config.R, then run:\n  Rscript R/06_main.R\n")
}
