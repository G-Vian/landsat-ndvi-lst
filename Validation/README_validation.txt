================================================================================
  LANDSAT VALIDATION SUITE - NDVI & LST - Santos, Brazil
================================================================================

  Five diagnostic scripts that read the outputs of the main Landsat pipeline
  and characterise the resulting NDVI and LST series before they are used as
  covariates in a statistical model.

--------------------------------------------------------------------------------
  FILES
--------------------------------------------------------------------------------

  00_common.R                  Shared paths, input discovery, column-name
                               normalisation, figure theme and saving.
                               ### THE ONLY FILE YOU NEED TO EDIT ###

  01_coverage_calendar.R       Which months of the series carry data, and where
                               the gaps are. Run this first: the coverage
                               structure conditions how every other result
                               should be read.

  02_climate_lst_correlation.R Satellite LST against ground air temperature.
                               Requires the monthly climate spreadsheet.

  03_lst_by_zone_temporal.R    How surface temperature varies between zones and
                               within each zone over time.

  04_ndvi_by_zone_temporal.R   The same for NDVI.

  05_lst_zone_ranking.R        Whether the same zones are consistently the
                               warmest. The key script for justifying a
                               zone-level spatial exposure.

--------------------------------------------------------------------------------
  BEFORE THE FIRST RUN
--------------------------------------------------------------------------------

  Open 00_common.R and check two things in Section 1:

    VALIDATION_ROOT          where the outputs go
    PIPELINE_DIR_CANDIDATES  where the main pipeline wrote its results

  The scripts locate the input CSVs by inspecting their COLUMN NAMES rather
  than by a fixed file name, so they work with both the Portuguese pipeline
  (tabelas/, columns ano/mes/nome_bairro) and the English one (tables/,
  columns year/month/zone_name) without any change.

  Script 02 additionally needs PATH_CLIMATE and the three COL_TEMP_* column
  names set at the top of that file. If a column name is wrong the script stops
  and prints every column it did find, so the fix is to copy the right name
  from that list.

  Optional but recommended: install the ragg and ggrepel packages. ragg renders
  figure text noticeably more cleanly on a headless Linux node; ggrepel places
  non-overlapping labels in the stability map of script 05. Both degrade
  gracefully if absent.

    Rscript -e 'install.packages(c("ragg","ggrepel"), lib="/home/g.vian/R_libs", repos="https://cloud.r-project.org")'

--------------------------------------------------------------------------------
  RUNNING
--------------------------------------------------------------------------------

  module add singularity
  cd /home/g.vian/Pesquisa_Epidemic/PROJETO_SANTOS/Dados_Satelite/files/Validation

  singularity exec --bind /home/g.vian:/home/g.vian \
    --env R_LIBS="/home/g.vian/R_libs:/usr/local/lib/R/site-library:/usr/local/lib/R/library" \
    /home/public/R_inla/r_inla.sif \
    Rscript 01_coverage_calendar.R

  Repeat for 02 through 05. They are independent of one another and can be run
  in any order; only 00_common.R must sit alongside them. Each takes well under
  a minute, since they read tables rather than rasters.

--------------------------------------------------------------------------------
  OUTPUTS
--------------------------------------------------------------------------------

  Each script writes into its own subfolder of VALIDATION_ROOT:

    01_coverage_calendar/   3 figures, 3 CSV tables, 1 text report
    02_climate_lst/         6 figures, 2 CSV tables, 1 text report
    03_lst_by_zone/         6 figures, 3 CSV tables, 1 text report
    04_ndvi_by_zone/        6 figures, 3 CSV tables, 1 text report
    05_lst_ranking/         6 figures, 3 CSV tables, 1 text report

  Every figure is written twice: a 320 dpi PNG and a vector PDF. Submit the PDF
  where the journal accepts vector artwork, since its text and lines stay sharp
  at any magnification. Adjust FIG_DPI or set FIG_SAVE_PDF to FALSE in
  00_common.R if you want fewer or lighter files.

  The text reports are written for a reader who has not seen the code: each one
  states what was computed, what the numbers mean, and where the caveats are.

--------------------------------------------------------------------------------
  A NOTE ON THE EXCLUDED MONTH
--------------------------------------------------------------------------------

  July 2021 is excluded from every statistic. The Landsat-7 scene of that month
  reports a city-wide mean surface temperature of 15.19 C, which is not
  plausible for a mid-winter daytime overpass in Santos, where 22-30 C is
  expected. The cause is diffuse cold bias from SLC-off striping. It passes
  every automatic filter in the pipeline because no pixel is clamped, so it has
  to be removed by hand.

  The exclusion is declared in one place, LOW_QUALITY_MONTHS in 00_common.R.
  Add rows there if further months are found to be unreliable; all five scripts
  pick the change up automatically.

================================================================================
