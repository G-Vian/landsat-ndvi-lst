# Landsat validation suite

Five diagnostic scripts that read the outputs of a Landsat NDVI/LST pipeline
and characterise the resulting series **before** they are used as covariates in
a statistical model.

They answer the questions a reviewer is likely to ask first. How much of the
record actually carries data, and where are the gaps? Does the satellite signal
track ground measurements? Do the spatial units differ from one another enough
to be worth modelling separately, and is that difference stable over time?

Nothing here is specific to a study area. Paths, labels, hemisphere and
excluded months are configured in one file.

---

## Contents

| File | Purpose |
|---|---|
| `00_common.R` | Shared configuration, input discovery, column normalisation, figure theme. **The only file you need to edit.** |
| `01_coverage_calendar.R` | Which months carry data and where the gaps are. **Run first** — the coverage structure conditions how everything else should be read. |
| `02_climate_lst_correlation.R` | Satellite LST against ground air temperature. Needs a monthly climate table. |
| `03_lst_by_zone_temporal.R` | How surface temperature varies between zones and within each zone over time. |
| `04_ndvi_by_zone_temporal.R` | The same for NDVI. |
| `05_lst_zone_ranking.R` | Whether the same zones are consistently the warmest. The key script for justifying a zone-level spatial exposure. |

Scripts 01–05 are independent and can be run in any order; only `00_common.R`
must sit alongside them.

---

## Requirements

**R ≥ 4.1.** Required packages: `readr`, `dplyr`, `tidyr`, `ggplot2`, `scales`.
Script 02 additionally needs `readxl` if your climate table is an Excel
workbook (not needed for CSV).

Optional but recommended:

```r
install.packages(c("ragg", "ggrepel"))
```

`ragg` renders figure text noticeably more cleanly on a headless Linux node;
`ggrepel` places non-overlapping labels in the stability map of script 05. Both
degrade gracefully if absent.

---

## Input expected

A completed run of an upstream Landsat pipeline that has written, somewhere
under a `tables/` or `tabelas/` folder:

- a **city-wide monthly table**, one row per year–month, with mean NDVI and LST;
- a **per-zone monthly table**, one row per zone–year–month, with the same.

**File names do not matter.** The suite locates each table by inspecting
column names rather than by a fixed path, and normalises whichever spelling it
finds onto one internal convention. It therefore works unchanged with a
pipeline that writes `ano`/`mes`/`nome_bairro` and with one that writes
`year`/`month`/`zone_name`.

Recognised spellings are listed in `COLUMN_ALIASES` in `00_common.R`. If your
pipeline uses a name that is not there, add it to that list — one line, and all
five scripts pick it up.

---

## Configuration

Everything is in **Section 1 of `00_common.R`**, marked `### EDIT ###`.

### Where the results are

```r
PIPELINE_DIR_CANDIDATES <- character(0)
```

Leave empty to search automatically: the suite looks for a directory containing
a `tables/` or `tabelas/` folder, beside these scripts and one or two levels up.
That covers the common case where the suite lives inside the project tree.

List explicit paths if the results are elsewhere. The first that exists is used,
so several machines or several runs can share one configuration:

```r
PIPELINE_DIR_CANDIDATES <- c(
  "~/my_project/results_landsat",
  "/scratch/user/project/Resultados_Landsat"
)
```

### Where the outputs go

```r
VALIDATION_ROOT <- NULL     # -> "validation_output" beside the scripts
```

### Study area label and period

```r
AOI_LABEL  <- "Study area"  # appears in figure titles and reports
YEAR_START <- NULL          # NULL = infer from the data
YEAR_END   <- NULL
```

Leaving the period `NULL` is usually right: the reports then describe the span
the pipeline actually produced, rather than a constant that may have gone stale.

### Hemisphere

```r
HEMISPHERE <- "south"       # "south" | "north" | NA
```

This controls the season labels. Getting it wrong silently mislabels every
seasonal figure, so it is worth a moment's check. Set `NA` near the equator,
where a four-season scheme does not describe the climate and the label would be
actively misleading — seasonal grouping is then disabled rather than faked.

### Excluded months

```r
LOW_QUALITY_MONTHS <- data.frame(
  year = integer(0), month = integer(0), reason = character(0)
)
```

Some months survive every automatic filter in a pipeline yet are clearly
unusable on inspection. The classic case is a diffuse cold bias from Landsat-7
SLC-off striping: no single pixel is extreme enough to be clamped, so nothing
trips, while the scene mean is badly wrong.

Start empty. Run script 01, inspect the city-wide series, and add any month
whose mean is physically implausible **for your climate**, with a reason you
would be willing to defend in print:

```r
LOW_QUALITY_MONTHS <- data.frame(
  year   = 2021,
  month  = 7,
  reason = "LE07 SLC-off, diffuse cold bias (city mean LST implausible)"
)
```

Declared once here, the exclusion propagates to all five scripts: those rows
are dropped from every statistic, and the exclusion is stated in figure
captions and text reports automatically.

### Script 02 only

The climate table path and its temperature column names are set at the top of
`02_climate_lst_correlation.R`. Both can be left `NULL`, in which case the
script searches for a climate-looking spreadsheet or CSV nearby and matches
columns by keyword — printing exactly what it chose. That is a convenience for
a first run; **for anything you intend to publish, write the names out
explicitly** so the choice is on the record.

If a column name does not match, the script stops and prints every column it
did find, so the fix is to copy the right one from that list.

---

## Running

```bash
Rscript 01_coverage_calendar.R
Rscript 02_climate_lst_correlation.R
Rscript 03_lst_by_zone_temporal.R
Rscript 04_ndvi_by_zone_temporal.R
Rscript 05_lst_zone_ranking.R
```

Inside a container:

```bash
singularity exec --bind $HOME:$HOME \
  --env R_LIBS="$HOME/R_libs:/usr/local/lib/R/site-library:/usr/local/lib/R/library" \
  /path/to/container.sif \
  Rscript 01_coverage_calendar.R
```

Each takes well under a minute: they read tables, not rasters.

---

## Outputs

Each script writes into its own subfolder of `VALIDATION_ROOT`:

```
01_coverage_calendar/   3 figures, 3 CSV tables, 1 text report
02_climate_lst/         6 figures, 2 CSV tables, 1 text report
03_lst_by_zone/         6 figures, 3 CSV tables, 1 text report
04_ndvi_by_zone/        6 figures, 3 CSV tables, 1 text report
05_lst_ranking/         6 figures, 3 CSV tables, 1 text report
```

Every figure is written twice: a 320 dpi PNG and a vector PDF. Submit the PDF
where the journal accepts vector artwork, since its text and lines stay sharp at
any magnification. Adjust `FIG_DPI`, or set `FIG_SAVE_PDF <- FALSE` in
`00_common.R`, for fewer or lighter files.

The text reports are written for a reader who has not seen the code: each states
what was computed, what the numbers mean, and where the caveats are.

---

## Reading the results

**Start with 01.** If half the months are empty, every subsequent statistic is
conditional on a thin and possibly non-random sample of dates, and should be
described that way.

**Check the pixel counts.** Zones with few valid pixels produce unstable means.
The per-zone tables carry the counts; a mean over a handful of pixels is not
comparable to one over tens of thousands.

**Script 05 is the one that justifies a spatial exposure.** If zone rankings
are unstable from year to year, a time-invariant zone-level covariate is
describing noise rather than a persistent spatial pattern, and that finding
matters more than any single correlation.

---

## Troubleshooting

| Symptom | Cause | Fix |
|---|---|---|
| `Could not find pipeline results automatically` | Results elsewhere | Set `PIPELINE_DIR_CANDIDATES` |
| `Could not find a city-wide monthly CSV` | Column names not recognised | Add the spelling to `COLUMN_ALIASES` |
| Seasons look wrong | Hemisphere | Set `HEMISPHERE` |
| `Fewer than N usable months` | Climate and satellite periods do not overlap | Check both cover the same years |
| Figure text looks rough | `ragg` absent | `install.packages("ragg")` |
| Month names not in English | Locale | The suite forces `LC_TIME=C`; report it if it recurs |

---

## Licence

MIT. Landsat data are courtesy of the U.S. Geological Survey and are in the
public domain.
