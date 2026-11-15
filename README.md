# Energy billing analysis

Monthly energy bills of 92 mobile access sites (January 2019 to November 2021). The analysis flags suspicious bills and forecasts next-month consumption for each site.

## Requirements

- R 4.5
- [Quarto](https://quarto.org) (the PDF report is rendered with its bundled Typst engine)

## Setup

```
Rscript -e "renv::restore()"
```

On Windows, install R 4.5 and Quarto with their standard installers, then run the same commands from the project root (PowerShell, cmd or the RStudio terminal). If `Rscript` is not found, add R's `bin` folder to the `PATH`, or run `renv::restore()` in the R console. Rtools is only needed if a package has to be built from source.

## Run

From the project root:

```
Rscript run_all.R                                  # tables, figures and results in outputs/
Rscript -e "testthat::test_dir('tests/testthat')"  # unit tests
quarto render report/rapport.qmd                   # report/rapport.pdf
```

On a MacBook Air M1 (8 GB), from a fresh clone: `renv::restore()` 9 s (packages already in the renv cache), `run_all.R` 50 min, tests 6 s, report 5 s. The calibration (`scripts/05_calibrate.R`) takes 46 of the 50 minutes: it uses up to 6 cores through forked workers, and a single core on Windows (same results, as seeds are set per task, but several hours). A fresh clone reproduced every table of `outputs/tables/` byte for byte.

The data path defaults to `data/TimeSeries.csv` and can be overridden:

```
DATA_PATH=/path/to/file.csv Rscript run_all.R      # macOS, Linux
$env:DATA_PATH="C:\path\to\file.csv"; Rscript run_all.R   # Windows PowerShell
```

## Pipeline

`run_all.R` runs the scripts in this order. Three of them reuse earlier results from `outputs/results/`: `04_suspects.R` (from 02 and 03), `08_select.R` (from 07) and `09_forecast.R` (from 08).

| Script | Step of `PLAN.md` | Main outputs |
|---|---|---|
| `01_descriptive.R` | 3, descriptive statistics | `tables/site_summary.csv` |
| `01b_course_view.R` | 3 bis, course tests on representative series | `tables/course_view_tests.csv` |
| `02_baseline.R` | 4, consistency checks, Hampel filter, naive forecast | `tables/check_periods.csv`, `tables/hampel_flags.csv` |
| `03_panel.R` | 5, expected value of each bill and level breaks | `tables/panel_breaks.csv` |
| `04_suspects.R` | 6, typology of the flagged bills | `tables/suspects.csv`, `tables/suspects_big.csv` |
| `05_calibrate.R` | 7, calibration by anomaly injection | `tables/calibration_grid.csv` |
| `06_baselines.R` | 8, rolling-origin backtest of the baselines | `tables/backtest_baselines.csv` |
| `07_models.R` | 9, course models site by site | `tables/models_accuracy.csv`, `tables/arima_diagnosis.csv` |
| `08_select.R` | 10, site classes and selection by class | `tables/selection_by_class.csv` |
| `09_forecast.R` | 11, test months and December 2021 forecast | `tables/forecast_dec2021.csv` |

Choice rules (thresholds, penalty, selection, intervals) are set in `R/config.R` and were written before each first run. The test months (June to November 2021) are used as forecast targets only in `09_forecast.R`; the forecasting scripts 06 to 08 stop before them.

## Layout

```
data/      raw bills (never modified)
R/         functions and configuration
scripts/   analysis steps, called by run_all.R
outputs/   generated tables, figures and results (not versioned)
report/    report source
```

## Workflow

Each step of `PLAN.md` is done on its own branch, created from an up-to-date `main`:

```
git switch main
git pull
git switch -c step-2-load
```

When the step's check passes (tests green, criterion of the plan met), push the branch and open a pull request to `main`. The other member reviews it before merging. The step's box in `PLAN.md` is ticked in that pull request. Nobody commits directly to `main`.

Only one person changes `renv.lock` at a time, in the same commit as the code that needs the new packages. After pulling a change to `renv.lock`, run `Rscript -e "renv::restore()"`.

Line endings are normalized to LF by `.gitattributes`, so files edited on Windows do not show spurious changes.
