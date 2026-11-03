# Energy billing analysis

Monthly energy bills of 92 mobile access sites (January 2019 to November 2021). The analysis flags suspicious bills and forecasts next-month consumption for each site.

## Requirements

- R 4.5
- [Quarto](https://quarto.org) (the PDF report is rendered with its bundled Typst engine)

## Setup

```
Rscript -e 'renv::restore()'
```

On Windows, install R 4.5 and Quarto with their standard installers, then run the same commands from the project root (PowerShell or the RStudio terminal). Rtools is only needed if a package has to be built from source.

## Run

From the project root:

```
Rscript run_all.R                  # writes tables, figures and results to outputs/
quarto render report/rapport.qmd   # renders the report to report/rapport.pdf
```

The data path defaults to `data/TimeSeries.csv` and can be overridden:

```
DATA_PATH=/path/to/file.csv Rscript run_all.R      # macOS, Linux
$env:DATA_PATH="C:\path\to\file.csv"; Rscript run_all.R   # Windows PowerShell
```

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

Only one person changes `renv.lock` at a time, in the same commit as the code that needs the new packages. After pulling a change to `renv.lock`, run `Rscript -e 'renv::restore()'`.

Line endings are normalised to LF by `.gitattributes`, so files edited on Windows do not show spurious changes.
