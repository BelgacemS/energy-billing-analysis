# Full analysis pipeline. Run from the project root: Rscript run_all.R
if (!file.exists("R/config.R")) stop("run_all.R must be run from the project root")

source("scripts/01_descriptive.R", local = new.env())
source("scripts/02_baseline.R", local = new.env())
