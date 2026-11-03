# Full analysis pipeline. Run from the project root: Rscript run_all.R
if (!file.exists("R/config.R")) stop("run_all.R must be run from the project root")
source("R/config.R")

bills <- read.csv(config$data_path, colClasses = "character")
cat("Rows read:", nrow(bills), "\n")

dir.create(config$results_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(list(n_rows = nrow(bills)), file.path(config$results_dir, "skeleton.rds"))
