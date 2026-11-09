# Step 4: minimal end-to-end answer to both questions.
# Q1: deterministic checks and Hampel filter. Q2: per-day naive forecast.
# Run from the project root: Rscript scripts/02_baseline.R
source("R/config.R")
source("R/load.R")
source("R/features.R")
source("R/checks.R")
source("R/hampel.R")
source("R/naive.R")

for (d in c(config$tables_dir, config$results_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}
bills <- load_bills(config$data_path)

# ---- Q1: deterministic checks ----------------------------------------------
periods <- check_periods(bills)
coverage <- check_coverage(bills)
values <- check_values(bills)
duplicates <- check_duplicates(bills)
repeats <- check_repeats(bills, config$repeat_k_max)

# ---- Q1: Hampel filter ------------------------------------------------------
hampel_all <- hampel_bills(bills, config$hampel_half_width, config$hampel_threshold)
hampel_flags <- hampel_all[hampel_all$flag, ]
hampel_flags <- hampel_flags[order(-abs(hampel_flags$gap_kwh)), ]
hampel_flags$direction <- ifelse(hampel_flags$gap_kwh > 0, "over", "under")

# ---- Q2: per-day naive forecast of December 2021 ---------------------------
forecast <- naive_per_day_forecast(bills, config$forecast_month)
stopifnot(nrow(forecast) == length(unique(bills$site)), !anyNA(forecast$kwh))

# ---- Saved tables and numbers for the report -------------------------------
write.csv(periods, file.path(config$tables_dir, "check_periods.csv"), row.names = FALSE)
write.csv(coverage, file.path(config$tables_dir, "check_coverage.csv"), row.names = FALSE)
write.csv(repeats$repeats, file.path(config$tables_dir, "check_repeats.csv"), row.names = FALSE)
write.csv(hampel_flags, file.path(config$tables_dir, "hampel_flags.csv"), row.names = FALSE)
write.csv(forecast, file.path(config$tables_dir, "forecast_naive_dec2021.csv"), row.names = FALSE)

results <- list(
  n_rows          = nrow(bills),
  n_sites         = length(unique(bills$site)),
  periods         = periods,
  coverage        = coverage,
  n_values        = nrow(values),
  n_duplicates    = nrow(duplicates),
  repeats         = repeats,
  hampel_flags    = hampel_flags,
  n_hampel_sites  = length(unique(hampel_flags$site)),
  n_hampel_edge   = sum(hampel_flags$edge),
  hampel_half_width = config$hampel_half_width,
  hampel_threshold  = config$hampel_threshold,
  forecast        = forecast,
  forecast_total  = sum(forecast$kwh)
)
saveRDS(results, file.path(config$results_dir, "baseline.rds"))
cat("Baseline written to", config$output_dir, "\n")
