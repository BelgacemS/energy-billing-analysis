# Step 8: rolling-origin backtest and baseline forecasts, validation only.
# Run from the project root: Rscript scripts/06_baselines.R
source("R/config.R")
source("R/load.R")
source("R/features.R")
source("R/naive.R")
source("R/backtest.R")
suppressPackageStartupMessages(library(ggplot2))

for (d in c(config$figures_dir, config$tables_dir, config$results_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}
save_fig <- function(plot, name, width = 7, height = 4.5) {
  ggsave(file.path(config$figures_dir, paste0(name, ".png")), plot,
         width = width, height = height, dpi = 300, bg = "white")
}
theme_set(theme_minimal(base_size = 10))

bills <- load_bills(config$data_path)
big <- big_sites(bills, config$big_site_kwh)

# ---- Backtest on the validation months --------------------------------------
targets <- seq(config$validation_start, add_months(config$test_start, -1), by = "month")
baselines <- list(naive_per_day = naive_per_day_forecast, naive_raw = naive_raw_forecast,
                  snaive_per_day = snaive_per_day_forecast)
bt <- backtest(bills, targets, baselines, config$test_start)
stopifnot(length(targets) == 11, max(bt$month) < config$test_start,
          nrow(bt) == length(baselines) * length(targets) * length(unique(bills$site)))
bt$big <- bt$site %in% big

accuracy <- rbind(cbind(scope = "all", backtest_accuracy(bt)),
                  cbind(scope = "ordinary", backtest_accuracy(bt[!bt$big, ])),
                  cbind(scope = "big", backtest_accuracy(bt[bt$big, ])))
rownames(accuracy) <- NULL
by_month <- do.call(rbind, lapply(split(bt[!bt$big, ], bt$month[!bt$big]), backtest_accuracy))
by_month$month <- rep(targets, each = length(baselines))
rownames(by_month) <- NULL

# Per-day naive against raw naive, bill by bill (ordinary sites).
wide <- reshape(bt[!bt$big, c("method", "month", "site", "rel_error")], idvar = c("month", "site"),
                timevar = "method", direction = "wide")
# Equal when the origin and target months have the same number of days.
per_day_vs_raw <- c(better = mean(wide$rel_error.naive_per_day < wide$rel_error.naive_raw),
                    equal = mean(wide$rel_error.naive_per_day == wide$rel_error.naive_raw))

# December 2020, the only validation month like the forecast target.
december <- by_month[by_month$month == as.Date("2020-12-01"), ]

# Selection rule applied to the baselines alone, ordinary sites: the method
# the models of steps 9 and 10 must beat.
ord_acc <- accuracy[accuracy$scope == "ordinary", ]
best_baseline <- select_method(setNames(ord_acc$mase, ord_acc$method), config$method_order,
                               config$selection_tolerance)

# Sites where the seasonal naive beats the per-day naive (mean MASE).
site_mase <- tapply(bt$scaled_error, list(bt$site, bt$method), mean)  # (sites, methods)
snaive_wins <- sum(site_mase[!rownames(site_mase) %in% big, "snaive_per_day"] <
                     site_mase[!rownames(site_mase) %in% big, "naive_per_day"])

# ---- Figure ------------------------------------------------------------------
labels <- c(naive_per_day = "Naïve par jour", naive_raw = "Naïve brute",
            snaive_per_day = "Saisonnière naïve par jour")
by_month$label <- factor(labels[by_month$method], levels = labels)
p <- ggplot(by_month, aes(month, mase, colour = label)) +
  geom_line() +
  geom_point(size = 1.2) +
  scale_colour_manual(values = c("#08519c", "#6baed6", "#e6550d"), name = NULL) +
  scale_x_date(date_labels = "%m/%Y", breaks = targets[c(TRUE, FALSE)]) +
  labs(x = "Mois prévu", y = "MASE moyen (sites ordinaires)") +
  theme(legend.position = "bottom")
save_fig(p, "backtest_baselines")
# Conclusion: on the 90 ordinary sites and the 11 validation months, the
# per-day naive has a mean MASE of 0.95 and a median relative error of
# 1.1 % (90th percentile 4.4 %). The mean relative error, 4.3 %, is half
# made of one bill: the 1240BR drop of May 2021, billed 138 kWh for 3,128
# forecast (2,167 %): an error relative to the billed value explodes on a
# small bill, hence the MASE for the choices. The per-day naive is better
# than the raw naive (2.37, 3.8 %) on 71 % of the bills and equal on 18 %
# (months of the same length), and better than the seasonal naive (4.43,
# 6.4 %) in every month, December 2020 included (0.64 against 1.42 and
# 4.29). A level a year old misses the
# changes of level (step 3), which weigh far more than the weak month
# effect: the seasonal naive is better on 3 ordinary sites only. The two big
# sites carry 248,412 of the 291,430 kWh of absolute error of the per-day
# naive and are under-forecast overall (+51,260 kWh), as expected from their
# growth. The rule keeps the per-day naive as the baseline to beat.

# ---- Saved tables and numbers for the report ---------------------------------
write.csv(accuracy, file.path(config$tables_dir, "backtest_baselines.csv"), row.names = FALSE)

results <- list(
  config         = config[c("validation_start", "test_start", "selection_tolerance", "method_order")],
  targets        = targets,
  test_targets   = seq(config$test_start, max(bills$month), by = "month"),  # not read here
  n_sites        = length(unique(bills$site)),
  big            = big,
  labels         = labels,
  backtest       = bt,
  accuracy       = accuracy,
  by_month       = by_month,
  december       = december,
  per_day_vs_raw = per_day_vs_raw,
  best_baseline  = best_baseline,
  snaive_wins    = snaive_wins
)
saveRDS(results, file.path(config$results_dir, "backtest_baselines.rds"))
cat("Baselines written to", config$output_dir, "\n")
