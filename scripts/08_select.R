# Step 10: classes of sites, extra candidates (month-adjusted naive, robust
# local level, combination) and selection of the method by class on the
# validation months. Needs the step 9 backtest (models.rds). About 1 minute.
# Run from the project root: Rscript scripts/08_select.R
source("R/config.R")
source("R/load.R")
source("R/features.R")
source("R/naive.R")
source("R/backtest.R")
source("R/segment.R")
source("R/panel.R")
source("R/robust.R")
source("R/classes.R")

for (d in c(config$figures_dir, config$tables_dir, config$results_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

bills <- load_bills(config$data_path)
pre <- bills[bills$month < config$test_start, ]  # test months are never read here
big <- big_sites(pre, config$big_site_kwh)
seg_args <- list(cost = config$segment_cost, penalty = config$segment_penalty,
                 plateau_min = config$plateau_min_months, K = config$segment_K)
mo <- readRDS(file.path(config$results_dir, "models.rds"))

# ---- Classes on the months before the validation -----------------------------
feat <- class_features(bills, big, config$validation_start, seg_args, config$level_shift_min)
feat_cols <- c("level", "slope", "noise", "breaks")
cl <- ward_classes(feat, feat_cols, config$class_k_range)
feat$class <- paste0("C", cl$class)
profile <- do.call(rbind, lapply(split(feat, feat$class), function(f) {
  z <- scale(as.matrix(feat[feat_cols]))[feat$class == f$class[1], , drop = FALSE]
  centre <- colMeans(z)
  near <- f$site[order(rowSums(sweep(z, 2, centre)^2))][1:min(3, nrow(f))]
  data.frame(class = f$class[1], n = nrow(f), level_kwh_day = exp(stats::median(f$level)),
             slope = exp(stats::median(f$slope)) - 1, noise = exp(stats::median(f$noise)),
             breaks = stats::median(f$breaks), share_break = mean(f$breaks > 0),
             examples = paste(near, collapse = ", "), stringsAsFactors = FALSE)
}))
set.seed(config$calibration_seed)
suffix_tab <- table(class = feat$class, suffix = feat$suffix)
suffix_p <- stats::fisher.test(suffix_tab, simulate.p.value = TRUE, B = 1e4)$p.value
site_class <- c(stats::setNames(feat$class, feat$site), stats::setNames(rep("big", length(big)), big))

# ---- Extra candidates on the validation months --------------------------------
targets <- seq(config$validation_start, add_months(config$test_start, -1), by = "month")
methods <- list(naive_per_day = naive_per_day_forecast, naive_adj = naive_adj_forecast,
                robust_level = robust_level_forecast)
methods$combination <- combination_forecaster(methods[config$combination_methods])
timing <- system.time(bt_new <- backtest(
  bills, targets, methods, config$test_start,
  robust_preprocess(big, seg_args, config$level_shift_min)))[["elapsed"]]

# Same per-day naive as step 9 (sanity check), then the step 9 models.
bt9 <- mo$backtest
key <- function(d) paste(d$site, d$month)
n9 <- bt9[bt9$method == "naive_per_day", ]
nn <- bt_new[bt_new$method == "naive_per_day", ]
stopifnot(isTRUE(all.equal(n9$forecast[order(key(n9))], nn$forecast[order(key(nn))])))
cols <- c("method", "month", "site", "billed", "forecast", "error", "rel_error", "scale",
          "scaled_error")
bt <- rbind(bt_new[cols], bt9[bt9$method != "naive_per_day", cols])
bt$class <- unname(site_class[bt$site])

# ---- Selection by class -------------------------------------------------------
mase <- tapply(bt$scaled_error, list(bt$class, bt$method), mean)  # (classes, methods)
ord_mase <- tapply(bt$scaled_error[bt$class != "big"], bt$method[bt$class != "big"], mean)
selection <- do.call(rbind, lapply(rownames(mase), function(k) {
  m <- mase[k, ]
  data.frame(class = k, best = names(which.min(m)), best_mase = min(m),
             naive_mase = m[["naive_per_day"]],
             chosen = select_method(m, config$method_order, config$selection_tolerance),
             stringsAsFactors = FALSE)
}))
selection$chosen_mase <- mase[cbind(selection$class, selection$chosen)]
pooled_choice <- select_method(ord_mase, config$method_order, config$selection_tolerance)
accuracy <- rbind(cbind(scope = "ordinary", backtest_accuracy(bt[bt$class != "big", ])),
                  cbind(scope = "big", backtest_accuracy(bt[bt$class == "big", ])))

labels <- c(naive_per_day = "Naïve par jour", naive_adj = "Naïve corrigée du mois",
            robust_level = "Niveau robuste", combination = "Combinaison",
            ses = "Lissage simple", ets = "ETS", arima_noseason = "auto.arima sans saison",
            arima = "auto.arima par défaut")
# Conclusion: 5 classes (mean silhouette 0.36, weak structure), not linked to
# the suffix (Fisher p = 0.66). On the 90 ordinary sites the per-day naive
# stays best (MASE 0.95), then the month-adjusted naive and the combination
# (1.01); the robust level is worst of the new candidates (1.22): it keeps the
# old level for the first months of a new one (1088KO, August 2020). The rule
# keeps the per-day naive for 4 of the 5 classes and the combination for C4
# (4 sites, 0.245 against 0.251, 2.4 % lower, just past the 2 % tolerance);
# auto.arima without seasonal difference for the big sites (0.79 against 0.95).

# ---- Saved tables and numbers for the report ---------------------------------
write.csv(feat, file.path(config$tables_dir, "site_classes.csv"), row.names = FALSE)
write.csv(selection, file.path(config$tables_dir, "selection_by_class.csv"), row.names = FALSE)

results <- list(
  config = config[c("class_k_range", "combination_methods", "method_order",
                    "selection_tolerance")],
  timing_seconds = timing,
  features = feat,
  feature_cols = feat_cols,
  k = cl$k,
  silhouette = cl$silhouette,
  profile = profile,
  suffix_tab = suffix_tab,
  suffix_p = suffix_p,
  site_class = site_class,
  backtest = bt,
  mase = mase,
  ordinary_mase = ord_mase,
  accuracy = accuracy,
  selection = selection,
  pooled_choice = pooled_choice,
  labels = labels
)
saveRDS(results, file.path(config$results_dir, "select.rds"))
cat("Selection written to", config$output_dir, "in", round(timing), "s of backtest\n")
