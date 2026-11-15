# Step 9: course models site by site (SES, ETS, auto.arima) on the validation
# months, and diagnosis of auto.arima. About 4 minutes.
# Run from the project root: Rscript scripts/07_models.R
source("R/config.R")
source("R/load.R")
source("R/features.R")
source("R/naive.R")
source("R/backtest.R")
source("R/models.R")
source("R/segment.R")
source("R/panel.R")
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
pre <- bills[bills$month < config$test_start, ]  # test months are never read here
big <- big_sites(pre, config$big_site_kwh)

# ---- Backtest on the validation months --------------------------------------
targets <- seq(config$validation_start, add_months(config$test_start, -1), by = "month")
methods <- list(
  naive_per_day  = naive_per_day_forecast,
  ses            = model_forecaster(fit_ses),
  ets            = model_forecaster(make_fit_ets(config$ets_additive_only)),
  arima          = model_forecaster(make_fit_arima(NA)),
  arima_noseason = model_forecaster(make_fit_arima(0, config$arima_season_sites))
)
timing <- system.time(bt <- backtest(bills, targets, methods, config$test_start))[["elapsed"]]
bt$big <- bt$site %in% big
bt$fallback <- grepl("^fallback", bt$info)

accuracy <- rbind(cbind(scope = "all", backtest_accuracy(bt)),
                  cbind(scope = "ordinary", backtest_accuracy(bt[!bt$big, ])),
                  cbind(scope = "big", backtest_accuracy(bt[bt$big, ])))
rownames(accuracy) <- NULL
by_month <- do.call(rbind, lapply(split(bt[!bt$big, ], bt$month[!bt$big]), backtest_accuracy))
by_month$month <- rep(targets, each = length(methods))
rownames(by_month) <- NULL

# Per site: mean validation MASE of each model against the per-day naive.
site_mase <- tapply(bt$scaled_error, list(bt$site, bt$method), mean)  # (sites, methods)
ord_sites <- !rownames(site_mase) %in% big
beats_naive <- colSums(site_mase[ord_sites, ] < site_mase[ord_sites, "naive_per_day"])

# Orders chosen by auto.arima at each origin (ordinary sites).
ar <- bt[bt$method %in% c("arima", "arima_noseason") & !bt$big & !bt$fallback, ]
ar$kind <- arima_kind(ar$info)
ar$degenerate <- grepl("^ARIMA\\(0,\\d,0\\)\\(0,\\d,0\\)", ar$info)
ar$naive_like <- ar$kind == "naive"
ar$seasonal_diff <- ar$kind == "seasonal_diff"
orders_share <- aggregate(cbind(degenerate, naive_like, seasonal_diff) ~ method, ar, mean)
# auto.arima sets D = 0 without testing when the history has 2 x 12 months
# or less (forecast 9.0.2), so a seasonal difference can only appear from the
# origin with 25 months: sites with D = 1 per target month, default variant.
seasonal_by_month <- tapply(ar$seasonal_diff[ar$method == "arima"], ar$month[ar$method == "arima"], sum)
top_orders <- lapply(split(ar$info, ar$method), function(x) head(sort(table(x), decreasing = TRUE), 5))

# MASE of each kind of auto.arima fit against the per-day naive on the same
# site-months: a random walk without drift is the naive itself, so the gap
# between the two comes from the other kinds.
nv <- bt[bt$method == "naive_per_day", c("month", "site", "scaled_error")]
ar <- merge(ar, nv, by = c("month", "site"), suffixes = c("", "_naive"))
kind_mase <- aggregate(cbind(scaled_error, scaled_error_naive) ~ method + kind, ar, mean)
kind_mase$n <- aggregate(scaled_error ~ method + kind, ar, length)$scaled_error

# Selection rule on all ordinary sites pooled, for information: the rule is
# applied by class at step 10.
ord_acc <- accuracy[accuracy$scope == "ordinary", ]
pooled_choice <- select_method(setNames(ord_acc$mase, ord_acc$method), config$method_order,
                               config$selection_tolerance)

# ---- Diagnosis of auto.arima on the months before the test ------------------
seg_args <- list(cost = config$segment_cost, penalty = config$segment_penalty,
                 plateau_min = config$plateau_min_months, K = config$segment_K)
panel_pre <- panel_fit(pre, month_effect(pre, big), big, seg_args)
large <- panel_pre$breaks[abs(panel_pre$breaks$change) > config$level_shift_min, ]
diagnosis <- do.call(rbind, lapply(split(pre, pre$site), function(d) {
  d <- d[order(d$month), ]
  y <- site_ts(d)
  reg <- regimes(nrow(d), match(large$month[large$site == d$site[1]], d$month))
  do.call(rbind, lapply(c("arima", "arima_noseason"), function(v) {
    D <- if (v == "arima" || d$site[1] %in% config$arima_season_sites) NA else 0
    fit <- suppressWarnings(forecast::auto.arima(y, D = D))
    cbind(site = d$site[1], variant = v, n_regimes = max(reg),
          arima_diagnosis(fit, y, reg, config$ljung_box_alpha))
  }))
}))
fails_bt <- unique(bt[bt$fallback & bt$method %in% c("arima", "arima_noseason"), c("site", "method")])
diagnosis$fit_failed <- paste(diagnosis$site, diagnosis$variant) %in%
  paste(fails_bt$site, fails_bt$method)
diagnosis$val_mase <- site_mase[cbind(diagnosis$site, diagnosis$variant)]
diagnosis$naive_mase <- site_mase[diagnosis$site, "naive_per_day"]
diagnosis$big <- diagnosis$site %in% big
reasons <- with(diagnosis, cbind(degenerate, break_diff, lb_reject, fit_failed))
diagnosis$reasons <- apply(reasons, 1, function(r) paste(colnames(reasons)[r], collapse = ", "))
diagnosis$flagged <- rowSums(reasons) > 0
rownames(diagnosis) <- NULL

# Counts per variant and reason (ordinary sites), and accuracy of flagged
# sites against the others.
diag_ord <- diagnosis[!diagnosis$big, ]
reason_counts <- aggregate(cbind(degenerate, break_diff, lb_reject, fit_failed, flagged) ~ variant,
                           diag_ord, sum)
seasonal_diff_sites <- tapply(diag_ord$D == 1, diag_ord$variant, sum)
flag_vs_mase <- aggregate(cbind(val_mase, naive_mase) ~ variant + flagged, diag_ord, mean)
degenerate_models <- table(diag_ord$model[diag_ord$degenerate & diag_ord$variant == "arima_noseason"])
n_pre_months <- length(unique(pre$month))
lb_lag <- unique(diag_ord$lb_lag)

# The break criterion on random walks (D = 0 variant), two nulls per
# ordinary site (see R/config.R): expected number of sites flagged although
# differencing is the right model.
set.seed(config$break_null_seed)
null_rate <- do.call(rbind, lapply(split(pre[!pre$site %in% big, ], pre$site[!pre$site %in% big]),
                                   function(d) {
  inc <- diff(log(d$kwh_per_day[order(d$month)]))
  run <- function(draw) break_diff_null_rate(draw, nrow(d), config$break_null_reps, seg_args,
                                             config$level_shift_min, config$ljung_box_alpha)
  c(gaussian_mad = run(function(n) stats::rnorm(n, sd = stats::mad(inc))),
    bootstrap = run(function(n) sample(inc, n, replace = TRUE)))
}))  # (sites, nulls)
break_null_expected <- colSums(null_rate)

# ---- Figure ------------------------------------------------------------------
labels <- c(naive_per_day = "Naïve par jour", ses = "Lissage exponentiel simple", ets = "ETS",
            arima = "auto.arima, réglages par défaut",
            arima_noseason = "auto.arima, sans différence saisonnière")
by_month$label <- factor(labels[by_month$method], levels = labels)
p <- ggplot(by_month, aes(month, mase, colour = label)) +
  geom_line() +
  geom_point(size = 1.2) +
  scale_colour_manual(values = c("#08519c", "#74c476", "#31a354", "#fd8d3c", "#e6550d"),
                      name = NULL) +
  scale_x_date(date_labels = "%m/%Y", breaks = targets[c(TRUE, FALSE)]) +
  labs(x = "Mois prévu", y = "MASE moyenne (sites ordinaires)") +
  guides(colour = guide_legend(nrow = 2)) +
  theme(legend.position = "bottom")
save_fig(p, "backtest_models")
# Conclusion: on the 90 ordinary sites, no course model beats the per-day
# naive (mean validation MASE 0.95): SES 1.02, ETS 1.04, auto.arima without
# seasonal difference 1.16, default auto.arima 1.27; the naive is best in
# 10 of 11 months (SES 1.420 against 1.437 in May 2021), and the pooled
# rule keeps it. The models mostly rebuild the naive: SES alpha >= 0.95 in
# 75 % of the fits, ETS(A,N,N) (that is SES) in 91 %, ARIMA(0,1,0)
# without drift (the per-day naive itself) in 37 % (default) and 42 %
# (D = 0) of the auto.arima fits. Every other kind
# of fit does worse than the naive on the same site-months (D = 0 variant,
# MASE of the model against the naive): reversion to a mean 1.78 / 1.31,
# drift 1.06 / 0.76, other AR or MA terms 0.63 / 0.47; in the default
# variant, a seasonal difference 1.49 / 0.68: auto.arima does not consider
# it with 24 months of history or less, then chooses it for 36 to 40 of the
# 90 sites each month from the February 2021 forecast. On the two big sites
# the reverse holds: auto.arima D = 0 (which keeps the seasonal difference
# for 0698BR) 0.79 and ETS 0.84 against 0.95.
# Diagnosis on January 2019 to May 2021 (D = 0 variant): 54 degenerate
# orders (44 random walks without drift, i.e. the naive, 6 with drift, 3
# constant means, 1 with d = 2), 42 sites whose differencing disappears
# once the breaks are removed, 5 with autocorrelated residuals (Ljung-Box on
# 5 lags has little power on 29 months), no failed fit. The break
# criterion cannot tell a break from a random walk: on random walks it
# would flag 8.2 sites with Gaussian increments of MAD scale and 35.4 with
# each site's own increments resampled, against 42 observed. 71 of the 90
# sites are flagged; auto.arima errs 1.25 times as much as the naive on
# them, 1.11 times on the others.

# ---- Saved tables and numbers for the report ---------------------------------
write.csv(accuracy, file.path(config$tables_dir, "models_accuracy.csv"), row.names = FALSE)
write.csv(diagnosis, file.path(config$tables_dir, "arima_diagnosis.csv"), row.names = FALSE)

results <- list(
  config          = config[c("ets_additive_only", "arima_season_sites", "ljung_box_alpha")],
  timing_seconds  = timing,
  targets         = targets,
  big             = big,
  labels          = labels,
  backtest        = bt,
  accuracy        = accuracy,
  by_month        = by_month,
  n_fallback      = tapply(bt$fallback, bt$method, sum),
  beats_naive     = beats_naive,
  n_ordinary      = sum(ord_sites),
  orders_share    = orders_share,
  seasonal_by_month = seasonal_by_month,
  top_orders      = top_orders,
  kind_mase       = kind_mase,
  pooled_choice   = pooled_choice,
  n_large_breaks  = nrow(large),
  diagnosis       = diagnosis,
  reason_counts   = reason_counts,
  seasonal_diff_sites = seasonal_diff_sites,
  flag_vs_mase    = flag_vs_mase,
  degenerate_models = degenerate_models,
  n_pre_months    = n_pre_months,
  lb_lag          = lb_lag,
  break_null_rate = null_rate,
  break_null_expected = break_null_expected
)
saveRDS(results, file.path(config$results_dir, "models.rds"))
cat("Models written to", config$output_dir, "in", round(timing), "s of backtest\n")
