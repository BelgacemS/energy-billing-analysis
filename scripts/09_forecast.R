# Step 11: single reading of the test months with the strategy of step 10,
# reliability flags, prediction intervals and forecast of December 2021.
# Needs select.rds (step 10). About 10 seconds.
# Run from the project root: Rscript scripts/09_forecast.R
source("R/config.R")
source("R/load.R")
source("R/features.R")
source("R/naive.R")
source("R/backtest.R")
source("R/models.R")
source("R/segment.R")
source("R/panel.R")
source("R/robust.R")
source("R/intervals.R")

for (d in c(config$tables_dir, config$results_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}
bills <- load_bills(config$data_path)
pre <- bills[bills$month < config$test_start, ]
big <- big_sites(pre, config$big_site_kwh)
sel <- readRDS(file.path(config$results_dir, "select.rds"))
seg_args <- list(cost = config$segment_cost, penalty = config$segment_penalty,
                 plateau_min = config$plateau_min_months, K = config$segment_K)
prep_fun <- robust_preprocess(big, seg_args, config$level_shift_min)
levels <- config$interval_levels

# ---- Strategy and static reliability criteria ---------------------------------
site_class <- sel$site_class
site_method <- setNames(sel$selection$chosen[match(site_class, sel$selection$class)],
                        names(site_class))
base <- list(naive_per_day = naive_per_day_forecast, naive_adj = naive_adj_forecast,
             robust_level = robust_level_forecast)
forecasters <- list(
  naive_per_day = naive_per_day_forecast,
  combination = combination_forecaster(base[config$combination_methods]),
  arima_noseason = model_forecaster(make_fit_arima(0, config$arima_season_sites))
)
stopifnot(all(site_method %in% names(forecasters)))
strategy <- strategy_forecaster(site_method, forecasters)

vb <- sel$backtest
vb <- vb[vb$method == site_method[vb$site], ]  # validation rows of the strategy
poor_fit <- tapply(vb$scaled_error, vb$site, mean) > config$reliability_mase_max
noise <- vapply(split(pre, pre$site), function(d) noise_scale(log(d$kwh_per_day[order(d$month)])),
                numeric(1))
noisy <- noise > stats::quantile(noise[!names(noise) %in% big], config$reliability_noise_quantile)

# End-of-series anomaly at the origin of each target month, and the sigma
# of the Gaussian interval (per-day naive), from the bills up to the origin.
origin_state <- function(targets) {
  do.call(rbind, lapply(targets, function(t) {
    train <- bills[bills$month < t, ]
    ea <- end_anomaly(prep_fun(train), config$edge_months)
    sd <- log_change_sd(train)
    data.frame(month = t, site = names(ea), end_anomaly = unname(ea),
               sigma = unname(sd[names(ea)]), stringsAsFactors = FALSE)
  }))
}
flag_rows <- function(d) {
  d$poor_fit <- unname(poor_fit[d$site])
  d$noisy <- unname(noisy[d$site])
  d$flagged <- d$poor_fit | d$end_anomaly | d$noisy
  d
}

# Validation errors of the strategy, pooled by flag.
val_targets <- sort(unique(vb$month))
vb <- flag_rows(merge(vb, origin_state(val_targets), by = c("month", "site")))
vb$log_err <- log(vb$billed / vb$forecast)

# ---- Test months, read once ---------------------------------------------------
test_targets <- seq(config$test_start, max(bills$month), by = "month")
tb <- backtest(bills, test_targets, list(strategy = strategy, naive_per_day = naive_per_day_forecast),
               config$test_start, prep_fun, allow_test = TRUE)
tb <- flag_rows(merge(tb, origin_state(test_targets), by = c("month", "site")))
tb$log_err <- log(tb$billed / tb$forecast)
ts <- tb[tb$method == "strategy", ]

# Intervals of each construction on the test bills.
q_val <- empirical_quantiles(vb$log_err, vb$flagged, levels)
arima_sites <- names(site_method)[site_method == "arima_noseason"]
arima_d <- function(s) if (s %in% config$arima_season_sites) NA else 0
interval_rows <- function(d, q, train_of) {
  do.call(rbind, lapply(levels, function(l) {
    qq <- q[q$level == l, ]
    k <- match(as.character(d$flagged), qq$group)
    g <- gaussian_bounds(d$forecast, d$sigma, l)
    for (i in which(d$site %in% arima_sites)) {
      g[i, ] <- arima_bounds(train_of(d$month[i])[train_of(d$month[i])$site == d$site[i], ],
                             d$month[i], l, arima_d(d$site[i]))
    }
    rbind(data.frame(site = d$site, month = d$month, construction = "empirical", level = l,
                     lo = d$forecast * exp(qq$lo[k]), hi = d$forecast * exp(qq$hi[k])),
          data.frame(site = d$site, month = d$month, construction = "gaussian", level = l,
                     lo = g[, "lo"], hi = g[, "hi"]))
  }))
}
train_of <- function(t) bills[bills$month < t, ]
iv <- interval_rows(ts, q_val, train_of)
iv <- merge(iv, ts[c("site", "month", "billed", "flagged")], by = c("site", "month"))
cov <- stats::aggregate(cbind(inside = billed >= lo & billed <= hi) ~ construction + level, iv, mean)
cov_flag <- stats::aggregate(cbind(inside = billed >= lo & billed <= hi) ~ construction + level + flagged,
                             iv, mean)
gap <- tapply(abs(cov$inside - cov$level), cov$construction, mean)
construction <- names(which.min(gap))

# ---- Accuracy on the test months and the reliability flag ---------------------
test_acc <- rbind(cbind(scope = "ordinary", backtest_accuracy(tb[!tb$site %in% big, ])),
                  cbind(scope = "big", backtest_accuracy(tb[tb$site %in% big, ])))
val_mase <- c(ordinary = mean(vb$scaled_error[!vb$site %in% big]),
              big = mean(vb$scaled_error[vb$site %in% big]))
flag_mase <- stats::aggregate(scaled_error ~ flagged, ts, mean)
flag_n <- stats::aggregate(scaled_error ~ flagged, ts, length)
criterion_mase <- do.call(rbind, lapply(c("poor_fit", "end_anomaly", "noisy"), function(k) {
  data.frame(criterion = k, n = sum(ts[[k]]), mase_yes = mean(ts$scaled_error[ts[[k]]]),
             mase_no = mean(ts$scaled_error[!ts[[k]]]))
}))

# ---- December 2021 --------------------------------------------------------------
target <- config$forecast_month
prep <- prep_fun(bills)
dec <- strategy(prep, target)
names(dec)[names(dec) == "info"] <- "method"
dec$class <- unname(site_class[dec$site])
state <- flag_rows(origin_state(target))
dec <- merge(dec, state[c("site", "end_anomaly", "poor_fit", "noisy", "flagged", "sigma")],
             by = "site")
dec$forecast <- dec$kwh
all_err <- rbind(vb[c("log_err", "flagged")], ts[c("log_err", "flagged")])
q_all <- empirical_quantiles(all_err$log_err, all_err$flagged, levels)
dec_iv <- interval_rows(transform(dec, month = target), q_all, function(t) bills)
dec_iv <- dec_iv[dec_iv$construction == construction, ]
for (l in levels) {
  w <- dec_iv[dec_iv$level == l, ]
  tag <- sprintf("%02d", round(100 * l))
  dec[[paste0("lo", tag)]] <- w$lo[match(dec$site, w$site)]
  dec[[paste0("hi", tag)]] <- w$hi[match(dec$site, w$site)]
}
temp <- temporary_scenario(prep, target, config$edge_months)
last <- naive_per_day_forecast(prep, target)
dec$scenario_temporary <- ifelse(dec$end_anomaly, temp$kwh[match(dec$site, temp$site)], NA)
dec$scenario_lasting <- ifelse(dec$end_anomaly, last$kwh[match(dec$site, last$site)], NA)
dec$reasons <- apply(dec[c("poor_fit", "end_anomaly", "noisy")], 1,
                     function(r) paste(c("poor_fit", "end_anomaly", "noisy")[r], collapse = ", "))
dec$unreliable <- dec$end_anomaly | dec$noisy  # criteria that held on the test
dec <- dec[order(dec$site), c("site", "class", "method", "forecast", "lo95", "lo80", "hi80", "hi95",
                              "flagged", "unreliable", "reasons", "scenario_temporary",
                              "scenario_lasting")]
rownames(dec) <- NULL
stopifnot(nrow(dec) == length(unique(bills$site)), !anyNA(dec[c("forecast", "lo95", "hi95")]),
          all(dec$lo95 <= dec$lo80 & dec$lo80 <= dec$forecast & dec$forecast <= dec$hi80 &
                dec$hi80 <= dec$hi95))

# Conclusion (test months read once): ordinary sites MASE 1.16 for the
# strategy and the per-day naive (0.95 in validation); big sites, auto.arima
# 1.05 against 0.94 for the naive, its validation gain did not hold.
# Coverage at 80 / 95 %: empirical 70 / 86 %, Gaussian 86 / 91 %: Gaussian
# kept. The flag does not separate the sites in MASE (1.16 both), which
# divides by each site's usual change; in relative error flagged bills err
# more (mean 3.7 % against 2.2 %), through end_anomaly (median 6.3 % against
# 1.5 %) and noisy (2.6 % against 1.3 %); poor_fit adds nothing (1.3 %
# against 1.6 %). December 2021: 1,005,532 kWh in all, 4 sites with two
# scenarios.

# ---- Saved tables and numbers for the report ---------------------------------
write.csv(dec, file.path(config$tables_dir, "forecast_dec2021.csv"), row.names = FALSE)
results <- list(
  config = config[c("interval_levels", "reliability_mase_max", "reliability_noise_quantile",
                    "edge_months")],
  test_targets = test_targets,
  site_method = site_method,
  big = big,
  test = tb,
  test_accuracy = test_acc,
  val_mase = val_mase,
  quantiles_val = q_val,
  quantiles_all = q_all,
  coverage = cov,
  coverage_flag = cov_flag,
  coverage_gap = gap,
  construction = construction,
  flag_mase = merge(flag_mase, setNames(flag_n, c("flagged", "n")), by = "flagged"),
  criterion_mase = criterion_mase,
  noise_threshold = stats::quantile(noise[!names(noise) %in% big], config$reliability_noise_quantile),
  december = dec
)
saveRDS(results, file.path(config$results_dir, "forecast.rds"))
cat("Forecast written to", config$output_dir, "\n")
