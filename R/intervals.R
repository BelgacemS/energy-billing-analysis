# Forecast strategy, reliability flags and prediction intervals (step 11).
# Rules are in R/config.R.

#' Forecaster that applies to each site the method of `site_method` (named
#' by site), taken from `forecasters` (named list sharing one preprocessing).
strategy_forecaster <- function(site_method, forecasters) {
  function(prep, target) {
    out <- do.call(rbind, lapply(unique(site_method), function(m) {
      sites <- names(site_method)[site_method == m]
      fc <- forecasters[[m]](prep[prep$site %in% sites, ], target)
      data.frame(site = fc$site, kwh = fc$kwh, info = m, stringsAsFactors = FALSE)
    }))
    out[order(out$site), ]
  }
}

#' End-of-series anomaly of each site at the origin of a robust_preprocess
#' output: last bill off its plateau, or a break above the minimal change
#' (start of the last plateau, `from_last`) in the last `edge` months.
end_anomaly <- function(prep, edge) {
  vapply(split(prep, prep$site), function(d) {
    d <- d[order(d$month), ]
    n <- nrow(d)
    off <- d$y[n] != d$y_clean[n]
    start <- min(which(d$from_last))
    off || (start > 1 && start > n - edge)
  }, logical(1))
}

#' Temporary-gap scenario of each site for `target`: back to the level
#' before the anomaly. After a break above the minimal change in the last
#' `edge` months, the cleaned level of the month before that break; otherwise
#' the robust level (last bill set back to its plateau). Plus the target
#' month effect, on the kWh scale.
temporary_scenario <- function(prep, target, edge) {
  eff <- target_effect(prep, target)
  robust <- robust_level_forecast(prep, target)
  kwh <- stats::setNames(robust$kwh, robust$site)
  for (d in split(prep, prep$site)) {
    d <- d[order(d$month), ]
    n <- nrow(d)
    start <- min(which(d$from_last))
    if (start > 1 && start > n - edge) {
      kwh[[d$site[1]]] <- exp(d$y_clean[start - 1] + eff[[d$site[1]]]) * days_in_month(target)
    }
  }
  data.frame(site = names(kwh), kwh = unname(kwh), stringsAsFactors = FALSE)
}

#' Standard deviation of the monthly changes of log kWh per day of each
#' site (sigma of the per-day naive, a random walk on the log scale).
log_change_sd <- function(bills) {
  b <- bills[order(bills$site, bills$month), ]
  vapply(split(log(b$kwh_per_day), b$site), function(x) stats::sd(diff(x)), numeric(1))
}

#' Equal-tailed quantiles of the errors `err` (log billed / forecast) by
#' `group`, for each level: one row per group and level, lo and hi.
empirical_quantiles <- function(err, group, levels) {
  g <- split(err, group)
  do.call(rbind, lapply(names(g), function(k) {
    data.frame(group = k, level = levels,
               lo = stats::quantile(g[[k]], (1 - levels) / 2, names = FALSE),
               hi = stats::quantile(g[[k]], (1 + levels) / 2, names = FALSE),
               n = length(g[[k]]), stringsAsFactors = FALSE)
  }))
}

#' Gaussian bounds on the kWh scale: forecast x exp(+- z sigma), z the
#' normal quantile of the central `level`.
gaussian_bounds <- function(forecast, sigma, level) {
  z <- stats::qnorm((1 + level) / 2)
  cbind(lo = forecast * exp(-z * sigma), hi = forecast * exp(z * sigma))
}

#' auto.arima interval of one site (log kWh per day, D as in make_fit_arima)
#' for the month after the last bill of `d`, on the kWh scale.
arima_bounds <- function(d, target, level, D) {
  m <- suppressWarnings(forecast::auto.arima(site_ts(d[order(d$month), ]), D = D))
  f <- forecast::forecast(m, h = 1, level = 100 * level)
  days <- days_in_month(target)
  c(lo = exp(as.numeric(f$lower)) * days, hi = exp(as.numeric(f$upper)) * days)
}

#' Share of bills inside their interval.
coverage <- function(billed, lo, hi) mean(billed >= lo & billed <= hi)
