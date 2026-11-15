# Course models fitted site by site (step 9): simple exponential smoothing,
# ETS and auto.arima on log kWh per calendar day, and the diagnosis of
# auto.arima. Settings are in R/config.R.

#' Monthly ts of log kWh per calendar day of one site (bills sorted by month).
site_ts <- function(d) {
  start <- as.integer(c(format(d$month[1], "%Y"), format(d$month[1], "%m")))
  stats::ts(log(d$kwh_per_day), start = start, frequency = 12)
}

#' Text label of an Arima fit, e.g. "ARIMA(0,1,1)(0,0,0)[12] drift".
arima_label <- function(fit) {
  o <- forecast::arimaorder(fit)
  if (length(o) < 7) o <- c(o, P = 0, D = 0, Q = 0, Frequency = 12)
  extra <- intersect(c("drift", "intercept"), names(stats::coef(fit)))
  paste0(sprintf("ARIMA(%d,%d,%d)(%d,%d,%d)[12]", o[["p"]], o[["d"]], o[["q"]], o[["P"]],
                 o[["D"]], o[["Q"]]), if (length(extra)) paste0(" ", extra[1]) else "")
}

#' Fitting functions: each takes a site ts and its name, and returns the
#' one-step forecast on the log scale (`mean`) and a label (`info`).
fit_ses <- function(y, site) {
  f <- forecast::ses(y, h = 1)
  list(mean = as.numeric(f$mean), info = sprintf("SES alpha=%.4f", f$model$par[["alpha"]]))
}
make_fit_ets <- function(additive_only) {
  function(y, site) {
    m <- forecast::ets(y, additive.only = additive_only)
    list(mean = as.numeric(forecast::forecast(m, h = 1)$mean), info = m$method)
  }
}
#' auto.arima with D fixed (NA: chosen by the seasonal test), except for the
#' sites of `season_sites`, which keep the default choice.
make_fit_arima <- function(D = NA, season_sites = character()) {
  function(y, site) {
    m <- forecast::auto.arima(y, D = if (site %in% season_sites) NA else D)
    list(mean = as.numeric(forecast::forecast(m, h = 1)$mean), info = arima_label(m))
  }
}

#' Forecaster (see R/backtest.R) that fits `fit_fun` on each site of `train`.
#'
#' Forecast = exp(log forecast) x days of the target month. A fit that fails
#' falls back to the per-day naive, with info "fallback: <message>".
model_forecaster <- function(fit_fun) {
  function(train, target) {
    rows <- lapply(split(train, train$site), function(d) {
      d <- d[order(d$month), ]
      if (add_months(max(d$month), 1) != target) stop("The last bill must precede the target")
      res <- tryCatch(suppressWarnings(fit_fun(site_ts(d), d$site[1])),
                      error = function(e) list(mean = log(d$kwh_per_day[nrow(d)]),
                                               info = paste("fallback:", conditionMessage(e))))
      data.frame(site = d$site[1], kwh = exp(res$mean) * days_in_month(target),
                 info = res$info, stringsAsFactors = FALSE)
    })
    out <- do.call(rbind, rows)
    rownames(out) <- NULL
    out
  }
}

#' Ljung-Box test of the residuals of an Arima fit: h = 10 lags (24 for a
#' seasonal model) capped at T / 5, fitdf = p + q + P + Q. NA when h <= fitdf.
arima_ljung_box <- function(fit) {
  o <- forecast::arimaorder(fit)
  get <- function(k) if (k %in% names(o)) o[[k]] else 0
  fitdf <- get("p") + get("q") + get("P") + get("Q")
  seasonal <- get("P") + get("D") + get("Q") > 0
  r <- stats::residuals(fit)
  h <- min(if (seasonal) 24 else 10, floor(length(r) / 5))
  if (h <= fitdf) return(c(lag = h, p_value = NA_real_))
  c(lag = h, p_value = stats::Box.test(r, lag = h, type = "Ljung-Box", fitdf = fitdf)$p.value)
}

#' Kind of an auto.arima fit from its label (arima_label): "naive" (random
#' walk without drift, the per-day naive itself), "seasonal_diff" (D = 1),
#' "drift", "mean" (stationary model with a mean) or "other_arma".
arima_kind <- function(info) {
  ifelse(grepl("^ARIMA\\(0,1,0\\)\\(0,0,0\\)\\[12\\]$", info), "naive",
  ifelse(grepl("\\(\\d,1,\\d\\)\\[12\\]", info), "seasonal_diff",
  ifelse(grepl("drift", info), "drift",
  ifelse(grepl("intercept", info), "mean", "other_arma"))))
}

#' Regime index of each month between level breaks: `break_index` gives the
#' first month of each new regime.
regimes <- function(n, break_index) cumsum(seq_len(n) %in% break_index) + 1L

#' Diagnosis of one auto.arima fit on series `y` (see R/config.R): order,
#' degenerate order, differencing caused by breaks (`regime`: regime of each
#' month between large breaks), Ljung-Box p-value.
arima_diagnosis <- function(fit, y, regime, alpha) {
  o <- forecast::arimaorder(fit)
  get <- function(k) if (k %in% names(o)) o[[k]] else 0L
  d <- get("d")
  D <- get("D")
  degenerate <- get("p") + get("q") + get("P") + get("Q") == 0
  break_diff <- FALSE
  if (length(unique(regime)) > 1 && d + D > 0) {
    r <- stats::ts(stats::residuals(stats::lm(as.numeric(y) ~ factor(regime))),
                   start = stats::start(y), frequency = 12)
    break_diff <- (d > 0 && forecast::ndiffs(r) < d) ||
      (D > 0 && forecast::nsdiffs(r) == 0)
  }
  lb <- arima_ljung_box(fit)
  data.frame(model = arima_label(fit), d = d, D = D, degenerate = degenerate,
             break_diff = break_diff, lb_lag = lb[["lag"]], lb_p = lb[["p_value"]],
             lb_reject = !is.na(lb[["p_value"]]) && lb[["p_value"]] < alpha,
             stringsAsFactors = FALSE)
}

#' Share of random walks flagged as break-driven differencing: `reps` walks
#' of length `n` whose increments come from `draw(n)`, segmented and
#' diagnosed as the sites are (D = 0 variant). Measures how often the
#' criterion fires when differencing is the right model.
break_diff_null_rate <- function(draw, n, reps, seg_args, min_change, alpha) {
  hits <- vapply(seq_len(reps), function(r) {
    y <- stats::ts(cumsum(draw(n)), frequency = 12)
    seg <- do.call(segment_series, c(list(as.numeric(y)), seg_args))
    br <- seg$breaks
    large <- br$index[abs(exp(br$level_after - br$level_before) - 1) > min_change]
    fit <- suppressWarnings(forecast::auto.arima(y, D = 0))
    arima_diagnosis(fit, y, regimes(n, large), alpha)$break_diff
  }, logical(1))
  mean(hits)
}
