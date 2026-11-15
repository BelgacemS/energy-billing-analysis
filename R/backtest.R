# Rolling-origin backtest, one month ahead (step 8), and the baseline
# forecasters. A forecaster is a function(train, target) returning one row per
# site of `train` with columns site and kwh (forecast of the bill of `target`).

#' Raw naive: the next bill equals the last one.
naive_raw_forecast <- function(train, target) {
  last <- train[train$month == add_months(target, -1), ]
  data.frame(site = last$site, kwh = last$kwh, stringsAsFactors = FALSE)
}

#' Seasonal naive per day: kWh per day of the same month a year earlier
#' times the days of the target month.
snaive_per_day_forecast <- function(train, target) {
  ly <- train[train$month == add_months(target, -12), ]
  data.frame(site = ly$site, kwh = ly$kwh_per_day * days_in_month(target),
             stringsAsFactors = FALSE)
}

#' Scale of the MASE of each site: mean absolute one-step error of the
#' per-day naive over the bills of `train` (in sample, kWh).
naive_scale <- function(train) {
  b <- train[order(train$site, train$month), ]
  n <- nrow(b)
  same <- b$site[-1] == b$site[-n] & b$month[-1] == add_months(b$month[-n], 1)
  err <- abs(b$kwh[-1] - b$kwh_per_day[-n] * b$month_days[-1])[same]
  scale <- tapply(err, b$site[-1][same], mean)
  if (!all(unique(b$site) %in% names(scale))) stop("No two consecutive bills for some sites")
  if (any(!(scale > 0))) stop("Zero MASE scale for some sites")
  scale
}

#' One-month-ahead forecasts of every target month with every method.
#'
#' For a target month t, the origin is t - 1: `train` holds the bills up to
#' the origin only, `preprocess(train)` is refitted on it and passed to each
#' forecaster of `methods` (named list). Test months (t >= test_start) are
#' refused unless allow_test = TRUE (step 11). Returns one row per method,
#' target month and site: billed, forecast, error = billed - forecast,
#' relative error |error| / billed, MASE scale (naive_scale) and scaled error
#' |error| / scale.
backtest <- function(bills, targets, methods, test_start, preprocess = identity,
                     allow_test = FALSE) {
  if (!allow_test && any(targets >= test_start)) {
    stop("Test months are read once, at step 11 (allow_test = TRUE)")
  }
  rows <- lapply(targets, function(target) {
    train <- bills[bills$month <= add_months(target, -1), ]
    scale <- naive_scale(train)
    prep <- preprocess(train)
    actual <- bills[bills$month == target, c("site", "kwh")]
    do.call(rbind, lapply(names(methods), function(m) {
      fc <- methods[[m]](prep, target)[c("site", "kwh")]
      if (anyDuplicated(fc$site) || !setequal(fc$site, actual$site) ||
          any(!is.finite(fc$kwh))) {
        stop("Method ", m, ": one finite forecast per site expected for ", target)
      }
      out <- merge(actual, fc, by = "site", suffixes = c("", "_fc"))
      err <- out$kwh - out$kwh_fc
      data.frame(method = m, origin = add_months(target, -1), month = target,
                 site = out$site, billed = out$kwh, forecast = out$kwh_fc, error = err,
                 rel_error = abs(err) / out$kwh, scale = unname(scale[out$site]),
                 scaled_error = abs(err) / unname(scale[out$site]),
                 stringsAsFactors = FALSE)
    }))
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

#' Accuracy of each method over the rows of a backtest: mean MASE, median,
#' mean and 90th percentile of the relative error, total absolute and signed
#' error in kWh (signed > 0: forecasts too low overall).
backtest_accuracy <- function(bt) {
  do.call(rbind, lapply(split(bt, bt$method), function(d) {
    data.frame(method = d$method[1], n = nrow(d), mase = mean(d$scaled_error),
               rel_median = stats::median(d$rel_error), rel_mean = mean(d$rel_error),
               rel_p90 = unname(stats::quantile(d$rel_error, 0.9)),
               abs_kwh = sum(abs(d$error)), signed_kwh = sum(d$error),
               stringsAsFactors = FALSE)
  }))
}

#' Selection rule of R/config.R: among the methods whose mean MASE is less
#' than `tolerance` (relative) above the lowest, the first one in `order`.
select_method <- function(mase, order, tolerance) {
  if (is.null(names(mase)) || !all(names(mase) %in% order)) {
    stop("Every method must be named and ranked in method_order")
  }
  ok <- names(mase)[mase < min(mase) * (1 + tolerance)]
  order[order %in% ok][1]
}
