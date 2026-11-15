# Extra forecasting candidates of step 10: per-day naive corrected by the
# month effect, robust local level, combination. Settings in R/config.R.

#' Month effect on log kWh per day from the complete calendar years of
#' `bills`: deviation from the site's median of the year, median over the
#' ordinary sites (group "common") or over the years of each big site (group =
#' site), centred to mean 0. Same format as month_effect (group, moy, effect).
#' With a single complete year, a big site's own effect would copy its
#' deviations exactly (flat corrected series, zero noise): the big sites then
#' take the common effect.
calendar_effect <- function(bills, big) {
  years <- complete_years(bills)
  if (length(years) == 0) stop("No complete calendar year to estimate the month effect")
  b <- bills[as.integer(format(bills$month, "%Y")) %in% years, ]
  dev <- seasonal_deviation(b, "kwh_per_day", "year")
  own <- if (length(years) >= 2) big else character()
  dev$group <- ifelse(dev$site %in% own, dev$site, "common")
  out <- stats::aggregate(dev ~ group + moy, dev, stats::median)
  for (s in setdiff(big, own)) {
    out <- rbind(out, transform(out[out$group == "common", ], group = s))
  }
  out$effect <- out$dev - stats::ave(out$dev, out$group)
  out <- out[order(out$group != "common", out$group, out$moy), c("group", "moy", "effect")]
  rownames(out) <- NULL
  out
}

#' Preprocessing of the backtest (refitted at each origin): the bills of
#' `train` sorted by site and month, plus the month effect of each bill
#' (`effect`), log kWh per day minus the effect (`y`), the same with the
#' months off their plateau replaced by the plateau level (`y_clean`) and the
#' first month of the last plateau after a break above `min_change`
#' (`from_last`, TRUE from that month on).
robust_preprocess <- function(big, seg_args, min_change) {
  function(train) {
    eff <- calendar_effect(train, big)
    pf <- panel_fit(train, eff, big, seg_args)
    b <- train[order(train$site, train$month), ]
    p <- pf$bills  # same order
    stopifnot(identical(b$site, p$site), identical(b$month, p$month))
    b$effect <- p$effect
    b$y <- p$y
    b$y_clean <- ifelse(p$excursion, p$level, p$y)
    large <- pf$breaks[abs(pf$breaks$change) > min_change, ]
    last <- tapply(large$month, large$site, max)
    start <- as.Date(unname(last[b$site]), origin = "1970-01-01")
    b$from_last <- is.na(start) | b$month >= start
    b
  }
}

# Effect of the target month of each site: the bill 12 months earlier has the
# same calendar month and the effect of the site's group.
target_effect <- function(prep, target) {
  ly <- prep[prep$month == add_months(target, -12), ]
  stats::setNames(ly$effect, ly$site)
}

#' Per-day naive corrected by the month effect (prep from robust_preprocess).
naive_adj_forecast <- function(prep, target) {
  last <- prep[prep$month == add_months(target, -1), ]
  f <- last$y + target_effect(prep, target)[last$site]
  data.frame(site = last$site, kwh = exp(unname(f)) * days_in_month(target),
             stringsAsFactors = FALSE)
}

#' Robust local level: simple exponential smoothing of y_clean from the
#' start of the last plateau, plus the target month effect.
robust_level_forecast <- function(prep, target) {
  eff <- target_effect(prep, target)
  rows <- lapply(split(prep, prep$site), function(d) {
    z <- d$y_clean[d$from_last]
    level <- tryCatch(as.numeric(suppressWarnings(forecast::ses(z, h = 1))$mean),
                      error = function(e) z[length(z)])
    data.frame(site = d$site[1], kwh = exp(level + eff[[d$site[1]]]) * days_in_month(target),
               stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

#' Equal-weight combination of the kWh forecasts of `forecasters` (named
#' list of forecasters sharing the same preprocessing).
combination_forecaster <- function(forecasters) {
  function(prep, target) {
    fc <- lapply(forecasters, function(f) {
      x <- f(prep, target)
      stats::setNames(x$kwh, x$site)
    })
    sites <- names(fc[[1]])
    kwh <- rowMeans(vapply(fc, function(x) x[sites], numeric(length(sites))))
    data.frame(site = sites, kwh = unname(kwh), stringsAsFactors = FALSE)
  }
}
