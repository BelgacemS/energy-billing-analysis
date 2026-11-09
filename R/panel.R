# Expected value of each bill by a robust panel (step 5): plateau of the site
# plus month effect, on log kWh per calendar day.

#' Month effect on log kWh per day, centred to mean 0 over the 12 months.
#'
#' Deviation of each bill from the trend of the classical additive
#' decomposition (stats::decompose): centred 2x12 moving average, which removes
#' a linear trend and gives no exact zeros (a rolling median does, for 8 to 37 %
#' of the bills, and pins the median of a small effect at 0). The first and
#' last 6 months of each site have no trend and are left out. The ordinary
#' sites share one effect (median over sites and years, group "common"); each
#' site of `big` gets its own (group = site name). The average smooths a break
#' over 12 months; the median over sites limits the bias as long as the sites
#' do not all break in the same direction around the same month.
month_effect <- function(bills, big) {
  b <- bills[order(bills$site, bills$month), ]
  x <- log(b$kwh_per_day)
  weights <- c(0.5, rep(1, 11), 0.5) / 12
  trend <- stats::ave(x, b$site, FUN = function(v) as.numeric(stats::filter(v, weights, sides = 2)))
  dev <- data.frame(site = b$site, moy = as.integer(format(b$month, "%m")), dev = x - trend,
                    stringsAsFactors = FALSE)
  dev <- dev[!is.na(dev$dev), ]
  dev$group <- ifelse(dev$site %in% big, dev$site, "common")
  out <- stats::aggregate(dev ~ group + moy, dev, stats::median)
  out$effect <- out$dev - stats::ave(out$dev, out$group)
  out <- out[order(out$group != "common", out$group, out$moy), c("group", "moy", "effect")]
  rownames(out) <- NULL
  out
}

#' Expected value, residual and robust z-score of every bill.
#'
#' For site i and month t, with x = log kWh per calendar day and s the month
#' effect of the site's group: y = x - s is segmented by segment_series, and
#' the expected value is level + s. resid = x - expected; z = resid / sigma,
#' sigma being the noise scale of y. `seg_args` is a list of arguments passed
#' to segment_series. Returns the bills table (`bills`) and the lasting
#' breaks of every site (`breaks`, dated at the first month of the new
#' plateau, change = relative change of kWh per day).
panel_fit <- function(bills, effect, big, seg_args = list()) {
  b <- bills[order(bills$site, bills$month), ]
  b$x <- log(b$kwh_per_day)
  b$group <- ifelse(b$site %in% big, b$site, "common")
  moy <- as.integer(format(b$month, "%m"))
  b$effect <- effect$effect[match(paste(b$group, moy), paste(effect$group, effect$moy))]
  if (anyNA(b$effect)) stop("Month effect missing for some bills")
  b$y <- b$x - b$effect

  fits <- lapply(split(b$y, b$site), function(y) do.call(segment_series, c(list(y), seg_args)))
  months <- do.call(rbind, lapply(fits, `[[`, "months"))
  b <- cbind(b, months)  # split() and b are both ordered by site, then month
  b$sigma <- unname(vapply(fits, `[[`, numeric(1), "sigma")[b$site])

  b$expected_log <- b$level + b$effect
  b$expected_kwh <- exp(b$expected_log) * b$month_days
  b$gap_kwh <- b$kwh - b$expected_kwh
  b$resid <- b$x - b$expected_log
  b$z <- b$resid / b$sigma

  site_months <- split(b$month, b$site)
  breaks <- do.call(rbind, lapply(names(fits), function(s) {
    br <- fits[[s]]$breaks
    if (nrow(br) == 0) return(NULL)
    data.frame(site = s, month = site_months[[s]][br$index],
               level_before = br$level_before, level_after = br$level_after,
               change = exp(br$level_after - br$level_before) - 1,
               stringsAsFactors = FALSE)
  }))
  if (is.null(breaks)) {
    breaks <- data.frame(site = character(), month = as.Date(character()),
                         level_before = numeric(), level_after = numeric(),
                         change = numeric())
  }
  rownames(b) <- NULL
  keep <- c("site", "suffix", "group", "month", "kwh", "month_days", "x", "effect", "y",
            "segment", "excursion", "plateau", "level", "sigma", "expected_log",
            "expected_kwh", "gap_kwh", "resid", "z")
  list(bills = b[keep], breaks = breaks)
}
