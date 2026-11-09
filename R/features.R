# Descriptive features of the billing panel (step 3).

#' Sites whose median monthly bill exceeds `threshold` kWh.
big_sites <- function(bills, threshold) {
  med <- tapply(bills$kwh, bills$site, median)
  sort(names(med)[med > threshold])
}

#' Same date shifted by `k` months (dates are first days of months).
add_months <- function(month, k) {
  idx <- as.integer(format(month, "%Y")) * 12 + as.integer(format(month, "%m")) - 1 + k
  as.Date(sprintf("%04d-%02d-01", idx %/% 12, idx %% 12 + 1))
}

#' Log deviation of each bill from a local level of its site.
#'
#' centre = "year": median of the site's calendar year. "rolling": running
#' median over 13 months, which follows a level change better; it needs the
#' bills sorted by site and month with no gap. Both keep slow growth out of the
#' month effects. With "year", an incomplete year (2021 has 11 months) puts one
#' exact zero per site on whichever month is the median: restrict to complete
#' years when estimating a profile. `value` is "kwh" or "kwh_per_day".
seasonal_deviation <- function(bills, value = c("kwh_per_day", "kwh"),
                               centre = c("year", "rolling")) {
  value <- match.arg(value)
  centre <- match.arg(centre)
  x <- log(bills[[value]])
  year <- as.integer(format(bills$month, "%Y"))
  level <- if (centre == "year") {
    stats::ave(x, bills$site, year, FUN = median)
  } else {
    idx <- year * 12L + as.integer(format(bills$month, "%m"))
    same_site <- bills$site[-1] == bills$site[-nrow(bills)]
    if (is.unsorted(bills$site) || any(diff(idx)[same_site] != 1L)) {
      stop("Rolling centre needs bills sorted by site and month, without gaps")
    }
    stats::ave(x, bills$site, FUN = function(v) stats::runmed(v, 13, endrule = "median"))
  }
  data.frame(
    site   = bills$site,
    suffix = bills$suffix,
    year   = year,
    moy    = as.integer(format(bills$month, "%m")),
    dev    = x - level,
    stringsAsFactors = FALSE
  )
}

#' Calendar years with 12 bills for every site.
complete_years <- function(bills) {
  year <- as.integer(format(bills$month, "%Y"))
  per_site <- table(year) / length(unique(bills$site))
  as.integer(names(per_site)[per_site == 12])
}

#' Common month profile: median deviation over sites and years, per calendar
#' month (and per group if `by` names a column). `ratio` is exp(median).
common_profile <- function(dev, by = NULL) {
  keys <- dev[c(by, "moy")]
  out <- stats::aggregate(dev["dev"], keys, median)
  out$ratio <- exp(out$dev)
  out[do.call(order, out[c(by, "moy")]), , drop = FALSE]
}

#' Site-specific month effect: deviation minus the common profile.
specific_effect <- function(dev) {
  prof <- common_profile(dev)
  dev$effect <- dev$dev - prof$dev[match(dev$moy, prof$moy)]
  dev
}

#' Rank correlation across sites of the deviation of each calendar month
#' between two years. Near 0 means the site's own month effect does not repeat
#' from one year to the next. At a fixed month the common profile is a constant
#' shift, so it does not change the rank correlation and is not removed.
effect_reproducibility <- function(dev, year1, year2) {
  months <- intersect(dev$moy[dev$year == year1], dev$moy[dev$year == year2])
  rho <- vapply(sort(months), function(m) {
    a <- dev[dev$year == year1 & dev$moy == m, c("site", "dev")]
    b <- dev[dev$year == year2 & dev$moy == m, c("site", "dev")]
    ab <- merge(a, b, by = "site")
    stats::cor(ab$dev.x, ab$dev.y, method = "spearman")
  }, numeric(1))
  data.frame(moy = sort(months), year1 = year1, year2 = year2, rho = rho)
}

#' Largest level change of each site.
#'
#' For every month t, compares the medians of the `window` months before and
#' from t on (log scale). A median window ignores isolated spikes, so a large
#' value means a change that lasts at least window / 2 months. Around a clean
#' step several months give the same median gap; among them the step is the
#' one with the largest gap of means, which also holds at the edges.
max_level_shift <- function(mat, window) {
  if (anyNA(mat)) stop("Missing values in the matrix")
  n <- ncol(mat)
  breaks <- (window + 1):(n - window + 1)
  gap <- function(f) {
    g <- vapply(breaks, function(t) {
      apply(mat[, t:(t + window - 1), drop = FALSE], 1, f) -
        apply(mat[, (t - window):(t - 1), drop = FALSE], 1, f)
    }, numeric(nrow(mat)))
    matrix(g, nrow = nrow(mat))  # (sites, candidate break months)
  }
  shifts <- gap(median)
  mean_gaps <- gap(mean)
  best <- vapply(seq_len(nrow(mat)), function(i) {
    s <- abs(shifts[i, ])
    tied <- which(s >= max(s) - 1e-12)
    tied[which.max(abs(mean_gaps[i, tied]))]
  }, integer(1))
  data.frame(
    site        = rownames(mat),
    shift       = shifts[cbind(seq_len(nrow(mat)), best)],
    shift_month = colnames(mat)[breaks[best]],
    stringsAsFactors = FALSE
  )
}

#' Growth of each site: ratio of median kWh per day over the months of
#' `year2` and the same calendar months of `year1`, minus 1.
site_growth <- function(bills, year1, year2) {
  year <- as.integer(format(bills$month, "%Y"))
  moy <- format(bills$month, "%m")
  common_moy <- intersect(moy[year == year1], moy[year == year2])
  keep <- moy %in% common_moy
  med <- function(y) {
    sel <- keep & year == y
    tapply(log(bills$kwh_per_day[sel]), bills$site[sel], median)
  }
  m1 <- med(year1)
  m2 <- med(year2)[names(m1)]
  data.frame(site = names(m1), growth = exp(m2 - m1) - 1, stringsAsFactors = FALSE)
}

#' One-month-ahead naive forecast errors for each bill of `months`.
#'
#' raw: next bill = last bill. per_day: next bill = last kWh per day times the
#' days of the next month. Errors are relative to the actual bill.
naive_errors <- function(bills, months) {
  prev <- bills[c("site", "month", "kwh", "kwh_per_day")]
  prev$month <- add_months(prev$month, 1)
  names(prev)[3:4] <- c("prev_kwh", "prev_kwh_per_day")
  cur <- bills[bills$month %in% months, c("site", "month", "kwh", "month_days")]
  out <- merge(cur, prev, by = c("site", "month"))
  out$err_raw <- abs(out$kwh - out$prev_kwh) / out$kwh
  out$err_per_day <- abs(out$kwh - out$prev_kwh_per_day * out$month_days) / out$kwh
  out[order(out$site, out$month), c("site", "month", "err_raw", "err_per_day")]
}
