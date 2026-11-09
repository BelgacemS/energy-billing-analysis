# Hampel filter by site (step 4): the reference detector the panel method of
# steps 5 to 7 has to beat.

#' Hampel filter on one series.
#'
#' Each point is compared with the median of the window of `half_width` points
#' on each side (truncated at the ends). The scale is the window MAD (x 1.4826,
#' consistent with the standard deviation under normality), floored by
#' `min_scale` so that a window of near-identical values does not flag every
#' small change. Returns the window median, the scale and the robust z-score.
hampel <- function(x, half_width, threshold, min_scale = 0) {
  n <- length(x)
  med <- scale <- numeric(n)
  for (i in seq_len(n)) {
    w <- x[max(1, i - half_width):min(n, i + half_width)]
    med[i] <- stats::median(w)
    scale[i] <- max(stats::mad(w), min_scale)
  }
  z <- (x - med) / scale
  data.frame(median = med, scale = scale, z = z, flag = abs(z) > threshold)
}

#' Hampel filter on log kWh per day of every site.
#'
#' The scale floor of a site is its robust noise level, the MAD of the first
#' differences divided by sqrt(2). The expected bill is the window median
#' times the days of the month; `gap_kwh` is billed minus expected (positive
#' means more billed than expected).
hampel_bills <- function(bills, half_width, threshold) {
  b <- bills[order(bills$site, bills$month), ]
  x <- log(b$kwh_per_day)
  groups <- split(seq_along(x), b$site)
  res <- do.call(rbind, lapply(groups, function(idx) {
    xs <- x[idx]
    hampel(xs, half_width, threshold, min_scale = stats::mad(diff(xs)) / sqrt(2))
  }))
  res <- res[order(unlist(groups)), ]  # back to the row order of b
  out <- cbind(b[c("site", "suffix", "month", "kwh", "month_days")], res)
  out$expected_kwh <- exp(out$median) * out$month_days
  out$gap_kwh <- out$kwh - out$expected_kwh
  rownames(out) <- NULL
  out
}
