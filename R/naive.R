# Naive forecasts (step 4, reused by the backtest of step 8).

#' Per-day naive forecast of `target_month` for every site.
#'
#' Last kWh per calendar day of the site times the days of the target month.
#' The last bill must be the month just before the target.
naive_per_day_forecast <- function(bills, target_month) {
  last_month <- max(bills$month)
  if (add_months(last_month, 1) != target_month) {
    stop("The last billed month must be the month before the target")
  }
  last <- bills[bills$month == last_month, c("site", "suffix", "kwh_per_day")]
  out <- data.frame(
    site   = last$site,
    suffix = last$suffix,
    month  = target_month,
    method = "naive_per_day",
    kwh    = last$kwh_per_day * days_in_month(target_month),
    stringsAsFactors = FALSE
  )
  out <- out[order(out$site), ]
  rownames(out) <- NULL
  out
}
