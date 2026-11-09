# Deterministic checks on the bills (step 4). Each returns the offending rows.

#' Bills whose period is not the full calendar month of BEGIN.
check_periods <- function(bills) {
  last_day <- add_months(bills$month, 1) - 1
  bad <- bills$begin != bills$month | bills$end != last_day
  out <- bills[bad, c("site", "month", "begin", "end", "billed_days", "month_days", "kwh")]
  out$missing_days <- out$month_days - out$billed_days
  rownames(out) <- NULL
  out
}

#' Gaps and overlaps between consecutive bills of a site.
#'
#' `days` is the number of uncovered days (gap, > 0) or of days billed twice
#' (overlap, < 0) between the end of a bill and the start of the next one.
#' A series that starts late or stops early is not detected here; on these
#' data the step 2 test of a complete 92 x 35 panel rules it out.
check_coverage <- function(bills) {
  b <- bills[order(bills$site, bills$begin), ]
  same_site <- b$site[-1] == b$site[-nrow(b)]
  days <- as.integer(b$begin[-1] - b$end[-nrow(b)]) - 1L
  bad <- which(same_site & days != 0)
  data.frame(
    site       = b$site[bad + 1],
    prev_end   = b$end[bad],
    next_begin = b$begin[bad + 1],
    days       = days[bad],
    type       = ifelse(days[bad] > 0, "gap", "overlap"),
    stringsAsFactors = FALSE
  )
}

#' Bills with a missing, zero or negative amount.
check_values <- function(bills) {
  bills[is.na(bills$kwh) | bills$kwh <= 0, c("site", "month", "kwh")]
}

#' Bills sharing a site and month with another bill.
check_duplicates <- function(bills) {
  key <- paste(bills$site, bills$month)
  bills[key %in% key[duplicated(key)], c("site", "month", "begin", "end", "kwh")]
}

#' Identical amounts two months in a row, and how often that happens by chance.
#'
#' The month-to-month change of a bill is spread smoothly around 0, so the
#' number of changes exactly equal to 0 should be close to the number equal to
#' any other small signed value k. The reference is the mean count of changes
#' equal to k for k in +-1..+-`k_max` kWh, and the p-value is the Poisson
#' probability of seeing at least the observed number of zeros.
#' Consecutive bills of a site are assumed to be consecutive months, as in
#' the complete panel tested at step 2.
check_repeats <- function(bills, k_max = 10) {
  b <- bills[order(bills$site, bills$month), ]
  same_site <- b$site[-1] == b$site[-nrow(b)]
  change <- diff(b$kwh)[same_site]
  rows <- which(same_site & diff(b$kwh) == 0) + 1
  ks <- c(-k_max:-1, 1:k_max)
  expected <- mean(vapply(ks, function(k) sum(change == k), numeric(1)))
  observed <- length(rows)
  list(
    repeats  = b[rows, c("site", "month", "kwh")],
    observed = observed,
    n_sites  = length(unique(b$site[rows])),
    expected = expected,
    p_value  = stats::ppois(observed - 1, expected, lower.tail = FALSE)
  )
}
