# Loading and normalisation of the raw bills.

#' Parse kWh values written with a space as thousands separator ("2 500").
parse_kwh <- function(x) {
  cleaned <- gsub("[[:space:]  ]", "", x)
  out <- suppressWarnings(as.numeric(cleaned))
  bad <- is.na(out)
  if (any(bad)) {
    stop("Unparsable kWh values: ", paste(head(unique(x[bad]), 5), collapse = ", "))
  }
  out
}

#' Number of days in the calendar month of each date.
days_in_month <- function(date) {
  y <- as.integer(format(date, "%Y"))
  m <- as.integer(format(date, "%m"))
  first <- as.Date(sprintf("%04d-%02d-01", y, m))
  next_first <- as.Date(sprintf("%04d-%02d-01", y + (m == 12), m %% 12 + 1))
  as.integer(next_first - first)
}

#' Read the raw CSV and return one row per bill.
#'
#' The billing month is the month of BEGIN. kwh_per_day divides by the days of
#' the calendar month, not by the billed days, so a partial period shows up as
#' a low value instead of being rescaled (irregular periods are handled by the
#' checks of step 4).
load_bills <- function(path) {
  raw <- utils::read.csv(path, colClasses = "character", encoding = "UTF-8",
                         check.names = FALSE)
  expected <- c("NAME", "BEGIN", "END", "KWH_PERIOD")
  missing_cols <- setdiff(expected, names(raw))
  if (length(missing_cols) > 0) {
    stop("Missing columns: ", paste(missing_cols, collapse = ", "))
  }

  begin <- as.Date(raw$BEGIN, format = "%d/%m/%Y")
  end <- as.Date(raw$END, format = "%d/%m/%Y")
  if (anyNA(begin) || anyNA(end)) stop("Unparsable dates in BEGIN or END")

  site <- trimws(raw$NAME)
  month <- as.Date(format(begin, "%Y-%m-01"))
  kwh <- parse_kwh(raw$KWH_PERIOD)
  month_days <- days_in_month(month)

  bills <- data.frame(
    site        = site,
    suffix      = substr(site, nchar(site) - 1, nchar(site)),
    month       = month,
    begin       = begin,
    end         = end,
    billed_days = as.integer(end - begin) + 1L,
    month_days  = month_days,
    kwh         = kwh,
    kwh_per_day = kwh / month_days,
    stringsAsFactors = FALSE
  )
  bills <- bills[order(bills$site, bills$month), ]
  rownames(bills) <- NULL
  bills
}

#' Site x month matrix of log(kWh per calendar day).
#'
#' Columns cover every month from the first to the last billing month, so a
#' missing bill appears as NA. Duplicated site-months and non-positive values
#' are errors: the log scale is undefined for them.
log_kwh_matrix <- function(bills) {
  key <- paste(bills$site, bills$month)
  if (anyDuplicated(key)) {
    stop("Several bills for the same site and month: ",
         paste(head(unique(key[duplicated(key)]), 5), collapse = ", "))
  }
  if (any(bills$kwh <= 0)) stop("Non-positive kWh values, log undefined")

  months <- seq(min(bills$month), max(bills$month), by = "month")
  sites <- sort(unique(bills$site))
  mat <- matrix(NA_real_, nrow = length(sites), ncol = length(months),
                dimnames = list(sites, format(months, "%Y-%m")))  # (sites, months)
  mat[cbind(match(bills$site, sites), match(bills$month, months))] <-
    log(bills$kwh_per_day)
  mat
}
