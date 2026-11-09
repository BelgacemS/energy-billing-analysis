# Hand-made bills: site A is regular, site B has a partial month, a gap, an
# overlap and a repeated amount.
make_bills <- function() {
  rows <- data.frame(
    site  = c("A", "A", "A", "B", "B", "B", "B"),
    begin = as.Date(c("2020-01-01", "2020-02-01", "2020-03-01",
                      "2020-01-01", "2020-02-10", "2020-03-01", "2020-03-25")),
    end   = as.Date(c("2020-01-31", "2020-02-29", "2020-03-31",
                      "2020-01-31", "2020-02-29", "2020-03-31", "2020-04-30")),
    kwh   = c(100, 120, 120, 50, 40, 60, 55),
    stringsAsFactors = FALSE
  )
  rows$suffix <- "BR"
  rows$month <- as.Date(format(rows$begin, "%Y-%m-01"))
  rows$billed_days <- as.integer(rows$end - rows$begin) + 1L
  rows$month_days <- days_in_month(rows$month)
  rows$kwh_per_day <- rows$kwh / rows$month_days
  rows
}

test_that("check_periods flags bills not covering their calendar month", {
  p <- check_periods(make_bills())
  expect_equal(p$site, c("B", "B"))
  expect_equal(p$month, as.Date(c("2020-02-01", "2020-03-01")))
  expect_equal(p$missing_days, c(9L, -6L))  # 20 of 29 days; 37 days for a 31-day month
})

test_that("check_coverage finds gaps and overlaps with their length", {
  cov <- check_coverage(make_bills())
  expect_equal(cov$type, c("gap", "overlap"))
  expect_equal(cov$days, c(9L, -7L))  # 1-9 February uncovered; 25-31 March billed twice
  expect_equal(cov$next_begin, as.Date(c("2020-02-10", "2020-03-25")))
})

test_that("check_values and check_duplicates", {
  b <- make_bills()
  expect_equal(nrow(check_values(b)), 0)
  b$kwh[2] <- 0
  b$kwh[3] <- -5
  expect_equal(check_values(b)$kwh, c(0, -5))
  # Two B bills start in March (1 and 25): same site and month.
  d <- check_duplicates(make_bills())
  expect_equal(nrow(d), 2)
  expect_true(all(d$site == "B" & d$month == as.Date("2020-03-01")))
})

test_that("check_repeats counts zeros and the chance reference", {
  b <- make_bills()[1:3, ]  # A: changes +20 and 0
  r <- check_repeats(b, k_max = 20)
  expect_equal(r$observed, 1)
  expect_equal(r$repeats$month, as.Date("2020-03-01"))
  expect_equal(r$expected, 1 / 40)  # one change equal to +20 among 40 values of k
  expect_equal(r$p_value, 1 - exp(-1 / 40))
})

test_that("real data: 2 irregular periods, 2 gaps, 15 repeats on 12 sites", {
  b <- load_bills(root_path(config$data_path))
  p <- check_periods(b)
  expect_equal(paste(p$site, format(p$month, "%Y-%m")), c("0089KO 2020-02", "0628BR 2021-09"))
  cov <- check_coverage(b)
  expect_equal(cov$type, c("gap", "gap"))
  expect_equal(cov$days, c(1L, 14L))
  expect_equal(nrow(check_values(b)), 0)
  expect_equal(nrow(check_duplicates(b)), 0)
  r <- check_repeats(b, config$repeat_k_max)
  expect_equal(c(r$observed, r$n_sites), c(15, 12))
  expect_equal(r$expected, 12.55)
})
