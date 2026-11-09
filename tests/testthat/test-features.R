# Small panel: 2 sites x 24 months, with known month effects.
make_panel <- function() {
  months <- seq(as.Date("2019-01-01"), by = "month", length.out = 24)
  days <- days_in_month(months)
  site_a <- 100 * days
  site_a[c(7, 19)] <- site_a[c(7, 19)] * 1.1  # July +10% each year
  site_b <- 50 * days
  data.frame(site = rep(c("A", "B"), each = 24), suffix = rep(c("BR", "KO"), each = 24),
             month = rep(months, 2), month_days = rep(days, 2),
             kwh = c(site_a, site_b), kwh_per_day = c(site_a, site_b) / rep(days, 2))
}

test_that("add_months crosses year boundaries", {
  d <- as.Date(c("2019-12-01", "2020-01-01", "2021-11-01"))
  expect_equal(add_months(d, 1), as.Date(c("2020-01-01", "2020-02-01", "2021-12-01")))
  expect_equal(add_months(d, -1), as.Date(c("2019-11-01", "2019-12-01", "2021-10-01")))
})

test_that("big_sites uses the median monthly bill", {
  b <- data.frame(site = c("A", "A", "A", "B", "B"), kwh = c(1, 1, 1e6, 5e5, 5e5))
  expect_equal(big_sites(b, 1e5), "B")
})

test_that("per-day deviation removes the calendar effect and keeps the real one", {
  b <- make_panel()
  dev <- seasonal_deviation(b, "kwh_per_day")
  expect_equal(dev$dev[dev$site == "B"], rep(0, 24))
  prof <- common_profile(dev)
  expect_equal(nrow(prof), 12)
  # Median of the two sites: (log(1.1) + 0) / 2 in July, 0 elsewhere.
  expect_equal(prof$dev[prof$moy == 7], log(1.1) / 2)
  expect_equal(prof$dev[prof$moy == 3], 0)
  # Raw kWh: February is below March by the ratio of their days.
  raw <- common_profile(seasonal_deviation(b, "kwh"))
  expect_equal(raw$ratio[2] / raw$ratio[3], 28 / 31, tolerance = 0.02)
})

test_that("effect_reproducibility is 1 for a site effect repeated each year", {
  months <- seq(as.Date("2019-01-01"), by = "month", length.out = 24)
  moy <- as.integer(format(months, "%m"))
  # Site k has month effect k * sin(moy), identical in both years.
  b <- do.call(rbind, lapply(1:3, function(k) {
    data.frame(site = paste0("S", k), suffix = "BR", month = months,
               kwh_per_day = exp(0.05 * k * sin(moy)))
  }))
  r <- effect_reproducibility(seasonal_deviation(b), 2019, 2020)
  expect_equal(nrow(r), 12)
  expect_equal(r$rho, rep(1, 12))
})

test_that("max_level_shift finds a step at the right month and ignores a spike", {
  x <- rep(0, 30)
  mat <- rbind(step = c(rep(0, 15), rep(log(0.5), 15)),
               spike = replace(x, 10, log(3)),
               flat = x)
  colnames(mat) <- format(seq(as.Date("2019-01-01"), by = "month", length.out = 30), "%Y-%m")
  s <- max_level_shift(mat, window = 6)
  expect_equal(s$shift[1], log(0.5))
  expect_equal(s$shift_month[1], colnames(mat)[16])
  expect_equal(s$shift[2:3], c(0, 0))
})

test_that("site_growth compares the same calendar months", {
  b <- make_panel()
  b$kwh_per_day[b$site == "A" & b$month >= as.Date("2020-01-01")] <- 120
  g <- site_growth(b, 2019, 2020)
  expect_equal(g$growth, c(0.2, 0))
})

test_that("naive_errors: per-day naive is exact on a constant daily rate", {
  b <- make_panel()
  e <- naive_errors(b[b$site == "B", ], as.Date(c("2019-02-01", "2019-03-01")))
  expect_equal(e$err_per_day, c(0, 0))
  expect_equal(e$err_raw, c(abs(28 - 31) / 28, abs(31 - 28) / 31))
})

test_that("rolling centre isolates a one-month effect on a stable site", {
  b <- make_panel()
  dev <- seasonal_deviation(b, "kwh_per_day", centre = "rolling")
  a <- dev[dev$site == "A", ]
  expect_equal(a$dev[c(7, 19)], rep(log(1.1), 2))
  expect_equal(a$dev[-c(7, 19)], rep(0, 22))
})
