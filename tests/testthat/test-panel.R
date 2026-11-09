make_bills <- function(rates, sites) {
  months <- seq(as.Date("2019-01-01"), by = "month", length.out = length(rates) / length(sites))
  b <- data.frame(site = rep(sites, each = length(months)), suffix = "BR",
                  month = rep(months, length(sites)), stringsAsFactors = FALSE)
  b$month_days <- days_in_month(b$month)
  b$kwh_per_day <- rates
  b$kwh <- b$kwh_per_day * b$month_days
  b
}

test_that("month_effect is centred and separates the big sites", {
  set.seed(1)
  n <- 36
  season <- rep(c(rep(-0.02, 6), rep(0.02, 6)), 3)
  rates <- c(exp(log(10) + season + rnorm(n, sd = 0.001)),
             exp(log(20) + season + rnorm(n, sd = 0.001)),
             exp(log(1e4) - 2 * season + rnorm(n, sd = 0.001)))
  b <- make_bills(rates, c("A", "B", "Z"))
  eff <- month_effect(b, big = "Z")
  expect_setequal(unique(eff$group), c("common", "Z"))
  expect_equal(as.numeric(tapply(eff$effect, eff$group, mean)), c(0, 0), tolerance = 1e-12)
  common <- eff$effect[eff$group == "common"]
  big <- eff$effect[eff$group == "Z"]
  expect_lt(mean(common[1:6]), 0)
  expect_gt(mean(big[1:6]), 0)
})

test_that("month_effect recovers a small season despite breaks and growth", {
  set.seed(3)
  n <- 35
  season <- 0.02 * sin(2 * pi * (1:n) / 12)  # starts in January 2019
  sites <- sprintf("S%02d", 1:40)
  rates <- unlist(lapply(seq_along(sites), function(i) {
    lvl <- log(1000 + 50 * i) + 0.002 * (1:n)           # slow growth
    lvl[(5 + i %% 25):n] <- lvl[(5 + i %% 25):n] + 0.1  # one break per site
    exp(lvl + season + rnorm(n, sd = 0.005))
  }))
  eff <- month_effect(make_bills(rates, sites), big = character())
  truth <- 0.02 * sin(2 * pi * (1:12) / 12)
  expect_equal(eff$effect, truth - mean(truth), tolerance = 0.004, scale = 1)
})

test_that("panel_fit gives the expected bill, the gap in kWh and the break", {
  set.seed(2)
  n <- 35
  ra <- exp(log(10) + rnorm(n, sd = 0.01))
  ra[12] <- ra[12] * 0.6                   # A: one bill 40 % low
  rb <- exp(log(20) + rnorm(n, sd = 0.01))
  rb[25:n] <- rb[25:n] * 1.5               # B: lasting +50 % from month 25
  b <- make_bills(c(ra, rb), c("A", "B"))
  eff <- data.frame(group = "common", moy = 1:12, effect = 0)
  fit <- panel_fit(b, eff, big = character())
  p <- fit$bills
  expect_equal(nrow(p), nrow(b))
  a12 <- p[p$site == "A", ][12, ]
  expect_true(a12$excursion)
  expect_lt(a12$z, -10)
  expect_equal(a12$gap_kwh, a12$kwh - a12$expected_kwh)
  expect_equal(a12$expected_kwh / a12$month_days, 10, tolerance = 0.02)
  expect_equal(fit$breaks$site, "B")
  expect_equal(fit$breaks$month, b$month[25])
  expect_equal(fit$breaks$change, 0.5, tolerance = 0.03)
})
