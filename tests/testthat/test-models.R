months <- function(n) seq(as.Date("2019-01-01"), by = "month", length.out = n)
one_site <- function(rate, site = "A") {
  m <- months(length(rate))
  data.frame(site = site, month = m, month_days = days_in_month(m),
             kwh = rate * days_in_month(m), kwh_per_day = rate, stringsAsFactors = FALSE)
}

test_that("site_ts is a monthly series of log kWh per day", {
  y <- site_ts(one_site(c(10, 20, 40))[2:3, ])
  expect_equal(stats::frequency(y), 12)
  expect_equal(stats::start(y), c(2019, 2))
  expect_equal(as.numeric(y), log(c(20, 40)))
})

test_that("model_forecaster back-transforms and falls back on errors", {
  b <- rbind(one_site(rep(10, 24)), one_site(rep(20, 24), "B"))
  target <- as.Date("2021-01-01")
  fixed <- model_forecaster(function(y, site) list(mean = log(5), info = "fixed"))
  out <- fixed(b, target)
  expect_equal(out$kwh, c(5, 5) * 31)
  expect_equal(out$info, c("fixed", "fixed"))
  broken <- model_forecaster(function(y, site) stop("no fit"))
  out <- broken(b, target)
  expect_equal(out$kwh, c(10, 20) * 31)  # per-day naive
  expect_equal(out$info, rep("fallback: no fit", 2))
  expect_error(fixed(b, as.Date("2021-02-01")), "precede the target")
})

test_that("the D = 0 variant never differences seasonally, except listed sites", {
  set.seed(3)
  rate <- exp(3 + 0.3 * sin(2 * pi * (1:36) / 12) + rnorm(36, sd = 0.01))
  y <- site_ts(one_site(rate))
  default <- make_fit_arima(NA)(y, "A")
  expect_match(default$info, "\\(\\d,1,\\d\\)\\[12\\]")  # seasonal difference chosen
  expect_match(make_fit_arima(0)(y, "A")$info, "\\(\\d,0,\\d\\)\\[12\\]")
  expect_equal(make_fit_arima(0, season_sites = "A")(y, "A")$info, default$info)
})

test_that("arima_label reports orders and drift", {
  set.seed(1)
  y <- stats::ts(cumsum(rnorm(30, mean = 1)), frequency = 12)
  fit <- forecast::Arima(y, order = c(0, 1, 1), include.drift = TRUE)
  expect_equal(arima_label(fit), "ARIMA(0,1,1)(0,0,0)[12] drift")
})

test_that("Ljung-Box: lags capped at T / 5, autocorrelation detected", {
  set.seed(2)
  wn <- stats::ts(rnorm(30), frequency = 12)
  lb <- arima_ljung_box(forecast::Arima(wn, order = c(0, 0, 0)))
  expect_equal(lb[["lag"]], 6)  # min(10, 30 / 5)
  ar <- stats::ts(stats::arima.sim(list(ar = 0.9), 30), frequency = 12)
  lb_ar <- arima_ljung_box(forecast::Arima(ar, order = c(0, 0, 0)))
  expect_lt(lb_ar[["p_value"]], 0.05)
  # too many parameters for the lags: no test
  expect_true(is.na(arima_ljung_box(forecast::Arima(ar, order = c(3, 0, 3)))[["p_value"]]))
})

test_that("arima_diagnosis flags degenerate orders and break-driven differencing", {
  set.seed(4)
  n <- 29
  y <- stats::ts(c(rep(0, 15), rep(3, 14)) + rnorm(n, sd = 0.3), frequency = 12)
  reg <- regimes(n, 16)
  expect_equal(reg, c(rep(1L, 15), rep(2L, 14)))
  fit <- forecast::Arima(y, order = c(0, 1, 1))
  dg <- arima_diagnosis(fit, y, reg, 0.05)
  expect_true(dg$break_diff)    # differencing only absorbs the step
  expect_false(dg$degenerate)
  expect_false(arima_diagnosis(fit, y, rep(1L, n), 0.05)$break_diff)  # no large break
  rw <- arima_diagnosis(forecast::Arima(y, order = c(0, 1, 0)), y, reg, 0.05)
  expect_true(rw$degenerate)
  expect_equal(rw$d, 1)
})
