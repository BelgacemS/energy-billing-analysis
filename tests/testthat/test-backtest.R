# Small hand-made panel: two sites, 24 months, known daily rates.
make_bills <- function(n = 24) {
  month <- seq(as.Date("2019-01-01"), by = "month", length.out = n)
  days <- days_in_month(month)
  rate <- list(A = 10 + seq_len(n), B = rep(c(20, 30), length.out = n))
  do.call(rbind, lapply(names(rate), function(s) {
    data.frame(site = s, suffix = "BR", month = month, month_days = days,
               kwh = rate[[s]] * days, kwh_per_day = rate[[s]], stringsAsFactors = FALSE)
  }))
}
baselines <- list(naive_raw = naive_raw_forecast, naive_per_day = naive_per_day_forecast,
                  snaive_per_day = snaive_per_day_forecast)
targets <- seq(as.Date("2020-07-01"), as.Date("2020-12-01"), by = "month")
test_start <- as.Date("2021-01-01")

test_that("baseline forecasters match hand computations", {
  b <- make_bills()
  target <- as.Date("2020-12-01")  # month 24; month 23 = November 2020, 30 days
  train <- b[b$month < target, ]
  expect_equal(naive_raw_forecast(train, target)$kwh, c(33 * 30, 20 * 30))
  expect_equal(naive_per_day_forecast(train, target)$kwh, c(33 * 31, 20 * 31))
  expect_equal(snaive_per_day_forecast(train, target)$kwh, c(22 * 31, 30 * 31))  # month 12
})

test_that("naive_scale is the in-sample MAE of the per-day naive", {
  b <- make_bills(4)  # January to April 2019
  sc <- naive_scale(b)
  d <- b$month_days[1:4]
  # site A: rate rises by 1 a day each month, so each error is the days of the month
  expect_equal(unname(sc["A"]), mean(d[2:4]))
  expect_equal(unname(sc["B"]), mean(10 * d[2:4]))
  flat <- b
  flat$kwh_per_day[flat$site == "A"] <- 5
  flat$kwh[flat$site == "A"] <- 5 * d
  expect_error(naive_scale(flat), "Zero MASE scale")
  expect_error(naive_scale(rbind(b, transform(b[1, ], site = "C"))), "No two consecutive bills")
})

test_that("backtest errors, relative and scaled errors are consistent", {
  b <- make_bills()
  bt <- backtest(b, targets, baselines, test_start)
  expect_equal(nrow(bt), 3 * length(targets) * 2)
  expect_equal(bt$error, bt$billed - bt$forecast)
  expect_equal(bt$rel_error, abs(bt$error) / bt$billed)
  dec <- bt[bt$month == as.Date("2020-12-01") & bt$site == "A" & bt$method == "naive_per_day", ]
  expect_equal(dec$error, 34 * 31 - 33 * 31)
  sc <- naive_scale(b[b$month <= as.Date("2020-11-01"), ])
  expect_equal(dec$scaled_error, 31 / unname(sc["A"]))
})

test_that("no leak: bills after the origin do not change its forecast", {
  b <- make_bills()
  seen <- as.Date(character())
  spy <- function(train) {  # preprocessing must only see the past
    seen <<- c(seen, max(train$month))
    train
  }
  bt <- backtest(b, targets, baselines, test_start, preprocess = spy)
  expect_equal(seen, add_months(targets, -1))
  for (target in as.list(targets)) {
    future <- b$month > add_months(target, -1)
    b2 <- b
    b2$kwh[future] <- b2$kwh[future] * 7
    b2$kwh_per_day[future] <- b2$kwh_per_day[future] * 7
    one <- function(x) {
      r <- backtest(x, target, baselines, test_start)
      r <- r[order(r$method, r$site), ]
      list(r$forecast, r$scale)
    }
    expect_equal(one(b2), one(b))
  }
})

test_that("test months are refused unless explicitly allowed", {
  b <- make_bills()
  expect_error(backtest(b, as.Date("2020-12-01"), baselines, as.Date("2020-12-01")), "step 11")
  expect_silent(backtest(b, as.Date("2020-12-01"), baselines, as.Date("2020-12-01"),
                         allow_test = TRUE))
})

test_that("a forecaster missing a site is an error", {
  b <- make_bills()
  bad <- list(bad = function(train, target) naive_raw_forecast(train, target)[1, ])
  expect_error(backtest(b, targets[1], bad, test_start), "one finite forecast per site")
})

test_that("backtest_accuracy summarises by method", {
  bt <- data.frame(method = c("m", "m", "m"), scaled_error = c(1, 2, 3),
                   rel_error = c(0.1, 0.2, 0.3), error = c(-5, 10, 20))
  acc <- backtest_accuracy(bt)
  expect_equal(acc$mase, 2)
  expect_equal(acc$rel_median, 0.2)
  expect_equal(acc$abs_kwh, 35)
  expect_equal(acc$signed_kwh, 25)
})

test_that("select_method keeps the simplest method within the tolerance", {
  ord <- c("naive_per_day", "naive_raw", "ets")
  expect_equal(select_method(c(ets = 0.90, naive_per_day = 1.00), ord, 0.02), "ets")
  expect_equal(select_method(c(ets = 0.99, naive_per_day = 1.00), ord, 0.02), "naive_per_day")
  expect_equal(select_method(c(ets = 1.00, naive_raw = 1.019), ord, 0.02), "naive_raw")
  expect_equal(select_method(c(ets = 1.00, naive_raw = 1.021), ord, 0.02), "ets")
  expect_error(select_method(c(other = 1), ord, 0.02), "method_order")
})
