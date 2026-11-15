season <- 0.03 * sin(2 * pi * (1:12) / 12)  # log effect of each calendar month

sim_bills <- function(n = 30, n_sites = 6, sd = 0.01, seed = 1, edit = NULL) {
  set.seed(seed)
  month <- seq(as.Date("2019-01-01"), by = "month", length.out = n)
  days <- days_in_month(month)
  moy <- as.integer(format(month, "%m"))
  out <- do.call(rbind, lapply(seq_len(n_sites), function(i) {
    x <- log(50 * i) + season[moy] + stats::rnorm(n, sd = sd)
    if (!is.null(edit) && i == 1) x <- edit(x)
    data.frame(site = sprintf("S%d", i), suffix = "KO", month = month, month_days = days,
               kwh = exp(x) * days, kwh_per_day = exp(x), stringsAsFactors = FALSE)
  }))
  out
}
seg <- list(cost = "biweight", penalty = "MBIC", plateau_min = 3, K = 6)
prep_fun <- robust_preprocess(character(), seg, 0.10)

test_that("calendar_effect recovers a shared month profile from complete years", {
  eff <- calendar_effect(sim_bills(sd = 0.001), character())
  expect_equal(nrow(eff), 12)
  expect_equal(sum(eff$effect), 0, tolerance = 1e-12)
  expect_equal(eff$effect, season - mean(season), tolerance = 0.005)
  expect_error(calendar_effect(sim_bills(n = 11), character()), "No complete calendar year")
})

test_that("naive_adj moves the last value by the month effect", {
  b <- sim_bills(sd = 1e-5)
  target <- add_months(max(b$month), 1)
  prep <- prep_fun(b)
  f <- naive_adj_forecast(prep, target)
  m <- as.integer(format(target, "%m"))
  expect_equal(f$kwh, 50 * seq_len(6) * exp(season[m]) * days_in_month(target), tolerance = 1e-4)
})

test_that("robust_level ignores a last-month spike and restarts at the last plateau", {
  spike <- sim_bills(edit = function(x) { x[30] <- x[30] + log(2); x })
  target <- add_months(max(spike$month), 1)
  f <- robust_level_forecast(prep_fun(spike), target)
  naive <- naive_per_day_forecast(spike, target)
  expect_lt(abs(log(f$kwh[1] / (50 * days_in_month(target)))), 0.05)
  expect_gt(naive$kwh[1] / f$kwh[1], 1.8)
  step <- sim_bills(edit = function(x) { x[18:30] <- x[18:30] + log(1.5); x })
  f <- robust_level_forecast(prep_fun(step), target)
  expect_lt(abs(log(f$kwh[1] / (75 * days_in_month(target)))), 0.05)
})

test_that("combination is the mean of its forecasters", {
  b <- sim_bills()
  target <- add_months(max(b$month), 1)
  prep <- prep_fun(b)
  comb <- combination_forecaster(list(a = naive_per_day_forecast, b = naive_adj_forecast))(prep, target)
  expect_equal(comb$kwh, (naive_per_day_forecast(prep, target)$kwh +
                            naive_adj_forecast(prep, target)$kwh) / 2)
})

test_that("no leak with the robust preprocessing plugged into the backtest", {
  b <- sim_bills()
  targets <- as.Date(c("2020-08-01", "2021-02-01"))
  methods <- list(naive_adj = naive_adj_forecast, robust_level = robust_level_forecast)
  run <- function(bills) backtest(bills, targets, methods, as.Date("2021-06-01"), prep_fun)
  ref <- run(b)
  future <- b
  late <- future$month >= targets[2]  # after both origins
  future$kwh[late] <- future$kwh[late] * 7
  future$kwh_per_day[late] <- future$kwh_per_day[late] * 7
  out <- run(future)
  expect_equal(out$forecast, ref$forecast)
})

test_that("ward_classes splits two clear groups and picks k = 2", {
  set.seed(2)
  f <- data.frame(a = c(stats::rnorm(10, 0, 0.1), stats::rnorm(10, 5, 0.1)),
                  b = c(stats::rnorm(10, 0, 0.1), stats::rnorm(10, 5, 0.1)))
  cl <- ward_classes(f, c("a", "b"), 2:5)
  expect_equal(cl$k, 2)
  expect_equal(length(unique(cl$class[1:10])), 1)
  expect_false(cl$class[1] == cl$class[11])
})
