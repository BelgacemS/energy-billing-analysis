sim <- function(edit = NULL, n = 30, n_sites = 5, seed = 1) {
  set.seed(seed)
  month <- seq(as.Date("2019-01-01"), by = "month", length.out = n)
  days <- days_in_month(month)
  do.call(rbind, lapply(seq_len(n_sites), function(i) {
    x <- log(40 * i) + stats::rnorm(n, sd = 0.01)
    if (!is.null(edit) && i == 1) x <- edit(x)
    data.frame(site = sprintf("S%d", i), suffix = "KO", month = month, month_days = days,
               kwh = exp(x) * days, kwh_per_day = exp(x), stringsAsFactors = FALSE)
  }))
}
seg <- list(cost = "biweight", penalty = "MBIC", plateau_min = 3, K = 6)
prep_fun <- robust_preprocess(character(), seg, 0.10)

test_that("strategy_forecaster applies each site's method", {
  b <- sim()
  target <- add_months(max(b$month), 1)
  double <- function(prep, target) transform(naive_per_day_forecast(prep, target), kwh = 2 * kwh)
  sm <- c(S1 = "a", S2 = "b", S3 = "a", S4 = "b", S5 = "a")
  out <- strategy_forecaster(sm, list(a = naive_per_day_forecast, b = double))(b, target)
  ref <- naive_per_day_forecast(b, target)
  expect_equal(out$site, ref$site)
  expect_equal(out$kwh, unname(ref$kwh * ifelse(sm[out$site] == "b", 2, 1)))
  expect_equal(out$info, unname(sm[out$site]))
})

test_that("end_anomaly flags a last-month spike and a recent break only", {
  expect_false(any(end_anomaly(prep_fun(sim()), 3)))
  spike <- end_anomaly(prep_fun(sim(function(x) { x[30] <- x[30] + log(2); x })), 3)
  expect_equal(unname(spike), c(TRUE, FALSE, FALSE, FALSE, FALSE))
  recent <- end_anomaly(prep_fun(sim(function(x) { x[28:30] <- x[28:30] + log(2); x })), 3)
  expect_true(recent[["S1"]])
  old <- end_anomaly(prep_fun(sim(function(x) { x[15:30] <- x[15:30] + log(2); x })), 3)
  expect_false(old[["S1"]])
})

test_that("empirical quantiles, Gaussian bounds and coverage", {
  q <- empirical_quantiles(c(-2, -1, 0, 1, 2, 10), c(rep("a", 5), "b"), 0.5)
  expect_equal(q$lo[q$group == "a"], -1)
  expect_equal(q$hi[q$group == "a"], 1)
  expect_equal(q$n, c(5, 1))
  g <- gaussian_bounds(100, 0.1, 0.95)
  expect_equal(unname(log(g[, "hi"] / 100)), unname(-log(g[, "lo"] / 100)))
  expect_equal(unname(log(g[, "hi"] / 100)), stats::qnorm(0.975) * 0.1)
  expect_equal(coverage(c(1, 5, 10), c(0, 0, 0), c(2, 5, 9)), 2 / 3)
})

test_that("temporary_scenario goes back to the level before a recent break", {
  b <- sim(function(x) { x[28:30] <- x[28:30] + log(2); x })
  target <- add_months(max(b$month), 1)
  prep <- prep_fun(b)
  tmp <- temporary_scenario(prep, target, 3)
  lasting <- naive_per_day_forecast(prep, target)
  expect_lt(abs(log(tmp$kwh[1] / (40 * days_in_month(target)))), 0.05)  # old level
  expect_gt(lasting$kwh[1] / tmp$kwh[1], 1.8)
  spike <- prep_fun(sim(function(x) { x[30] <- x[30] + log(2); x }))
  expect_equal(temporary_scenario(spike, target, 3)$kwh, robust_level_forecast(spike, target)$kwh)
})

test_that("log_change_sd is the sd of monthly log changes", {
  b <- sim()
  x <- log(b$kwh_per_day[b$site == "S2"])
  expect_equal(log_change_sd(b)[["S2"]], stats::sd(diff(x)))
})
