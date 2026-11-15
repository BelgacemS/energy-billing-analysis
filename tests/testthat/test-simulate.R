cfg <- config
cfg$segment_K <- 3

flat_parts <- function(n_sites = 4, n = 35, sd = 0.01, break_at = NULL, seed = 1) {
  set.seed(seed)
  m <- seq(as.Date("2019-01-01"), by = "month", length.out = n)
  parts <- lapply(seq_len(n_sites), function(i) {
    level <- rep(log(100), n)
    if (!is.null(break_at)) level[break_at:n] <- log(150)
    list(site = sprintf("S%d", i), suffix = "KO", month = m, month_days = days_in_month(m),
         effect = rep(0, n), level = level, resid = stats::rnorm(n, sd = sd),
         break_index = if (is.null(break_at)) integer() else break_at)
  })
  setNames(parts, vapply(parts, `[[`, "", "site"))
}

test_that("biweight_costs matches biweight_fit on every segment", {
  set.seed(2)
  z <- c(stats::rnorm(12), stats::rnorm(10) + 5)
  z[4] <- 20
  C <- biweight_costs(z, 3)
  for (e in seq_along(z)) for (k in seq_len(e)) {
    expect_equal(C[k, e], biweight_fit(z[k:e], 3)$cost, tolerance = 1e-10)
  }
  expect_true(all(is.na(C[lower.tri(C)])))
})

test_that("the biweight penalty is MBIC (beta = 3) or a number beta", {
  z <- c(rep(0, 15), rep(4, 15)) + rep(c(-0.3, 0.3), 15)
  expect_equal(find_changepoints(z, penalty = "MBIC", minseglen = 3),
               find_changepoints(z, penalty = 3, minseglen = 3))
  expect_equal(find_changepoints(z, penalty = 1000, minseglen = 3), 30L)  # no break pays
  expect_error(find_changepoints(z, penalty = "BIC"), "MBIC")
})

test_that("block_bootstrap keeps runs of consecutive values", {
  set.seed(3)
  x <- block_bootstrap(1:30, 35, 6)
  expect_length(x, 35)
  expect_true(all(diff(x)[-(6 * (1:5))] == 1))  # consecutive inside each block
})

test_that("anomaly factors and injection", {
  expect_equal(anomaly_factors("spike", 0.5, -1, 3, 10), 1 / 1.5)
  expect_equal(anomaly_factors("pair", 0.2, 1, 3, 10), c(0.8, 0.8, 1.4))
  expect_equal(anomaly_factors("pair", 0.2, -1, 3, 10), c(1.2, 1.2, 0.6))
  expect_equal(anomaly_factors("level_break", 1, 1, 3, 4), rep(2, 4))
  b <- simulate_panel(flat_parts(2), 6)
  start <- c(S1 = 5L, S2 = 10L)
  inj <- inject_anomaly(b, "temporary", 0.1, 1, start, 3)
  ratio <- inj$kwh / b$kwh
  expect_equal(which(ratio != 1), c(5:7, 35 + 10:12))
  expect_equal(inj$kwh / inj$kwh_per_day, b$month_days)
  expect_equal(anomaly_start(c(0, 0.999), 35, 3, 3), c(4, 30))
})

test_that("simulate_panel is reproducible and centred on the level", {
  parts <- flat_parts(3)
  set.seed(9)
  a <- simulate_panel(parts, 6)
  set.seed(9)
  expect_identical(simulate_panel(parts, 6), a)
  expect_equal(stats::median(log(a$kwh_per_day)), log(100), tolerance = 0.01)
})

test_that("a large spike and a large break are found, the background is quiet", {
  parts <- flat_parts(6, sd = 0.01)
  set.seed(4)
  b <- simulate_panel(parts, 6)
  eff <- month_effect(b, character())
  start <- setNames(rep(15L, 6), names(parts))
  spike <- inject_anomaly(b, "spike", 0.5, 1, start, 3)
  sc <- score_injection(detect_panel(spike, eff, 3, 3, cfg), spike, "spike", 1, start, 3, 1, 0.1)
  expect_true(all(sc$detected))
  expect_true(all(sc$type == "point"))
  brk <- inject_anomaly(b, "level_break", 0.5, -1, start, 3)
  sc <- score_injection(detect_panel(brk, eff, 3, 3, cfg), brk, "level_break", -1, start, 3, 1, 0.1)
  expect_true(all(sc$detected))
  expect_false(any(sc$as_spike))
  # the wrong direction is not a detection
  sc_up <- score_injection(detect_panel(brk, eff, 3, 3, cfg), brk, "level_break", 1, start, 3, 1, 0.1)
  expect_false(any(sc_up$detected))
  expect_equal(false_alerts(detect_panel(b, eff, 6, 12, cfg), parts, 1), 0)
})

test_that("false_alerts does not count a true break found one month off", {
  parts <- flat_parts(1, break_at = 20)
  m <- parts$S1$month
  fl <- data.frame(site = "S1", month = m[c(21, 25)], suspect = TRUE, episode = c(NA, "S1-1"),
                   change = c(0.5, NA), type = c("lasting_change", "point"))
  expect_equal(false_alerts(fl, parts, 1), 1)  # only the point
  fl$change[1] <- -0.5                          # wrong direction
  expect_equal(false_alerts(fl, parts, 1), 2)
})

test_that("choose_setting follows the rule", {
  fa <- c(0.02, 0.008, 0.005, 0.009)
  score <- c(0.9, 0.70, 0.698, 0.6)
  expect_equal(choose_setting(fa, score, 0.01, 0.005), 3)  # tie, lower false alerts
  expect_equal(choose_setting(fa, score, 0.01, 0.001), 2)
  expect_equal(choose_setting(fa, score, 0.001, 0.005), 3)  # none admissible
})
