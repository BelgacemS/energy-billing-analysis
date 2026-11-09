test_that("site_noise is the MAD of first differences over sqrt(2)", {
  mat <- rbind(a = c(0, 0.1, 0, 0.1, 0), b = rep(1, 5))
  expect_equal(unname(site_noise(mat)), c(stats::mad(diff(mat[1, ])) / sqrt(2), 0))
})

test_that("stl_strength separates a seasonal series from noise", {
  set.seed(1)
  n <- 36
  season <- ts(0.1 * sin(2 * pi * (1:n) / 12) + rnorm(n, sd = 0.01), frequency = 12)
  noise <- ts(rnorm(n, sd = 0.01), frequency = 12)
  expect_gt(stl_strength(season)["season"], 0.9)
  expect_lt(stl_strength(noise)["season"], 0.5)
  expect_gt(mstl_season_strength(season), 0.9)
})

test_that("sd_mean_slope recovers a proportional spread", {
  m <- c(10, 20, 40, 80, 160)
  s <- sd_mean_slope(m, 0.05 * m * exp(c(0.01, -0.01, 0.02, -0.02, 0)))
  expect_equal(unname(s["slope"]), 1, tolerance = 0.02)
  expect_lt(s["lower"], 1)
  expect_gt(s["upper"], 1)
})

test_that("pick_typical_sites returns three distinct sites by the written rule", {
  months <- seq(as.Date("2019-01-01"), by = "month", length.out = 30)
  set.seed(2)
  base <- function(sd) 0.01 * sd * rnorm(30)
  mat <- rbind(
    quiet = base(0.2),
    noisy = base(2),
    step = c(rep(0, 15), rep(log(0.5), 15)) + base(1),
    spike = replace(base(1), 15, log(1.5))
  )
  colnames(mat) <- format(months, "%Y-%m")
  bills <- data.frame(site = rep(rownames(mat), each = 30), suffix = "BR",
                      month = rep(months, 4), month_days = 30)
  bills$kwh_per_day <- exp(as.vector(t(mat)))
  bills$kwh <- bills$kwh_per_day * 30
  p <- pick_typical_sites(bills, mat, window = 6, min_shift = 0.1, half_width = 3, threshold = 3)
  expect_equal(unname(p), c("quiet", "step", "spike"))
})
