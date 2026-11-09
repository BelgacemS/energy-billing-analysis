noise <- function(n, seed) {
  set.seed(seed)
  stats::rnorm(n)
}

test_that("noise_scale ignores a step and a spike", {
  x <- noise(200, 1)
  base <- noise_scale(x)
  x[101:200] <- x[101:200] + 10
  x[50] <- x[50] + 20
  expect_equal(noise_scale(x), base, tolerance = 0.1)
})

test_that("find_changepoints returns the last month of each segment", {
  x <- c(rep(0, 20), rep(5, 15)) + rep(c(-0.1, 0.1), length.out = 35)
  expect_equal(find_changepoints(x / 0.1, "normal"), c(20L, 35L))
  expect_equal(find_changepoints(x, "np"), c(20L, 35L))
})

test_that("a clear step is found at the right month", {
  y <- noise(35, 2) * 0.03
  y[21:35] <- y[21:35] + 0.2
  fit <- segment_series(y)
  expect_equal(fit$breaks$index, 21L)
  expect_false(any(fit$months$excursion))
  expect_equal(fit$breaks$level_after - fit$breaks$level_before, 0.2, tolerance = 0.05)
})

test_that("an isolated spike is an excursion, not a lasting break", {
  y <- noise(35, 3) * 0.03
  y[15] <- y[15] + 0.3
  fit <- segment_series(y)
  expect_equal(nrow(fit$breaks), 0)
  expect_equal(which(fit$months$excursion), 15L)
  # The spike keeps the level of its plateau, so its residual is the spike.
  expect_equal(fit$months$level[15], median(y[-15]))
})

test_that("a constant noisy series has no break", {
  for (seed in 4:8) {
    fit <- segment_series(noise(35, seed) * 0.03)
    expect_equal(nrow(fit$breaks), 0)
  }
})

test_that("excursions take the previous plateau, or the next one at the start", {
  y <- noise(35, 9) * 0.03
  y[1] <- y[1] - 0.5          # excursion before the first plateau
  y[21:35] <- y[21:35] + 0.2  # step at 21
  y[34:35] <- y[34:35] + 0.5  # two-month excursion at the end
  fit <- segment_series(y)
  m <- fit$months
  expect_equal(which(m$excursion), c(1L, 34L, 35L))
  expect_equal(fit$breaks$index, 21L)
  expect_equal(m$plateau, c(rep(1L, 20), rep(2L, 15)))
  expect_equal(m$level[35], median(y[21:33]))
})

test_that("a noisy final regime is a plateau, not a run of excursions", {
  y <- noise(35, 10) * 0.02
  set.seed(11)
  y[29:35] <- -3 + stats::rnorm(7, sd = 0.2)  # drop with ten times the noise
  fit <- segment_series(y)
  expect_false(any(fit$months$excursion[29:35]))
  expect_equal(fit$breaks$index[1], 29L)
})

test_that("a cluster of spikes around the level creates no break", {
  y <- noise(35, 12) * 0.03
  y[10] <- y[10] - 1.5   # under-billing,
  y[12] <- y[12] + 1     # catch-up two months later, as 0089KO in 2019
  fit <- segment_series(y)
  expect_equal(nrow(fit$breaks), 0)
  expect_lt(abs(fit$months$level[10] - median(y)), 0.02)
})

test_that("segment_series refuses a series without noise", {
  expect_error(segment_series(rep(1, 35)), "Zero noise")
})
