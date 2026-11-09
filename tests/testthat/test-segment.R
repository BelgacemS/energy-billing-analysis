noise <- function(n, seed) {
  set.seed(seed)
  stats::rnorm(n)
}

# All partitions of 1..n into segments of at least m points, as end indices.
all_partitions <- function(n, m) {
  if (n < m) return(list())
  out <- if (n < 2 * m) list() else
    unlist(lapply(m:(n - m), function(e) lapply(all_partitions(e, m), function(p) c(p, n))),
           recursive = FALSE)
  c(list(n), out)
}

test_that("noise_scale ignores a step and a spike", {
  x <- noise(200, 1)
  base <- noise_scale(x)
  x[101:200] <- x[101:200] + 10
  x[50] <- x[50] + 20
  expect_equal(noise_scale(x), base, tolerance = 0.1)
})

test_that("biweight_fit matches a fine grid search", {
  grid <- seq(-10, 20, by = 1e-3)
  for (seed in 1:5) {
    z <- c(noise(8, seed), 12, -6)
    brute <- vapply(grid, function(m) sum(pmin((z - m)^2, 9)), numeric(1))
    fit <- biweight_fit(z, K = 3)
    expect_equal(fit$cost, min(brute), tolerance = 1e-5)
    expect_equal(fit$level, mean(z[abs(z - fit$level) < 3]))
  }
  expect_equal(biweight_fit(c(0, 1, 50), K = 3)$level, 0.5)
})

test_that("biweight_partition matches an exhaustive search", {
  K <- 3
  pen <- 4
  for (seed in 1:4) {
    z <- noise(9, seed) + c(rep(0, 4), rep(5, 5))
    z[2] <- z[2] + 10
    for (m in 1:3) {
      parts <- all_partitions(length(z), m)
      total <- vapply(parts, function(e) {
        s <- c(1, head(e, -1) + 1)
        sum(mapply(function(a, b) biweight_fit(z[a:b], K)$cost, s, e)) + pen * length(e)
      }, numeric(1))
      expect_equal(biweight_partition(z, K, pen, m), parts[[which.min(total)]])
    }
  }
})

test_that("cpt.mean and cpt.np return the last month of each segment", {
  x <- c(rep(0, 20), rep(5, 15)) + rep(c(-0.1, 0.1), length.out = 35)
  expect_equal(find_changepoints(x / 0.1, "normal"), c(20L, 35L))
  expect_equal(find_changepoints(x, "np"), c(20L, 35L))
  expect_equal(find_changepoints(x, "normal", minseglen = 20), 35L)  # too short
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
})

# Probabilistic properties are tested as rates over 100 series: a fixed seed
# can land in the few percent of false alarms the method is expected to give.
break_rate <- function(make, reps = 100) {
  mean(vapply(seq_len(reps), function(r) {
    set.seed(r)
    nrow(segment_series(make())$breaks) > 0
  }, logical(1)))
}

test_that("a constant noisy series rarely gets a break", {
  expect_lt(break_rate(function() stats::rnorm(35) * 0.03), 0.10)  # about 4 %
})

test_that("an under-billing over two months and its catch-up create no plateau", {
  pair <- function() {
    y <- stats::rnorm(35) * 0.03
    y[10:12] <- y[10:12] + c(-0.5, -0.5, 1)
    y
  }
  expect_lt(break_rate(pair), 0.10)
  y <- noise(35, 9) * 0.03
  y[10:12] <- y[10:12] + c(-0.5, -0.5, 1)
  fit <- segment_series(y)
  expect_equal(which(fit$months$excursion), 10:12)
  # The gaps are measured from the usual level, not from a false plateau.
  expect_lt(abs(fit$months$level[11] - median(y[-(10:12)])), 0.02)
})

test_that("normal months before a ramp keep the old level (1911BB in 2019)", {
  y <- noise(35, 14) * 0.03
  y[11:12] <- y[11:12] + c(0.4, 0.9)  # ramp
  y[13:35] <- y[13:35] + 1.1          # new level
  fit <- segment_series(y)
  expect_equal(nrow(fit$breaks), 1)
  expect_true(fit$breaks$index %in% 11:13)
  expect_false(any(fit$months$excursion[1:10]))
  expect_true(all(fit$months$excursion[11:12]))
})

test_that("a noisy final regime is found at its first month", {
  y <- noise(35, 10) * 0.02
  set.seed(11)
  y[29:35] <- -3 + stats::rnorm(7, sd = 0.2)  # drop with ten times the noise
  fit <- segment_series(y)
  # The regime may be cut in several plateaus, all far below the old level.
  expect_equal(fit$breaks$index[1], 29L)
  expect_true(all(fit$months$level[29:35] < -2.5))
})

test_that("a change in the last two months is not a plateau yet", {
  y <- noise(35, 15) * 0.03
  y[34:35] <- y[34:35] + 0.5
  fit <- segment_series(y)
  expect_equal(nrow(fit$breaks), 0)
  expect_equal(which(fit$months$excursion), 34:35)
})

test_that("segment_series refuses a series without noise", {
  expect_error(segment_series(rep(1, 35)), "Zero noise")
})
