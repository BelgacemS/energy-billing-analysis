test_that("hampel flags an isolated spike and nothing on a smooth series", {
  x <- sin(seq_len(30) / 3) * 0.01
  h <- hampel(x, half_width = 3, threshold = 3)
  expect_false(any(h$flag))
  x[15] <- x[15] + 1
  h <- hampel(x, half_width = 3, threshold = 3)
  expect_equal(which(h$flag), 15L)
  expect_equal(h$median[15], median(x[12:18]))
})

test_that("hampel truncates the window at the ends", {
  x <- c(10, 1:9)
  h <- hampel(x, half_width = 2, threshold = 3)
  expect_equal(h$median[1], median(x[1:3]))
  expect_true(h$flag[1])
})

test_that("min_scale stops a flat window from flagging small changes", {
  x <- c(rep(0, 10), 0.001, rep(0, 10))
  expect_true(hampel(x, 3, 3)$flag[11])  # MAD is 0
  expect_false(hampel(x, 3, 3, min_scale = 0.01)$flag[11])
})

test_that("hampel_bills gives the gap in kWh with its sign, site by site", {
  months <- seq(as.Date("2020-01-01"), by = "month", length.out = 12)
  days <- days_in_month(months)
  rate <- c(10 * (1 + 0.01 * sin(1:12)), 20 * (1 + 0.01 * cos(1:12)))
  b <- data.frame(site = rep(c("B", "A"), each = 12), suffix = "BR",
                  month = rep(months, 2), month_days = rep(days, 2))
  b$kwh_per_day <- rate
  b$kwh_per_day[18] <- b$kwh_per_day[18] * 0.5  # A, June: half the usual bill
  b$kwh <- b$kwh_per_day * b$month_days
  h <- hampel_bills(b, half_width = 3, threshold = 3)
  expect_equal(h$site, rep(c("A", "B"), each = 12))
  flagged <- h[h$flag, ]
  expect_equal(nrow(flagged), 1)
  expect_equal(flagged$site, "A")
  expect_equal(flagged$month, as.Date("2020-06-01"))
  expect_lt(flagged$gap_kwh, 0)
  expect_equal(flagged$gap_kwh, flagged$kwh - flagged$expected_kwh)
})

test_that("hampel marks truncated windows at both ends", {
  h <- hampel(1:10 / 100, half_width = 3, threshold = 3)
  expect_equal(h$edge, c(TRUE, TRUE, TRUE, rep(FALSE, 4), TRUE, TRUE, TRUE))
})
