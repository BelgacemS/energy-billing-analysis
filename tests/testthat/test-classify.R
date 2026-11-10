# Hand-made panels: level and deviation (in noise units) of every month are
# given, so the expected type of each flagged bill is known.
n_months <- 35
base <- log(100)

make_site <- function(site, level = rep(base, n_months), dev = numeric(n_months),
                      sigma = 0.01) {
  month <- seq(as.Date("2019-01-01"), by = "month", length.out = length(level))
  days <- days_in_month(month)
  y <- level + dev * sigma
  data.frame(site = site, suffix = "BR", group = "common", month = month,
             kwh = exp(y) * days, month_days = days, x = y, effect = 0, y = y,
             plateau = cumsum(c(1, diff(level) != 0)), level = level,
             excursion = abs(dev) > config$segment_K, sigma = sigma, expected_log = level,
             expected_kwh = exp(level) * days, gap_kwh = exp(y) * days - exp(level) * days,
             resid = dev * sigma, z = dev, stringsAsFactors = FALSE)
}
make_breaks <- function(p) {
  do.call(rbind, c(list(data.frame(site = character(), month = as.Date(character()),
                                   level_before = numeric(), level_after = numeric(),
                                   change = numeric())),
                   lapply(split(p, p$site), function(d) {
    i <- which(diff(d$level) != 0) + 1
    data.frame(site = d$site[i], month = d$month[i], level_before = d$level[i - 1],
               level_after = d$level[i], change = exp(d$level[i] - d$level[i - 1]) - 1,
               stringsAsFactors = FALSE)
  })))
}
no_periods <- data.frame(site = character(), month = as.Date(character()))
classify <- function(p, periods = no_periods, big = character()) {
  classify_bills(p, make_breaks(p), periods, big, config)
}
step_level <- function(at, change) c(rep(base, at - 1), rep(base + log(1 + change), n_months - at + 1))

test_that("find_episodes merges runs separated by at most one normal month", {
  ep <- find_episodes(c(FALSE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, TRUE), gap = 1)
  expect_equal(ep$start, c(2, 5, 8))
  expect_equal(ep$end, c(3, 5, 8))
  expect_equal(ep$episode, c(1, 1, 2))
  expect_equal(nrow(find_episodes(rep(FALSE, 5), 1)), 0)
})

test_that("a quiet series has no flagged bill", {
  expect_equal(nrow(classify(make_site("A"))), 0)
})

test_that("an irregular period is flagged even on its plateau", {
  p <- make_site("A")
  out <- classify(p, periods = data.frame(site = "A", month = p$month[10]))
  expect_equal(out$type, "irregular_period")
  expect_equal(out$month, p$month[10])
  expect_true(out$suspect)
  expect_equal(as.character(out$confidence), "high")  # deterministic check
  out <- classify(p, periods = data.frame(site = "A", month = p$month[34]))
  expect_equal(as.character(out$confidence), "high")  # even in the edge months
})

test_that("a point gap, its kWh, direction and range (hand computed)", {
  dev <- numeric(n_months)
  dev[12] <- 8
  dev[20] <- -5
  p <- make_site("A", dev = dev)
  out <- classify(p)
  expect_equal(out$type, c("point", "point"))
  expect_equal(out$episode, c("A-1", "A-2"))
  pt <- out[out$month == p$month[12], ]
  expected <- 100 * 31  # December 2019, 100 kWh a day
  expect_equal(pt$expected_kwh, expected)
  expect_equal(pt$kwh, expected * exp(0.08))
  expect_equal(pt$gap_kwh, expected * (exp(0.08) - 1))
  expect_equal(pt$gap_low, expected * (exp(0.08) - exp(0.02)))
  expect_equal(pt$gap_high, expected * (exp(0.08) - exp(-0.02)))
  expect_equal(pt$direction, "over")
  expect_equal(as.character(pt$confidence), "high")    # |z| = 8 >= 6
  low <- out[out$month == p$month[20], ]
  expect_equal(low$direction, "under")
  expect_equal(as.character(low$confidence), "medium")  # |z| = 5 < 6
})

test_that("a bill more than 3 times the expected value is implausible", {
  dev <- numeric(n_months)
  dev[12] <- log(4) / 0.01
  out <- classify(make_site("A", dev = dev))
  expect_equal(out$type, "implausible")
})

test_that("a few months of the same sign then back are a temporary gap", {
  dev <- numeric(n_months)
  dev[10:12] <- 5
  out <- classify(make_site("A", dev = dev))
  expect_equal(out$type, rep("temporary", 3))
  expect_equal(unique(out$episode), "A-1")
})

test_that("under-billing then catch-up summing to the level is a compensated pair", {
  dev <- numeric(n_months)
  dev[15] <- log(0.5) / 0.01  # half the usual bill
  dev[17] <- log(1.5) / 0.01  # catch-up after a normal month: 0.5 + 1 + 1.5 = 3
  out <- classify(make_site("A", dev = dev))
  expect_equal(out$type, rep("compensated_pair", 2))
  expect_equal(unique(out$episode), "A-1")
  expect_true(all(out$suspect))
})

test_that("a ramp between two levels is a transition, the break a lasting change", {
  level <- step_level(21, 1)  # doubles from month 21
  dev <- numeric(n_months)
  dev[19:20] <- log(c(1.3, 1.6)) / 0.01
  p <- make_site("A", level, dev)
  out <- classify(p)
  ramp <- out[out$type == "transition", ]
  expect_equal(ramp$month, p$month[19:20])
  expect_false(any(ramp$suspect))
  brk <- out[out$type == "lasting_change", ]
  expect_equal(brk$month, p$month[21])
  # the cumulative gap starts with the ramp, at the old level
  expect_equal(brk$cumul_from, p$month[19])
  d <- p$month_days
  expect_equal(brk$cumul_kwh, 100 * (0.3 * d[19] + 0.6 * d[20] + sum(d[21:n_months])))
  expect_true(brk$suspect)
})

test_that("a lasting change carries its cumulative gap since the break", {
  p <- make_site("A", step_level(21, 0.5))
  out <- classify(p)
  expect_equal(out$type, "lasting_change")
  expect_equal(out$change, 0.5)
  days <- p$month_days
  expect_equal(out$expected_kwh, 100 * days[21])
  expect_equal(out$gap_kwh, 50 * days[21])
  expect_equal(out$cumul_kwh, 50 * sum(days[21:n_months]))
  expect_equal(out$direction, "over")
})

test_that("a break shared by three sites is a probable network change", {
  p <- rbind(make_site("A", step_level(21, 0.5)), make_site("B", step_level(21, 0.3)),
             make_site("C", step_level(21, 0.4)), make_site("D", step_level(21, -0.4)))
  out <- classify(p)
  expect_equal(out$type[out$site != "D"], rep("network_change", 3))
  expect_false(any(out$suspect[out$site != "D"]))
  expect_true(all(out$confidence[out$site != "D"] == "medium"))  # never high
  expect_equal(out$type[out$site == "D"], "lasting_change")       # opposite direction
})

test_that("breaks of 10 % or less are ignored", {
  dev <- numeric(n_months)
  dev[15] <- 8
  out <- classify(make_site("A", step_level(16, 0.05), dev))
  expect_equal(out$type, "point")  # back to the same level despite the small break
})

test_that("an episode or a break near the end is indeterminate", {
  dev <- numeric(n_months)
  dev[34] <- 8
  out <- classify(make_site("A", dev = dev))
  expect_equal(out$type, "indeterminate")
  expect_equal(as.character(out$confidence), "low")
  out <- classify(make_site("A", step_level(33, 0.5)))
  expect_equal(out$type, "indeterminate")
})

test_that("a gap at a break that is no ramp is indeterminate", {
  dev <- numeric(n_months)
  dev[21] <- -100  # below the old level, outside the ramp band
  out <- classify(make_site("A", step_level(21, 1), dev))
  expect_equal(out$type[!is.na(out$episode)], "indeterminate")
})

test_that("edge months get a low confidence and big sites are marked", {
  dev <- numeric(n_months)
  dev[2] <- 10
  out <- classify(make_site("Z", dev = dev), big = "Z")
  expect_equal(as.character(out$confidence), "low")
  expect_true(out$big)
})

test_that("rank_flagged puts suspects first, then confidence, then |gap|", {
  tab <- data.frame(suspect = c(FALSE, TRUE, TRUE, TRUE),
                    confidence = factor(c("high", "low", "high", "high"),
                                        levels = c("high", "medium", "low")),
                    gap_kwh = c(1e6, 1e5, -10, 500))
  expect_equal(rownames(rank_flagged(tab)), c("4", "3", "2", "1"))
})

test_that("0089KO 2019: under-billing in June and catch-up in August form a pair", {
  bills <- load_bills(root_path(config$data_path))
  big <- big_sites(bills, config$big_site_kwh)
  eff <- month_effect(bills, big)
  fit <- panel_fit(bills[bills$site == "0089KO", ], eff, big,
                   list(cost = config$segment_cost, penalty = config$segment_penalty,
                        plateau_min = config$plateau_min_months, K = config$segment_K))
  out <- classify_bills(fit$bills, fit$breaks, check_periods(bills), big, config)
  pair <- out[out$type == "compensated_pair", ]
  expect_equal(pair$month, as.Date(c("2019-06-01", "2019-08-01")))
  expect_equal(pair$direction, c("under", "over"))
  expect_equal(unique(pair$pair_order), "under_then_over")
  expect_equal(unique(pair$pair_net_kwh), -519.4, tolerance = 0.1, scale = 1)
  expect_equal(out$type[out$month == as.Date("2020-02-01")], "irregular_period")
})

pair_dev <- function(months, ratios) {
  dev <- numeric(n_months)
  dev[months] <- log(ratios) / 0.01
  dev
}

test_that("a pair may start with the bill too high, with its order and net gap", {
  p <- make_site("A", dev = pair_dev(15:16, c(1.5, 0.5)))  # March 2020 (31 d), April (30 d)
  out <- classify(p)
  expect_equal(out$type, rep("compensated_pair", 2))
  expect_equal(unique(out$pair_order), "over_then_under")
  expect_equal(unique(out$pair_net_kwh), 100 * (1.5 * 31 + 0.5 * 30) - 100 * (31 + 30))
})

test_that("a pair may have three bills or more", {
  out <- classify(make_site("A", dev = pair_dev(10:12, c(0.5, 0.5, 2))))
  expect_equal(out$type, rep("compensated_pair", 3))
  expect_equal(unique(out$pair_order), "under_then_over")
  expect_equal(unique(out$pair_net_kwh), 100 * (0.5 * 31 + 0.5 * 30 + 2 * 31) - 100 * 92)
})

test_that("the pair tolerance is a sharp limit", {
  # months 15 and 17 (31 days) around a normal April (30 days): sum / expected - 1
  # = (31 x - 46.5) / 92, so x = 1.767 gives +9 % and x = 1.826 gives +11 %
  inside <- classify(make_site("A", dev = pair_dev(c(15, 17), c(0.5, 1.767))))
  expect_equal(inside$type, rep("compensated_pair", 2))
  outside <- classify(make_site("A", dev = pair_dev(c(15, 17), c(0.5, 1.826))))
  expect_equal(outside$type, c("point", "point"))
  expect_true(all(is.na(outside$pair_net_kwh)))
})

test_that("mixed signs that do not compensate are indeterminate", {
  out <- classify(make_site("A", dev = pair_dev(10:11, c(2, 0.9))))  # +46 % over the two months
  expect_equal(out$type, rep("indeterminate", 2))
})

test_that("longest_temporary measures the guard of temporary_max_months", {
  out <- classify(make_site("A", dev = pair_dev(10:14, rep(1.1, 5))))
  expect_equal(out$type, rep("temporary", 5))  # longer runs that come back stay temporary
  expect_equal(longest_temporary(out), 5L)
  expect_gt(longest_temporary(out), config$temporary_max_months)
  expect_equal(longest_temporary(out[0, ]), 0L)
})

test_that("testable months leave out the first and last plateau_min - 1 months", {
  months <- seq(as.Date("2019-01-01"), by = "month", length.out = 35)
  t <- testable_months(months, 3)
  expect_equal(length(t), 30)
  expect_equal(range(t), months[c(4, 33)])
})

test_that("a shared month must beat chance, not only reach 3 sites", {
  testable <- seq(as.Date("2019-04-01"), by = "month", length.out = 30)
  # 27 rises: 8 in month 21, 3 in month 10, one in each of 16 other months
  # (lambda = 0.9 a month); 1 fall
  singles <- testable[setdiff(seq_len(30), c(7, 18))][1:16]
  brk <- data.frame(site = sprintf("S%02d", 1:28),
                    month = c(rep(testable[18], 8), rep(testable[7], 3), singles, testable[18]),
                    change = c(rep(0.5, 27), -0.5))
  out <- shared_break_months(brk, big = character(), testable, 0.10, 3, 0.05)
  up <- out[out$direction == "up", ]
  expect_equal(unique(up$lambda), 27 / 30)
  m21 <- up[up$month == testable[18], ]
  expect_equal(m21$p_value, stats::ppois(7, 0.9, lower.tail = FALSE))
  expect_true(m21$shared)
  m10 <- up[up$month == testable[7], ]
  expect_equal(m10$n, 3)
  expect_gt(m10$p_value, 0.05)  # 3 sites at lambda = 0.9: P(X >= 3) = 0.063
  expect_false(m10$shared)
  expect_false(out$shared[out$direction == "down"])
  expect_false(any(up$shared[up$n == 1]))
})
