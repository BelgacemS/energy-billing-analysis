# Typology of the flagged bills (step 6). A bill is flagged when its period
# is irregular, when it lies off its plateau (excursion) or when it opens a
# new level; each flagged bill gets a type, a gap in kWh with its range and a
# confidence level. Rules and thresholds are in R/config.R.

#' Months where a break can be dated: every plateau lasts at least
#' `plateau_min` months, so a new plateau starts between month plateau_min + 1
#' and month n - plateau_min + 1 of the series.
testable_months <- function(months, plateau_min) {
  months <- sort(unique(months))
  months[(plateau_min + 1):(length(months) - plateau_min + 1)]
}

#' Breaks larger than `min_change` counted by month and direction over the
#' ordinary sites, against chance.
#'
#' Same logic as check_repeats: with lambda = mean number of breaks in that
#' direction per testable month, p_value = P(X >= n) for X ~ Poisson(lambda).
#' A month is shared in a direction when n >= min_sites and p_value < alpha.
#' A site has at most one break a month, so n counts sites. Returns one row
#' per month and direction with at least one break.
shared_break_months <- function(breaks, big, testable, min_change, min_sites, alpha) {
  b <- breaks[abs(breaks$change) > min_change & !breaks$site %in% big, ]
  if (nrow(b) == 0) {
    return(data.frame(month = as.Date(character()), direction = character(), n = integer(),
                      lambda = numeric(), p_value = numeric(), shared = logical()))
  }
  b$direction <- ifelse(b$change > 0, "up", "down")
  lambda <- table(factor(b$direction, c("up", "down"))) / length(testable)
  b$n <- stats::ave(seq_len(nrow(b)), b$month, b$direction, FUN = length)
  out <- unique(b[c("month", "direction", "n")])
  out$lambda <- as.vector(lambda[out$direction])
  out$p_value <- stats::ppois(out$n - 1, out$lambda, lower.tail = FALSE)
  out$shared <- out$n >= min_sites & out$p_value < alpha
  out <- out[order(out$month, out$direction), ]
  rownames(out) <- NULL
  out
}

#' Runs of TRUE in a logical vector: first and last index of each run, and
#' its episode once runs separated by at most `gap` FALSE are merged.
find_episodes <- function(flag, gap) {
  r <- rle(flag)
  len <- r$lengths[r$values]
  end <- cumsum(r$lengths)[r$values]
  if (length(end) == 0) return(data.frame(start = integer(), end = integer(), episode = integer()))
  start <- end - len + 1L
  episode <- cumsum(c(TRUE, start[-1] - end[-length(end)] - 1L > gap))
  data.frame(start = start, end = end, episode = episode)
}

#' Type of every excursion bill of one site.
#'
#' `d`: bills of the site sorted by month (panel_fit columns); `flag`:
#' excursion months to classify (irregular periods left out); `brk`: breaks of
#' the site above the minimal change, with `index` = row of their first month.
#' Months i..j are "back to the same level" when the months just before and
#' just after them lie between the same two large breaks (smaller breaks are
#' ignored until step 7). Returns one row per bill: index, type, episode, and
#' for a compensated pair its order (the first bill below or above) and its
#' net gap (sum billed - sum expected over the episode, at the level that
#' validated the pair).
classify_episodes <- function(d, flag, brk, cfg) {
  n <- nrow(d)
  K <- cfg$segment_K
  sigma <- d$sigma[1]
  regime <- cumsum(seq_len(n) %in% brk$index)

  # Ramp: touches a large break (from its first month to one month after its
  # end) and every month lies between the two levels, K noise units of slack.
  is_transition <- function(s, e) {
    near <- brk[brk$index >= s & brk$index <= e + 1, ]
    any(vapply(seq_len(nrow(near)), function(k) {
      lv <- c(near$level_before[k], near$level_after[k])
      all(d$y[s:e] >= min(lv) - K * sigma & d$y[s:e] <= max(lv) + K * sigma)
    }, logical(1)))
  }
  # Bills below and above whose sum matches the level before or after (the
  # closer of the two is kept); NULL when there is no such pair.
  pair_check <- function(s, e) {
    z <- d$z[s:e][flag[s:e]]
    if (!(any(z < 0) && any(z > 0))) return(NULL)
    levels <- c(if (s > 1) d$level[s - 1], if (e < n) d$level[e + 1])
    expected <- vapply(levels, function(l) sum(exp(l + d$effect[s:e]) * d$month_days[s:e]),
                       numeric(1))
    miss <- abs(sum(d$kwh[s:e]) / expected - 1)
    if (!any(miss <= cfg$pair_tolerance)) return(NULL)
    list(order = if (z[1] < 0) "under_then_over" else "over_then_under",
         net = sum(d$kwh[s:e]) - expected[which.min(miss)])
  }
  comes_back <- function(s, e) e < n && regime[max(s - 1, 1)] == regime[e + 1]
  classify_run <- function(s, e) {
    if (is_transition(s, e)) return("transition")
    if (!comes_back(s, e)) return("indeterminate")
    if (s == e) {
      ratio <- d$kwh[s] / d$expected_kwh[s]
      r <- cfg$implausible_ratio
      return(if (ratio > r || ratio < 1 / r) "implausible" else "point")
    }
    # Runs longer than temporary_max_months that come back are temporary too.
    if (length(unique(sign(d$z[s:e]))) == 1) "temporary" else "indeterminate"
  }

  runs <- find_episodes(flag, cfg$merge_gap_months)
  rows <- lapply(split(runs, runs$episode), function(g) {
    s <- min(g$start)
    e <- max(g$end)
    at_end <- e > n - cfg$edge_months
    ramp <- !at_end && is_transition(s, e)
    pair <- if (!at_end && !ramp) pair_check(s, e)
    type <- if (at_end) "indeterminate"
      else if (ramp) "transition"
      else if (!is.null(pair)) "compensated_pair"
      else rep(mapply(classify_run, g$start, g$end), g$end - g$start + 1)
    data.frame(index = unlist(Map(seq, g$start, g$end)), type = type, episode = g$episode[1],
               pair_order = if (is.null(pair)) NA_character_ else pair$order,
               pair_net_kwh = if (is.null(pair)) NA_real_ else pair$net,
               stringsAsFactors = FALSE)
  })
  do.call(rbind, c(list(data.frame(index = integer(), type = character(), episode = integer(),
                                   pair_order = character(), pair_net_kwh = numeric())),
                   rows))
}

#' One row per large break of one site, at its first month.
#'
#' The gap and z-score are taken against the level before the break
#' (expected = old level + month effect). The cumulative gap sums billed -
#' expected at the old level from the start of the change to the end of the
#' series: the change starts at the first bill of a transition episode that
#' touches the break (the ramp is part of the change), at the break otherwise
#' (`cumul_from`). `episodes`: episode rows of the site (classify_episodes).
#' Type: lasting change, probable network change if shared, indeterminate in
#' the last `edge_months` months.
break_rows <- function(d, brk, episodes, cfg) {
  n <- nrow(d)
  ramps <- episodes[episodes$type == "transition", ]
  ramp_start <- tapply(ramps$index, ramps$episode, min)
  ramp_end <- tapply(ramps$index, ramps$episode, max)
  do.call(rbind, c(list(data.frame(index = integer(), type = character(), expected_kwh = numeric(),
                                   z = numeric(), shared = logical(), cumul_kwh = numeric(),
                                   cumul_from = as.Date(character()), change = numeric())),
                   lapply(seq_len(nrow(brk)), function(k) {
    b <- brk$index[k]
    from <- min(c(b, ramp_start[ramp_start <= b & ramp_end >= b - 1]))
    old <- exp(brk$level_before[k] + d$effect) * d$month_days
    type <- if (b > n - cfg$edge_months) "indeterminate"
      else if (brk$shared[k]) "network_change" else "lasting_change"
    data.frame(index = b, type = type, expected_kwh = old[b],
               z = (d$y[b] - brk$level_before[k]) / d$sigma[1], shared = brk$shared[k],
               cumul_kwh = sum(d$kwh[from:n] - old[from:n]), cumul_from = d$month[from],
               change = brk$change[k], stringsAsFactors = FALSE)
  })))
}

#' Confidence of a flagged bill: high when `certain` (deterministic check),
#' else low in the edge months, high when |z| >= high_z outside a shared
#' break, medium otherwise.
confidence_level <- function(z, edge, shared, high_z, certain = FALSE) {
  out <- ifelse(edge, "low", ifelse(abs(z) >= high_z & !shared, "high", "medium"))
  out[certain] <- "high"
  factor(out, levels = c("high", "medium", "low"))
}

#' Table of the flagged bills of the panel.
#'
#' `bills`: panel_fit()$bills; `breaks`: panel_fit()$breaks; `periods`:
#' check_periods() output; `big`: names of the big sites; `cfg`: the config
#' list. One row per flagged bill, plus one row per large break at its first
#' month (a bill of an episode can therefore also carry a break row). Gap =
#' billed - expected, range = billed - expected x exp(+-range_sigmas x noise).
#' Shared breaks: shared_break_months over the testable months.
classify_bills <- function(bills, breaks, periods, big, cfg) {
  months <- shared_break_months(breaks, big, testable_months(bills$month, cfg$plateau_min_months),
                                cfg$level_shift_min, cfg$shared_break_min_sites,
                                cfg$shared_break_alpha)
  brk_all <- breaks[abs(breaks$change) > cfg$level_shift_min, ]
  brk_all$shared <- paste(brk_all$month, ifelse(brk_all$change > 0, "up", "down")) %in%
    paste(months$month, months$direction)[months$shared]
  bills <- bills[order(bills$site, bills$month), ]
  irregular_key <- paste(periods$site, periods$month)
  out <- do.call(rbind, lapply(split(bills, bills$site), function(d) {
    n <- nrow(d)
    brk <- brk_all[brk_all$site == d$site[1], ]
    brk$index <- match(brk$month, d$month)
    irregular <- paste(d$site, d$month) %in% irregular_key
    k <- sum(irregular)
    episodes <- classify_episodes(d, d$excursion & !irregular, brk, cfg)
    rows <- rbind(data.frame(index = which(irregular), type = rep("irregular_period", k),
                             episode = rep(NA_integer_, k), pair_order = rep(NA_character_, k),
                             pair_net_kwh = rep(NA_real_, k), stringsAsFactors = FALSE),
                  episodes)
    m <- nrow(rows)
    rows <- cbind(rows, expected_kwh = d$expected_kwh[rows$index], z = d$z[rows$index],
                  shared = rep(FALSE, m), cumul_kwh = rep(NA_real_, m),
                  cumul_from = rep(as.Date(NA), m), change = rep(NA_real_, m))
    brows <- break_rows(d, brk, episodes, cfg)
    brows[c("episode", "pair_order", "pair_net_kwh")] <-
      list(rep(NA_integer_, nrow(brows)), rep(NA_character_, nrow(brows)), rep(NA_real_, nrow(brows)))
    rows <- rbind(rows, brows)  # rbind matches the columns by name
    i <- rows$index
    data.frame(site = rep(d$site[1], nrow(rows)), suffix = d$suffix[i], month = d$month[i],
               kwh = d$kwh[i], sigma = d$sigma[i],
               edge = i <= cfg$edge_months | i > n - cfg$edge_months, rows[-1],
               stringsAsFactors = FALSE)
  }))
  rownames(out) <- NULL
  r <- cfg$range_sigmas
  data.frame(
    site = out$site, suffix = out$suffix, big = out$site %in% big, month = out$month,
    kwh = out$kwh, expected_kwh = out$expected_kwh, gap_kwh = out$kwh - out$expected_kwh,
    direction = ifelse(out$kwh > out$expected_kwh, "over", "under"),
    gap_low = out$kwh - out$expected_kwh * exp(r * out$sigma),
    gap_high = out$kwh - out$expected_kwh * exp(-r * out$sigma),
    z = out$z, type = out$type,
    episode = ifelse(is.na(out$episode), NA_character_, paste0(out$site, "-", out$episode)),
    edge = out$edge,
    confidence = confidence_level(out$z, out$edge, out$shared, cfg$high_confidence_z,
                                  certain = out$type == "irregular_period"),
    suspect = !out$type %in% c("transition", "network_change"),
    pair_order = out$pair_order, pair_net_kwh = out$pair_net_kwh,
    cumul_kwh = out$cumul_kwh, cumul_from = out$cumul_from, change = out$change,
    stringsAsFactors = FALSE
  )
}

#' Suspect bills first, then by confidence, then by decreasing |gap in kWh|.
rank_flagged <- function(tab) {
  tab[order(!tab$suspect, tab$confidence, -abs(tab$gap_kwh)), ]
}

#' Largest number of bills in one temporary episode (0 if none): the guard
#' of temporary_max_months.
longest_temporary <- function(tab) {
  ep <- tab$episode[tab$type == "temporary"]
  if (length(ep) == 0) 0L else max(table(ep))
}
