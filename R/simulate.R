# Simulated panels with injected anomalies, to measure detection and false
# alerts of the panel method and of the Hampel filter (step 7). Settings and
# the choice rule are in R/config.R.

#' Ingredients of the simulated panels, one list per ordinary site.
#'
#' `pb`: panel_fit()$bills of the reference segmentation; `breaks`: its
#' breaks. Level = median of y (log kWh per day minus month effect) over the
#' months on their plateau between two breaks above `min_change`; residuals =
#' y - level on those months, in time order. `break_index`: first month of
#' each true break.
simulation_parts <- function(pb, breaks, big, min_change) {
  ord <- pb[!pb$site %in% big, ]
  ord <- ord[order(ord$site, ord$month), ]
  large <- breaks[abs(breaks$change) > min_change, ]
  lapply(split(ord, ord$site), function(d) {
    idx <- match(large$month[large$site == d$site[1]], d$month)
    regime <- cumsum(seq_len(nrow(d)) %in% idx)
    ok <- !d$excursion
    level <- stats::ave(ifelse(ok, d$y, NA), regime, FUN = function(v) stats::median(v, na.rm = TRUE))
    list(site = d$site[1], suffix = d$suffix[1], month = d$month, month_days = d$month_days,
         effect = d$effect, level = level, resid = (d$y - level)[ok], break_index = sort(idx))
  })
}

#' Moving block bootstrap: blocks of `block` consecutive values of `r` with
#' uniform start, concatenated and cut to length n.
block_bootstrap <- function(r, n, block) {
  block <- min(block, length(r))
  starts <- sample.int(length(r) - block + 1, ceiling(n / block), replace = TRUE)
  r[as.vector(outer(0:(block - 1), starts, "+"))][seq_len(n)]
}

#' One simulated panel without anomaly: level + resampled residuals + month
#' effect, as a bills table (site, suffix, month, month_days, kwh,
#' kwh_per_day).
simulate_panel <- function(parts, block) {
  do.call(rbind, lapply(parts, function(p) {
    n <- length(p$month)
    rate <- exp(p$level + block_bootstrap(p$resid, n, block) + p$effect)
    data.frame(site = p$site, suffix = p$suffix, month = p$month, month_days = p$month_days,
               kwh = rate * p$month_days, kwh_per_day = rate, stringsAsFactors = FALSE)
  }))
}

#' Multiplicative factors of an anomaly from its start month: spike (1 month),
#' temporary (`months` months), pair (`months` - 1 months x (1 - a) then one
#' x (1 + (months - 1) a), or the reverse order) and break (to the end of the
#' series, `length_left` months). `sign` = +1 or -1 (for a pair: +1 starts
#' below, as an under-billing then its catch-up).
anomaly_factors <- function(type, size, sign, months, length_left) {
  switch(type,
    spike = (1 + size)^sign,
    temporary = rep((1 + size)^sign, months),
    pair = c(rep(1 - sign * size, months - 1), 1 + sign * (months - 1) * size),
    level_break = rep((1 + size)^sign, length_left),
    stop("Unknown anomaly type: ", type))
}

#' Inject one anomaly per site of a bills table (sorted by site, then month)
#' at the start indices `start` (named by site). Returns the bills.
inject_anomaly <- function(bills, type, size, sign, start, months) {
  for (s in names(start)) {
    rows <- which(bills$site == s)
    n <- length(rows)
    f <- anomaly_factors(type, size, sign, months, n - start[[s]] + 1)
    i <- rows[start[[s]] + seq_along(f) - 1L]
    bills$kwh[i] <- bills$kwh[i] * f
    bills$kwh_per_day[i] <- bills$kwh_per_day[i] * f
  }
  bills
}

#' Start index of the anomaly of each site from uniform draws `u`: between
#' month edge + 1 and the last start that keeps `span` months outside the
#' last `edge` months.
anomaly_start <- function(u, n, edge, span) {
  lo <- edge + 1L
  hi <- n - edge - span + 1L
  lo + floor(u * (hi - lo + 1L))
}

#' Panel detector on a bills table for one setting: month effect (shared,
#' estimated once by the caller), segmentation with K and beta, typology.
#' Returns the classify_bills table.
detect_panel <- function(bills, effect, K, beta, cfg) {
  seg <- list(cost = "biweight", penalty = beta, plateau_min = cfg$plateau_min_months, K = K)
  pf <- panel_fit(bills, effect, character(), seg)
  cfg$segment_K <- K
  periods <- data.frame(site = character(), month = as.Date(character()))
  classify_bills(pf$bills, pf$breaks, periods, character(), cfg)
}

#' Result of each injection in a classify_bills table: one row per site with
#' `detected`, `type` (type given to the injected months: "break" when a
#' break row lies on them, else the most frequent excursion type, NA when
#' not flagged) and `as_spike` (injected break not found, but its first or
#' second month flagged as an excursion). `start`: start index by site;
#' `months`: injected span for spikes, gaps and pairs.
score_injection <- function(flagged, bills, type, sign, start, months, window, min_change) {
  months_of <- split(bills$month, bills$site)
  sus <- flagged[flagged$suspect | flagged$type == "network_change", ]
  is_break <- is.na(sus$episode)
  out <- lapply(names(start), function(s) {
    m <- months_of[[s]]
    st <- start[[s]]
    f <- sus[sus$site == s, ]
    brk <- is_break[sus$site == s]
    if (type == "level_break") {
      near <- brk & abs(match(f$month, m) - st) <= window & sign * f$change > 0 &
        abs(f$change) > min_change
      spike_like <- !brk & f$month %in% m[st + 0:1]
      return(data.frame(site = s, detected = any(near), type = if (any(near)) "break" else NA,
                        as_spike = !any(near) && any(spike_like)))
    }
    span <- m[st + seq_len(if (type == "spike") 1 else months) - 1L]
    hit <- f$month %in% span & f$suspect
    t <- if (!any(hit)) NA_character_ else if (any(brk[hit])) "break" else
      names(which.max(table(f$type[hit & !brk])))
    data.frame(site = s, detected = any(hit), type = t, as_spike = NA)
  })
  do.call(rbind, out)
}

#' False alerts of a classify_bills table on a background panel: suspect rows
#' that are not a true break (`parts`) found within `window` months in the
#' same direction. Returns the count of false alert rows.
false_alerts <- function(flagged, parts, window) {
  sus <- flagged[flagged$suspect, ]
  true_break <- vapply(seq_len(nrow(sus)), function(i) {
    if (!is.na(sus$episode[i]) || is.na(sus$change[i])) return(FALSE)
    p <- parts[[sus$site[i]]]
    j <- match(sus$month[i], p$month)
    tb <- p$break_index[abs(p$break_index - j) <= window]
    any(sign(p$level[tb] - p$level[tb - 1]) == sign(sus$change[i]))
  }, logical(1))
  sum(!true_break)
}

#' Hampel flags of the injected months for each threshold: share of sites
#' whose injected span has |z| above the threshold (spike, gap, pair), or
#' whose first two months of a break do (`as_spike`).
hampel_hits <- function(hz, bills, type, start, months, thresholds) {
  sites <- names(start)
  span_len <- if (type == "spike") 1L else if (type == "level_break") 2L else months
  zmax <- vapply(sites, function(s) {
    z <- hz[bills$site == s]
    max(abs(z[start[[s]] + seq_len(span_len) - 1L]))
  }, numeric(1))
  vapply(thresholds, function(h) sum(zmax > h), numeric(1))
}

#' Index of the setting chosen by the rule of R/config.R: among the settings
#' with false alert rate `fa` at most `fa_max`, the highest `score`; scores
#' within `tie` of the best are tied and the lowest `fa` wins (then the first
#' setting). Without any admissible setting, the lowest `fa`.
choose_setting <- function(fa, score, fa_max, tie) {
  ok <- which(fa <= fa_max)
  if (length(ok) == 0) return(which.min(fa))
  top <- ok[score[ok] >= max(score[ok]) - tie]
  top[which.min(fa[top])]
}
