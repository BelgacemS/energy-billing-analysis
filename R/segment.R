# Segmentation of one series into plateaus robust to outliers (step 5).

#' Robust noise scale of a series: MAD of the first differences / sqrt(2).
#' Differencing removes the plateaus and the MAD ignores the few large jumps
#' caused by breaks and spikes.
noise_scale <- function(x) stats::mad(diff(x)) / sqrt(2)

#' Biweight fit of a constant level: min over mu of sum(min((z - mu)^2, K^2))
#' (Fearnhead and Rigaill, 2019). Returns the minimal cost and the level.
#'
#' The loss is quadratic in mu between consecutive breakpoints z_i +- K, where
#' the set of points within K of mu is fixed: there the optimum is the mean of
#' that set, clamped to the interval. Taking the best interval gives the exact
#' minimum.
biweight_fit <- function(z, K) {
  bp <- sort(unique(c(z - K, z + K)))
  lo <- bp[-length(bp)]
  hi <- bp[-1]
  inl <- abs(outer((lo + hi) / 2, z, "-")) < K  # (intervals, points)
  n_in <- rowSums(inl)
  mu <- (inl %*% z)[, 1] / pmax(n_in, 1)
  mu <- pmin(pmax(mu, lo), hi)
  cost <- rowSums(pmin(outer(mu, z, "-")^2, K^2))
  cost[n_in == 0] <- Inf
  best <- which.min(cost)
  list(cost = cost[best], level = mu[best])
}

#' Optimal partitioning (exact dynamic programming) with the biweight cost:
#' minimises sum over segments of their biweight cost + penalty per segment,
#' with segments of at least `minseglen` points. Returns the end index of
#' each segment (the last one is length(z)). O(n^2) segment costs, fine for a
#' few dozen months.
biweight_partition <- function(z, K, penalty, minseglen = 1) {
  n <- length(z)
  if (n < 2 * minseglen) return(n)
  best <- c(0, rep(Inf, n))  # best[e + 1]: optimal cost of z[1:e]
  start <- integer(n)        # first point of the last segment of z[1:e]
  for (e in seq(minseglen, n)) {
    s <- seq_len(e - minseglen + 1)
    s <- s[s == 1 | s > minseglen]  # the previous segment must be long enough
    cand <- vapply(s, function(k) best[k] + biweight_fit(z[k:e], K)$cost + penalty, numeric(1))
    start[e] <- s[which.min(cand)]
    best[e + 1] <- min(cand)
  }
  ends <- n
  while (start[ends[1]] > 1) ends <- c(start[ends[1]] - 1L, ends)
  as.integer(ends)
}

#' End index of each segment of a standardised series (the last is length(z)).
#'
#' "biweight": biweight_partition with penalty 3 log(n) per change, the
#' penalty value of changepoint's MBIC for a change in mean, without the
#' log(segment length) terms its MBIC cost adds. "normal" (cpt.mean, which assumes
#' unit variance) and "np" (cpt.np, scale free) are kept for comparison; both
#' return the last observation of each segment (checked on a clean step, the
#' cpt.np help says otherwise) and fail when n < 2 * minseglen, where no change
#' is possible anyway.
find_changepoints <- function(z, cost = c("biweight", "normal", "np"), penalty = "MBIC",
                              minseglen = 1, K = 3, nquantiles = 10) {
  cost <- match.arg(cost)
  n <- length(z)
  if (n < 2 * minseglen) return(n)
  if (cost == "biweight") {
    if (penalty != "MBIC") stop("Only the MBIC penalty is implemented for the biweight cost")
    return(biweight_partition(z, K, 3 * log(n), minseglen))
  }
  fit <- if (cost == "normal") {
    changepoint::cpt.mean(z, penalty = penalty, method = "PELT",
                          test.stat = "Normal", minseglen = minseglen)
  } else {
    changepoint.np::cpt.np(z, penalty = penalty, method = "PELT",
                           minseglen = minseglen, nquantiles = nquantiles)
  }
  c(changepoint::cpts(fit), n)
}

#' Segment id (1, 2, ...) of each point from the segment end indices.
segment_ids <- function(ends) rep(seq_along(ends), diff(c(0, ends)))

#' Segmentation of a series into plateaus of at least `plateau_min` months.
#'
#' The series is divided by its noise scale, then segmented. With the
#' biweight cost a month weighs at most K^2, so an isolated bill, a few
#' scattered ones or a short ramp cannot create or move a plateau. The level of
#' a plateau is its biweight level (mean of the months within K noise units of
#' it). A month is an `excursion` when it lies more than K noise units from its
#' plateau, i.e. when the cost treats it as an outlier. Returns `months`
#' (plateau id, level, excursion flag), `breaks` (first month of each new
#' plateau, levels before and after) and `sigma`.
segment_series <- function(y, cost = "biweight", penalty = "MBIC", plateau_min = 3,
                           K = 3, nquantiles = 10) {
  if (anyNA(y)) stop("Missing values in the series")
  sigma <- noise_scale(y)
  if (!(sigma > 0)) stop("Zero noise scale: the series cannot be standardised")
  z <- y / sigma

  ends <- find_changepoints(z, cost, penalty, plateau_min, K, nquantiles)
  plateau <- segment_ids(ends)
  level_z <- vapply(split(z, plateau), function(v) biweight_fit(v, K)$level, numeric(1))
  level <- unname(level_z[plateau]) * sigma
  starts <- head(ends, -1) + 1L
  list(
    months = data.frame(plateau = plateau, level = level,
                        excursion = abs(y - level) / sigma > K),
    breaks = data.frame(index = starts,
                        level_before = unname(level_z[plateau[starts] - 1L]) * sigma,
                        level_after = unname(level_z[plateau[starts]]) * sigma),
    sigma = sigma
  )
}

#' Share of runs over synthetic series meeting each criterion of step 5.
#'
#' Scenarios, all with unit noise (Gaussian, or Student t with 3 df scaled to
#' unit variance): constant level; step of each size of `step_sizes` from
#' `step_month`; spike of each size of `spike_sizes` at `anomaly_month`;
#' compensated pair of each size a of `spike_sizes`: two months billed a
#' below the level from `anomaly_month`, then the catch-up 2a above. Criteria:
#' any lasting break (false break except for "step"), break found exactly at
#' the step month, all anomalous months flagged as excursions.
synthetic_segmentation <- function(n, reps, seed, step_sizes, spike_sizes,
                                   step_month, anomaly_month, seg_args) {
  set.seed(seed)
  grid <- merge(rbind(
    data.frame(scenario = "constant", size = 0),
    data.frame(scenario = "step", size = step_sizes),
    data.frame(scenario = "spike", size = spike_sizes),
    data.frame(scenario = "pair", size = spike_sizes)
  ), data.frame(noise = c("gaussian", "t3")), by = NULL)
  out <- lapply(seq_len(nrow(grid)), function(k) {
    g <- grid[k, ]
    mu <- numeric(n)
    anomalous <- integer()
    if (g$scenario == "step") mu[step_month:n] <- g$size
    if (g$scenario == "spike") anomalous <- anomaly_month
    if (g$scenario == "pair") anomalous <- anomaly_month + 0:2
    if (g$scenario == "spike") mu[anomalous] <- g$size
    if (g$scenario == "pair") mu[anomalous] <- c(-1, -1, 2) * g$size
    runs <- vapply(seq_len(reps), function(r) {
      e <- if (g$noise == "gaussian") stats::rnorm(n) else stats::rt(n, 3) / sqrt(3)
      fit <- do.call(segment_series, c(list(mu + e), seg_args))
      c(any_break = nrow(fit$breaks) > 0,
        break_at_step = identical(fit$breaks$index, as.integer(step_month)),
        anomaly_flagged = length(anomalous) > 0 && all(fit$months$excursion[anomalous]))
    }, logical(3))
    cbind(g, t(rowMeans(runs)))
  })
  do.call(rbind, out)
}
