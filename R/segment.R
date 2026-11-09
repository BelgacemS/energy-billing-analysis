# Segmentation of one series into plateaus and short excursions (step 5).

#' Robust noise scale of a series: MAD of the first differences / sqrt(2).
#' Differencing removes the plateaus and the MAD ignores the few large jumps
#' caused by breaks and spikes.
noise_scale <- function(x) stats::mad(diff(x)) / sqrt(2)

#' End index of each segment found by PELT (the last one is length(z)).
#'
#' `z` must be standardised: the Normal cost of cpt.mean assumes unit
#' variance. cpt.np uses the empirical distribution and is scale free. Both
#' return the last observation of each segment except the final one (checked
#' on a clean step, the cpt.np help says otherwise).
find_changepoints <- function(z, cost = c("normal", "np"), penalty = "MBIC",
                              minseglen = 1, nquantiles = 10) {
  cost <- match.arg(cost)
  fit <- if (cost == "normal") {
    changepoint::cpt.mean(z, penalty = penalty, method = "PELT",
                          test.stat = "Normal", minseglen = minseglen)
  } else {
    changepoint.np::cpt.np(z, penalty = penalty, method = "PELT",
                           minseglen = minseglen, nquantiles = nquantiles)
  }
  c(changepoint::cpts(fit), length(z))
}

#' Segment id (1, 2, ...) of each point from the segment end indices.
segment_ids <- function(ends) rep(seq_along(ends), diff(c(0, ends)))

#' Two-pass segmentation of a series into plateaus.
#'
#' Pass 1 (minimal segment length `minseglen`) lets a spike form its own
#' segment instead of moving a plateau; runs of consecutive segments shorter
#' than `plateau_min`, lasting less than `plateau_min` months in all, are
#' marked as excursions; a longer run enters pass 2 through its median.
#' Pass 2 segments the remaining months with
#' segments of at least `plateau_min` months: its segments are the plateaus.
#' The level of a plateau is the median of its months (robust to the small
#' spikes pass 1 left inside). An excursion month takes the level of the
#' previous plateau, or of the next one before the first plateau, so that
#' its residual measures the gap to the usual level.
#'
#' Returns `months` (one row per point: segment of pass 1, excursion flag,
#' plateau id, level), `breaks` (first month of each new plateau with the
#' levels before and after) and `sigma`.
segment_series <- function(y, cost = "normal", penalty = "MBIC", minseglen = 1,
                           plateau_min = 3, nquantiles = 10) {
  if (anyNA(y)) stop("Missing values in the series")
  sigma <- noise_scale(y)
  if (!(sigma > 0)) stop("Zero noise scale: the series cannot be standardised")
  z <- y / sigma
  n <- length(y)

  seg1 <- segment_ids(find_changepoints(z, cost, penalty, minseglen, nquantiles))
  short <- tabulate(seg1)[seg1] < plateau_min
  # Consecutive short segments form one run; a short run is an excursion. A
  # long run is a stretch where no level holds for plateau_min months: a
  # regime noisier than the site's usual noise (1240BR after its drop) or a
  # cluster of spikes around the usual level (0089KO in mid-2019). Pass 2
  # sees it through its median, which keeps the level of a regime and
  # ignores the spikes of a cluster as long as they are fewer than half of
  # its months.
  run <- cumsum(c(TRUE, diff(short) != 0))
  excursion <- short & tabulate(run)[run] < plateau_min
  long_run <- short & !excursion
  z2 <- z
  z2[long_run] <- stats::ave(z[long_run], run[long_run], FUN = stats::median)
  keep <- which(!excursion)
  if (length(keep) < plateau_min) stop("Too few months outside excursions")

  ends2 <- find_changepoints(z2[keep], cost, penalty, plateau_min, nquantiles)
  plateau <- rep(NA_integer_, n)
  plateau[keep] <- segment_ids(ends2)
  # Excursions: previous plateau, or the first one at the start of the series.
  filled <- cummax(ifelse(is.na(plateau), 0L, plateau))
  plateau <- ifelse(filled == 0L, 1L, filled)

  level_of <- vapply(split(y[keep], plateau[keep]), stats::median, numeric(1))
  starts <- keep[c(FALSE, diff(plateau[keep]) != 0)]
  list(
    months = data.frame(segment = seg1, excursion = excursion, plateau = plateau,
                        level = unname(level_of[plateau])),
    breaks = data.frame(index = starts,
                        level_before = unname(level_of[plateau[starts] - 1L]),
                        level_after = unname(level_of[plateau[starts]])),
    sigma = sigma
  )
}

#' Share of runs over synthetic series meeting each criterion of step 5.
#'
#' Scenarios: constant level, step of each size of `step_sizes` starting at
#' `step_month`, spike of each size of `spike_sizes` at `spike_month`, all
#' with unit noise (Gaussian, or Student t with 3 df scaled to unit
#' variance). Criteria: any lasting break (false break for "constant" and
#' "spike"), break found exactly at the step month, spike month marked as an
#' excursion.
synthetic_segmentation <- function(n, reps, seed, step_sizes, spike_sizes,
                                   step_month, spike_month, seg_args) {
  set.seed(seed)
  scen <- rbind(
    data.frame(scenario = "constant", size = 0),
    data.frame(scenario = "step", size = step_sizes),
    data.frame(scenario = "spike", size = spike_sizes)
  )
  grid <- merge(scen, data.frame(noise = c("gaussian", "t3")), by = NULL)
  out <- lapply(seq_len(nrow(grid)), function(k) {
    g <- grid[k, ]
    mu <- numeric(n)
    if (g$scenario == "step") mu[step_month:n] <- g$size
    if (g$scenario == "spike") mu[spike_month] <- g$size
    runs <- vapply(seq_len(reps), function(r) {
      e <- if (g$noise == "gaussian") stats::rnorm(n) else stats::rt(n, 3) / sqrt(3)
      fit <- do.call(segment_series, c(list(mu + e), seg_args))
      c(any_break = nrow(fit$breaks) > 0,
        break_at_step = identical(fit$breaks$index, as.integer(step_month)),
        spike_excursion = fit$months$excursion[spike_month])
    }, logical(3))
    cbind(g, t(rowMeans(runs)))
  })
  do.call(rbind, out)
}
