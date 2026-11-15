# Project paths and parameters. All paths are relative to the project root.
# The data path can be overridden with the DATA_PATH environment variable.
config <- list(
  data_path   = Sys.getenv("DATA_PATH", "data/TimeSeries.csv"),
  output_dir  = "outputs",
  results_dir = file.path("outputs", "results"),
  figures_dir = file.path("outputs", "figures"),
  tables_dir  = file.path("outputs", "tables"),

  # Descriptive statistics (step 3).
  # Median monthly kWh above which a site is treated apart. The panel has a gap
  # between about 6,000 and 350,000 kWh, so any value in between gives 2 sites.
  big_site_kwh = 1e5,
  # Level change: largest gap between the medians of the `window` months before
  # and after a month, on log kWh per day, flagged above log(1 + min).
  level_shift_window = 6,
  level_shift_min    = 0.10,

  # First month of the test period of the backtest (June to November 2021).
  # Forecast accuracy figures computed before step 11 stop the month before.
  test_start = as.Date("2021-06-01"),

  # Baseline (step 4).
  # Identical bills two months in a row are compared with the mean count of
  # changes equal to k kWh, k in +-1..+-repeat_k_max.
  repeat_k_max = 10,
  # Hampel filter: window of 3 months on each side, flag beyond 3 robust
  # standard deviations (the usual default values).
  hampel_half_width = 3,
  hampel_threshold  = 3,
  forecast_month    = as.Date("2021-12-01"),

  # Course view (step 3 bis).
  # Three distinct ordinary sites, picked before looking at their plots:
  # stable = lowest noise (MAD of first differences / sqrt(2)) among sites
  # without level change; break = largest level change; spike = largest
  # Hampel |z| outside the window edges, among sites without level change.
  # Level changes and Hampel scores use the step 3 and step 4 parameters.
  season_period = 12,
  stl_s_window  = "periodic",  # 35 months: too short to let the season evolve
  # Non-robust STL, as R and nsdiffs do by default: the robust version sends a
  # break near the end of a series (1240BR) into the remainder, which wipes out
  # the trend and season strengths.
  stl_robust    = FALSE,
  # Season strength above which nsdiffs' default test asks for D = 1 (fixed in
  # the forecast package).
  nsdiffs_seas_threshold = 0.64,

  # Panel and segmentation (step 5). K and the penalty were set by the
  # calibration of step 7 (scripts/05_calibrate.R checks they match its rule).
  # Month effect: median deviation from the centred 2x12 moving average (trend
  # of decompose), shared by the ordinary sites; each big site keeps its own
  # (step 3).
  # Segmentation of log kWh per day minus the month effect, divided by the
  # site noise (MAD of first differences / sqrt(2)), into plateaus of at least
  # plateau_min_months. Biweight cost (Fearnhead and Rigaill, 2019): a month
  # weighs at most K^2, so isolated bills, compensated pairs and short ramps
  # cannot create or move a plateau. Exact optimal partitioning, penalty
  # 3 log(n) per change (the penalty value of changepoint's MBIC for a mean).
  # A Gaussian two-pass version (one-month segments, then plateaus) was
  # dropped: an under-billing
  # over two months followed by its catch-up created a false plateau in every
  # simulated series. "normal" (cpt.mean) and "np" (cpt.np) stay available for
  # the synthetic comparison.
  segment_cost       = "biweight",
  segment_penalty    = "MBIC",
  segment_K          = 6,         # outlier threshold, in noise units (3 until step 7)
  plateau_min_months = 3,         # a plateau lasts at least a quarter
  segment_nquantiles = 10,        # cpt.np default
  # Synthetic check of the three costs (step 5): 35 months, unit noise
  # (Gaussian or Student t with 3 df scaled to unit variance), amplitudes in
  # noise units. Spikes and compensated pairs (two months a below the level,
  # then 2a above) use synthetic_spike_sizes.
  synthetic_reps          = 200,
  synthetic_seed          = 2021,
  synthetic_step_sizes    = c(1, 2, 3),
  synthetic_spike_sizes   = c(3, 5, 8),
  synthetic_step_month    = 21,   # first month of the new level
  synthetic_anomaly_month = 15,

  # Typology of the flagged bills (step 6). Rules written before the first
  # run, from the step 5 outputs only. A bill is flagged when its period is
  # irregular (step 4 check), when it lies more than segment_K noise units from
  # its plateau (excursion), or when it is the first month of a plateau that
  # differs by more than level_shift_min from the previous one. Breaks of
  # level_shift_min or less are ignored until step 7: plateaus separated only
  # by such breaks count as one level.
  # Episodes: runs of consecutive excursion months, merged when separated by at
  # most merge_gap_months normal months. Types, in order of precedence:
  # irregular period; indeterminate (episode reaching the last edge_months
  # months); transition (episode touching a break, ending at most one month
  # before it, all months between the old and the new level within segment_K
  # noise units: a ramp, shown but not suspect); compensated pair (bills below
  # and above, in either order and possibly more than two, |sum billed - sum
  # expected| / sum expected within pair_tolerance at the level before or after
  # the episode: a billing error corrected later; its net gap is taken at the
  # closer of the two levels); then for each
  # run: point (1 month, back to the same level), implausible (a point billed
  # above implausible_ratio or below 1 / implausible_ratio times the expected
  # value), temporary (2 to temporary_max_months months of the same sign, back
  # to the same level; longer runs that come back stay temporary, so the cap
  # does not change the type: it is only a readability guard, checked by
  # scripts/04_suspects.R, and no step 5 run exceeds it), indeterminate
  # otherwise. A break is a lasting change (its cumulative gap at the old level
  # starts with the transition that leads to it, if any), a probable network
  # change when shared (not suspect), or indeterminate in the last edge_months
  # months.
  # Shared break, against chance (same logic as check_repeats): for each
  # direction, lambda = number of breaks above level_shift_min of the ordinary
  # sites / number of testable months (months where a new plateau can start
  # given plateau_min_months: months plateau_min_months + 1 to n -
  # plateau_min_months + 1, 30 of 35). A month is shared in a direction when
  # its n breaks satisfy P(X >= n) < shared_break_alpha, X ~ Poisson(lambda),
  # and n >= shared_break_min_sites. At least 3 sites in the same month
  # matched chance (about 7 months expected, 7 observed).
  pair_tolerance       = 0.10,
  implausible_ratio    = 3,
  shared_break_min_sites = 3,
  shared_break_alpha   = 0.05,
  edge_months          = 3,   # first and last months, where a plateau is less reliable
  merge_gap_months     = 1,
  temporary_max_months = 4,
  # Confidence: high for an irregular period (deterministic check); otherwise
  # high when |z| >= high_confidence_z outside the edge months and outside a
  # shared break, low in the edge months, medium otherwise. The
  # range of a gap is billed - expected x exp(+-range_sigmas x noise).
  high_confidence_z    = 6,
  range_sigmas         = 2,

  # Backtest (step 8). Rolling origin, one month ahead: the forecast of month
  # t uses the bills up to t - 1 only, and any preprocessing is refitted at
  # each origin. Validation targets run from validation_start to the month
  # before test_start (July 2020 to May 2021, 11 forecasts per site, December
  # 2020 included); the test targets (June to November 2021) are read once,
  # at step 11.
  # Error of a forecast: e = billed - forecast, in kWh. MASE (Hyndman and
  # Koehler, 2006): |e| / mean |e| of the per-day naive forecast over the
  # bills available at the origin (one-step, in sample), so it is free of
  # scale across sites and the per-day naive sits near 1.
  # Selection rule, written before any model is run: lowest mean validation
  # MASE; among the methods less than selection_tolerance (relative) above
  # that lowest value, the simplest one, simplicity following method_order. The
  # per-day naive comes first: it is the reference of the MASE and works in
  # the unit of analysis. The candidates of step 10 (class model, robust local
  # level, combination) go after those of step 9. The rule is applied by class
  # at step 10, never site by site (11 validation forecasts per site are too
  # few).
  validation_start    = as.Date("2020-07-01"),
  selection_tolerance = 0.02,
  method_order = c("naive_per_day", "naive_raw", "snaive_per_day", "ses", "ets",
                   "arima_noseason", "arima", "naive_adj", "robust_level", "combination"),

  # Course models (step 9), fitted site by site at each origin on log kWh per
  # calendar day (monthly ts). The one-step forecast is exp(log forecast) x
  # days of the target month: the median of the forecast distribution, which
  # suits absolute errors (MASE). A fit that fails falls back to the per-day
  # naive and is counted.
  # ses: forecast::ses. ets: forecast::ets, additive components only (on the
  # log scale a multiplicative component would duplicate the log).
  # arima: auto.arima with its defaults (seasonal difference chosen by the
  # "seas" test). arima_noseason: same with D = 0 (no seasonal difference,
  # seasonal AR and MA terms still allowed) for every site except
  # arima_season_sites, where the default choice is kept (0698BR, the only
  # series where the course tests agree on a seasonal difference, step 3 bis).
  ets_additive_only  = TRUE,
  arima_season_sites = "0698BR",
  # Diagnosis of auto.arima, both variants, fitted once per site on the months
  # before test_start. A site is listed when:
  # - degenerate order: p = q = P = Q = 0, no AR or MA term: a random walk
  #   (d = 1; without drift it is the per-day naive on the log scale), a
  #   seasonal random walk (D = 1), a constant mean (d = 0) or d = 2;
  # - differencing caused by breaks: the site has a level break above
  #   level_shift_min (step 5 segmentation refitted on the same months) and,
  #   once the mean of each regime between such breaks is removed, the KPSS
  #   test (ndiffs) asks for fewer differences than d, or nsdiffs for none
  #   when D = 1;
  # - autocorrelated residuals: Ljung-Box p-value below ljung_box_alpha, with
  #   h = 10 lags (2 x 12 for a seasonal model) capped at T / 5 as in the
  #   course note (FPP caps h at T / 5; rounded down here), and fitdf =
  #   p + q + P + Q (no test when h <= fitdf);
  # - a fit that failed during the backtest.
  ljung_box_alpha    = 0.05,
  # Added after the review of step 9, rule above unchanged: rate at which the
  # break criterion fires on random walks, for each ordinary site, under two
  # nulls that bracket the truth: Gaussian increments with sd = MAD of its
  # monthly log differences (the MAD drops the large jumps, so too few
  # breaks) and its own increments resampled with replacement (the jumps of
  # real breaks count as ordinary increments, so too many).
  break_null_reps = 10,
  break_null_seed = 2021,

  # Calibration by anomaly injection (step 7), rule written before the first
  # run. Simulated panels of the 90 ordinary sites: for each site, levels
  # between its breaks above level_shift_min and month effect from the
  # reference segmentation (the step 5 settings below), plus its residuals
  # around those levels (months off their plateau left out, as they may be
  # real anomalies) resampled in moving blocks of calibration_block months.
  # Only breaks above level_shift_min are true breaks: the smaller ones are
  # what the calibration has to judge, so they stay in the residuals, whose
  # lag-1 autocorrelation is then about 0.75 (0.09 at lag 5, hence 6-month
  # blocks). Each replication draws one panel per seed (calibration_seed +
  # replication) and one start month per site, shared by every anomaly and
  # setting (common random numbers); the anomaly is injected in every site,
  # outside the first and last edge_months months. Anomalies, multiplying the
  # bills: spike (1 month x (1 + a) or / (1 + a)), temporary gap
  # (calibration_episode_months months, same factor), compensated pair (two
  # months x (1 - a) then one x (1 + 2a), or the reverse order), level break
  # (x (1 + a) or / (1 + a) from the start month on).
  # Measures, for each setting (K, beta: penalty beta log n per change):
  # - false alert rate: suspect rows of the background panel (no injection)
  #   that are not one of its true breaks found within calibration_break_window
  #   months in the same direction, over the number of bills;
  # - detection: spike, temporary gap and pair found when one of their months
  #   is a suspect row; break found when a break row above level_shift_min
  #   lies within calibration_break_window months, same direction;
  # - typology: share of detections with the injected type (break: lasting
  #   or network change).
  # Rule: among the settings whose false alert rate is at most
  # calibration_fa_max, the highest score, score = mean of the mean detection
  # over the spike cells and the mean over the break cells (sizes and both
  # signs, equal weight); scores within calibration_tie of the best are tied
  # and the lowest false alert rate wins. If no setting is admissible, the
  # lowest false alert rate. Temporary gaps and pairs are measured at the
  # chosen setting only. Hampel (step 4) runs on the same panels for each of
  # hampel_thresholds; it is compared at its default threshold and at the
  # smallest threshold whose false alert rate does not exceed the panel's.
  calibration_reference   = list(K = 3, penalty = "MBIC"),
  calibration_K           = c(2.5, 3, 4, 5, 6),
  calibration_beta        = c(3, 5, 8, 12),
  calibration_reps        = 20,
  calibration_seed        = 2021,
  calibration_block       = 6,
  calibration_fa_max      = 0.01,
  calibration_tie         = 0.005,
  calibration_spike_sizes = c(0.01, 0.02, 0.05, 0.10, 0.20, 0.50, 1.00),
  calibration_break_sizes = c(0.05, 0.10, 0.20, 0.50, 1.00),
  calibration_episode_sizes  = c(0.05, 0.10, 0.20, 0.40),  # pair needs 2a < 1
  calibration_episode_months = 3,
  calibration_break_window   = 1,
  hampel_thresholds       = c(3, 4, 5, 6, 8),

  # Classes and selection (step 10), rules written before the first run.
  # Classes of the ordinary sites by Ward clustering (hclust "ward.D2",
  # Euclidean distance on standardised features) on the months before
  # validation_start only, so the classes use no validation or test month.
  # Features of log kWh per day minus the month effect: median level, OLS
  # slope per year, log noise scale (step 5), number of breaks above
  # level_shift_min (segmentation of R/config.R). Number of classes: highest
  # mean silhouette over class_k_range, the smaller k on a tie. The two big
  # sites form their own class.
  # Extra candidates, refitted at each origin on the bills up to the origin
  # (preprocessing of the backtest): month effect = median over the ordinary
  # sites (each big site its own) of the deviation from the site's median of
  # each complete calendar year (the decomposition trend of step 5 needs more
  # months than the first origins have). Added after the first run, which
  # stopped: with one complete year a big site's own effect copies its data,
  # so the big sites take the common effect until two years are complete.
  # - naive_adj: per-day naive corrected by the month effect, last log kWh
  #   per day - effect of its month + effect of the target month;
  # - robust_level: months off their plateau replaced by the plateau level
  #   (segmentation of R/config.R on log kWh per day minus the month effect),
  #   simple exponential smoothing of that cleaned series from the start of
  #   the last plateau after a break above level_shift_min, plus the month
  #   effect of the target month; the last cleaned value if the fit fails;
  # - combination: mean of the kWh forecasts of combination_methods (equal
  #   weights, Bates and Granger, 1969), chosen before the run: the per-day
  #   naive and the two candidates above.
  # Selection rule of step 8 (method_order, selection_tolerance) applied to
  # the mean validation MASE of each class, step 9 methods included.
  class_k_range       = 2:6,
  combination_methods = c("naive_per_day", "naive_adj", "robust_level")
)
