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

  # Panel and segmentation (step 5). Provisional values: penalty and cost are
  # chosen at step 7 by anomaly injection.
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
  segment_K          = 3,         # outlier threshold, in noise units
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
  range_sigmas         = 2
)
