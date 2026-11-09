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
  # Segmentation by PELT on log kWh per day minus the month effect, divided by
  # the site noise (MAD of first differences / sqrt(2)): the Normal cost of
  # cpt.mean assumes unit variance. Pass 1 uses one-month segments so that a
  # spike forms its own segment; runs of short segments lasting less than
  # plateau_min_months in all are excursions (a longer run, noisier regime or
  # cluster of spikes, enters pass 2 through its median). Pass 2 segments the
  # series without the excursions: its segments are the plateaus and its
  # changes the lasting breaks.
  segment_cost       = "normal",  # "normal" (cpt.mean) or "np" (cpt.np)
  segment_penalty    = "MBIC",    # changepoint default, 3 log(n) per change
  segment_minseglen  = 1,
  plateau_min_months = 3,         # a plateau lasts at least a quarter
  segment_nquantiles = 10,        # cpt.np default
  # Synthetic check of both costs (step 5): 35 months, unit noise (Gaussian or
  # Student t with 3 df scaled to unit variance), amplitudes in noise units.
  synthetic_reps        = 200,
  synthetic_seed        = 2021,
  synthetic_step_sizes  = c(1, 2, 3),
  synthetic_spike_sizes = c(3, 5, 8),
  synthetic_step_month  = 21,     # first month of the new level
  synthetic_spike_month = 15
)
