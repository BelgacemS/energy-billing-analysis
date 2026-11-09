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
  stl_s_window  = "periodic"  # 35 months: too short to let the season evolve
)
