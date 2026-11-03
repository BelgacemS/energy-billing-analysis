# Project paths and parameters. All paths are relative to the project root.
# The data path can be overridden with the DATA_PATH environment variable.
config <- list(
  data_path   = Sys.getenv("DATA_PATH", "data/TimeSeries.csv"),
  output_dir  = "outputs",
  results_dir = file.path("outputs", "results")
)
