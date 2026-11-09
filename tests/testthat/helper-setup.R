# test_dir() runs from tests/testthat, so paths are resolved from the project root.
project_root <- normalizePath(file.path("..", ".."))

for (f in list.files(file.path(project_root, "R"), pattern = "\\.R$", full.names = TRUE)) {
  source(f, encoding = "UTF-8")
}

# config paths are relative to the project root unless absolute (DATA_PATH).
root_path <- function(path) {
  if (grepl("^(/|~|[A-Za-z]:)", path)) path else file.path(project_root, path)
}
