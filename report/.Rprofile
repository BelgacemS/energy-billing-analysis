# R reads the .Rprofile of the directory Quarto is launched from. When the
# report is rendered from report/ (e.g. from an IDE), activate the project's
# renv library from the root.
local({
  owd <- setwd("..")
  on.exit(setwd(owd))
  source("renv/activate.R")
})
