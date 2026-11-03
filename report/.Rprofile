# Quarto starts R from report/, so the root .Rprofile is not read.
# Activate the project's renv library explicitly.
local({
  owd <- setwd("..")
  on.exit(setwd(owd))
  source("renv/activate.R")
})
