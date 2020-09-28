require(styler)

# Format the package code according to the tidyverse guide
exclusions <- list("R/RcppExports.R", "R/stanmodels.R")
styler:::style_pkg(exclude_files = exclusions)
