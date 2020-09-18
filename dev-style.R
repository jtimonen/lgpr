# This file is only for development purposes and not a part of the package
# itself. Run this in the package root directory to style the package code.
require(styler)

# Format the package code according to the tidyverse guide
styler:::style_pkg(exclude_files = 'R/stanmodels.R')
