# This file is only for development purposes and not a part of the package
# itself. Run this in the package root directory to (style and) lint the
# package.
#
# Author: Juho Timonen
#
require(lintr)
require(styler)

# Style the package according to the tidyverse guide
# styler::style_pkg()

# Specify linters
linters <- lintr::with_defaults(
  object_name_linter = NULL,
  open_curly_linter = NULL
  )

# Files that are not linted
exclusions <- list("R/RcppExports.R", "R/stanmodels.R")

# Lint the package
lout <- lintr::lint_package(linters = linters, exclusions = exclusions)
show(summary(lout))

# Lint a single file
#lintr::lint('R/all-classes.R', linters = linters)
