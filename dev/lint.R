require(lintr)

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
