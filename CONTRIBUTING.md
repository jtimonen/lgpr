# lgpr development

## Organization
- All package R code should be in `./R`
- All package Stan code should be in `./inst/stan`

## Style
- In general, follow the `tidyverse` style guide.
- Use snake_case for variable and function names. An exception is one-letter
  variable names, which can be capitalized. In some cases, also longer all caps    variable names are OK.
- Use snake_case also for the code file names. Some exceptions are files which
  are auto-generated, like `RcppExports.R`, and files that define only S4
  classes or their methods, like `AllClasses.R`.
- Line length should not exceed 80 characters, with the exception of long URLs.

## Workflow
- Document each new function with roxygen
- After making changes to `.R` code, first run `styler::style_pkg` to
  correct most style problems. Then run `dev-lint.R` to check for other
  style problems (uses the `lintr` package) and correct them manually if found.
- After making changes to `.stan` code, run `dev-cpp.R` and rebuild the C++
  code (for example by `clean and rebuild` in Rstudio).
- Before making new commits, make sure to also run `devtools::document()` and 
  `R CMD check`

