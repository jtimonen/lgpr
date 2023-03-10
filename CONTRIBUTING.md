# lgpr development

## Organization
- All package R code should be in `./R`
- All package Stan code should be in `./inst/stan`

## Style
- In general, follow the `tidyverse` style guide.
- Use mostly snake_case for variable, function and file names.
- Line length should not exceed 80 characters, with the exception of long URLs.

## Workflow
- Document each new function with roxygen comments
- After making changes to `.R` code,
  1. run `devtools::document()` to update the man pages
  2. source `dev/style.R` to correct most style problems
  3. source `dev/lint.R` to check for other style problems
  4. correct remaining problems manually if found
- After making changes to `.stan` code
  1. rebuild the whole package

