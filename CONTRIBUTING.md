## lgpr development

- Style follows mostly the `tidyverse` style guide
- Use snake_case for variable and function names
- Line length should not exceed 80 characters
- Run `lint.R` to check for style problems (uses `lintr` and `styler` packages)
- The `lintr` part only checks for problems
- The `styler` part formats code
- All package R code should be in `./R`
- All package Stan code should be in `./inst/stan`
