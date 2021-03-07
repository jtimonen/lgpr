# https://github.com/stan-dev/rstanarm/blob/master/tests/testthat/helpers/SW.R
SW <- function(expr) {
  capture.output(suppressWarnings(expr))
}
