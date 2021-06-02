# A model for the longitudinal proteomics data set (Liu et. al, 2018)
library(lgpr)
library(nlme)

# Load data
setup_data <- function() {
  stop("Data path not set!")
}

# Create model
setup_model <- function(...) {
  dat <- setup_data()
  model <- create_model(distance ~ age + age | Subject + age | Sex, dat, ...)
  return(model)
}
