# A model for the sleepstudy data
library(lgpr)
library(lme4)

# Load data
setup_data <- function() {
  lme4::sleepstudy
}

# Create model
setup_model <- function(...) {
  dat <- setup_data()
  model <- create_model(Reaction ~ Days + Days | Subject, dat, ...)
  return(model)
}

expected_relevances <- function() {
  c(0.3160955, 0.5442927, 0.1396118)
}
