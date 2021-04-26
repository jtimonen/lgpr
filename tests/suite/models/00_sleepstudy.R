# A model with two age components
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
