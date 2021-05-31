# A model with two age components
library(lgpr)
library(nlme)

# Load data
setup_data <- function() {
  as.data.frame(nlme::Orthodont)
}

# Create model
setup_model <- function(...) {
  dat <- setup_data()
  model <- create_model(distance ~ age + age | Subject + age | Sex, dat, ...)
  return(model)
}
