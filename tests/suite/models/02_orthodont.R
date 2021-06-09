# A model for Orthodont data
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

expected_relevances <- function() {
  c(0.3091031, 0.3736700, 0.2065741, 0.1106528)
}
