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
  c(0.2818898, 0.4012653, 0.1506231, 0.1662218)
}
