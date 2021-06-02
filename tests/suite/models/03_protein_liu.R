# A model for the longitudinal proteomics data set (Liu et. al, 2018)
library(lgpr)

# Load data
setup_data <- function() {
  a <- read_proteomics_data(protein = 450, verbose = TRUE)
}

# Create model
setup_model <- function(...) {
  dat <- setup_data()
  model <- create_model(y ~ age + age | id + age | sex, dat, ...)
  return(model)
}
