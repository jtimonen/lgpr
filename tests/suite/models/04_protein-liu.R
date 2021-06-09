# A model for the longitudinal proteomics data set (Liu et. al, 2018)
library(lgpr)

# Load data
setup_data <- function(verbose) {
  protein <- "Q8WZA1"
  a <- read_proteomics_data(protein = protein, verbose = verbose)
}

# Create model
setup_model <- function(...) {
  dat <- setup_data(verbose = FALSE)
  f <- y ~ gp(age) * zs(id) + gp(age) + gp_vm(diseaseAge) +
    zs(sex) * gp(age) + zs(group)
  model <- create_model(formula = f, dat, ...)
  return(model)
}

expected_relevances <- function() {
  c(0.28581617, 0.14307953, 0.01981455, 0.07246394, 0.01652509, 0.46230072)
}
