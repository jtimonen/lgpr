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
  f <- y ~ gp(age) * zs(id) + gp(age) + het(id) * gp_vm(diseaseAge) +
    zs(sex) * gp(age) + zs(group)
  model <- create_model(formula = f, dat, ...)
  return(model)
}

expected_relevances <- function() {
  c(0.270712119, 0.171311656, 0.008459786, 0.049614305, 0.013833720, 0.486068414)
}
