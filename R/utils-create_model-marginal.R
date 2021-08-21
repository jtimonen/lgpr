# Create a marginal GP model from base model
create_model.marginal <- function(base_model, prior, verbose) {
  si_base <- get_stan_input(base_model)
  parsed <- standata_marginal(base_model, prior, verbose)
  si_add <- dollar(parsed, "to_stan")
  si <- c(si_base, si_add)
  new("MarginalGPModel",
    base_model,
    parsed_input = si,
    y_scaling = dollar(parsed, "y_scaling"),
    info = creation_info()
  )
}

# Create additional Stan input
standata_marginal <- function(base_model, prior, verbose) {

  # Prior for sigma
  filled <- fill_prior(prior, 0, "sigma")
  defaulting_info(filled, verbose)
  raw <- dollar(filled, "prior")
  stan_prior <- parse_prior_pos(raw, 1, "sigma")

  # Check that data contains the response variable which is numeric
  dat <- get_data(base_model)
  y_name <- get_y_name(base_model)
  check_in_data(y_name, dat, "data")
  Y_RAW <- dollar(dat, y_name)
  check_response(Y_RAW, 1)

  # Normalize response and format for Stan input
  normalizer <- create_scaling(Y_RAW, y_name) # create scaling
  y <- apply_scaling(normalizer, Y_RAW)
  N <- length(y)
  y <- array(y, dim = c(1, N))

  # Return
  list(
    to_stan = c(stan_prior, list(y = y)),
    y_scaling = normalizer
  )
}
