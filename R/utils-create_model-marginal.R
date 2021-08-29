# Create a marginal GP model from base model
create_model.marginal <- function(base_model, prior, verbose) {
  si <- get_stan_input(base_model)
  parsed <- parse_obs_model.marginal(base_model)
  si_y <- dollar(parsed, "to_stan")
  si_prior <- parse_prior.marginal(base_model, prior, verbose)
  si <- c(si, si_y, si_prior)
  new("MarginalGPModel",
    base_model,
    y_scaling = dollar(parsed, "y_scaling"),
    info = creation_info()
  )
}

# Create additional Stan input
parse_prior.marginal <- function(base_model, prior, verbose) {
  filled <- fill_prior(prior, 0, "sigma")
  defaulting_info(filled, verbose)
  raw <- dollar(filled, "prior")
  parse_prior_pos(raw, 1, "sigma")
}

# Create additional Stan input
parse_obs_model.marginal <- function(base_model) {
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
  list(to_stan = list(y = y), y_scaling = normalizer)
}
