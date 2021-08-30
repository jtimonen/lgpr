# Create a marginal GP model from base model
create_model.marginal <- function(base_model, prior, verbose) {
  obs_model <- parse_obs_model.marginal(base_model)
  pri <- parse_prior.marginal(base_model, prior, verbose)
  new("MarginalGPModel",
    base_model,
    y = dollar(obs_model, "y"),
    y_scaling = dollar(obs_model, "y_scaling"),
    prior_sigma = dollar(pri, "prior_sigma"),
    hyper_sigma = dollar(pri, "hyper_sigma"),
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
  list(y = y, y_scaling = normalizer)
}
