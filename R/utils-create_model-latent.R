# Create a latent GP model from base model
create_model.latent <- function(base_model, likelihood, prior,
                                c_hat, num_trials, approx, verbose) {
  stan_input <- list(approx = approx)
  si <- get_stan_input(base_model)

  # Base latent model
  si_add <- standata_latent(base_model, prior, c_hat, num_trials, verbose)
  si <- c(si, si_add)
  model <- new("LatentGPModel",
    base_model,
    parsed_input = stan_input,
    info = creation_info()
  )

  # Approximate model
  if (!is.null(approx)) {
    model <- create_model.latent_approx(model, approx)
  }
  return(model)
}

# Parse Stan data for latent model
standata_latent <- function(base_model, likelihood, prior,
                            c_hat, num_trials, approx, verbose) {
  return(list(moi = "joo"))
}

# Create a latent approximate GP model from latent GP model
create_model.latent_approx <- function(latent_model, approx) {
  si <- get_stan_input(latent_model)
  si_add <- approx # TODO: parse
  si <- c(si, si_add)
  new("LatentGPModelApprox",
    latent_model,
    parsed_input = si,
    info = creation_info()
  )
}
