# Create a marginal GP model from base model
create_model.latent <- function(base_model, likelihood, prior, approx, verbose) {
  stan_input <- list(approx = approx)
  new("LatentGPModel",
    base_model,
    stan_input = stan_input,
    info = creation_info()
  )
}
