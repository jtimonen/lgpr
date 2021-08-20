# Create a marginal GP model from base model
create_model.marginal <- function(base_model, likelihood, prior, verbose) {
  stan_input <- list(joo = "moi")
  new("MarginalGPModel",
    base_model,
    stan_input = stan_input,
    info = creation_info()
  )
}
