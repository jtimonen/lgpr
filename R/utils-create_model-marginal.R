# Create a marginal GP model from base model
create_model.marginal <- function(base_model, likelihood, prior, verbose) {

  new("MarginalGPModel",
    info = creation_info()
  )
}
