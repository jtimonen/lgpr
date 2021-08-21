# Create a marginal GP model from base model
create_model.marginal <- function(base_model, prior, verbose) {
  si <- base_model@parsed_input
  si_add <- standata_marginal(base_model, prior, verbose)
  si <- c(si, si_add)
  new("MarginalGPModel",
    base_model,
    parsed_input = si,
    info = creation_info()
  )
}

# Create additional Stan input
standata_marginal <- function(base_model, prior, verbose) {
  fp <- fill_prior(prior, 0, "sigma")
  print(fp)
  return(fp)
}
