
#' Create an lgp model
#'
#' @export
#' @description Creates an object of class \code{lgpmodel}
#' @inheritParams parse_formula
#' @inheritParams parse_data
#' @inheritParams parse_likelihood
#' @inheritParams parse_prior
#' @inheritParams parse_options
#' @inheritParams parse_disease_options
#' @param id_variable Name of the unique subject identifier variable
#' (default = \code{"id"}).
#' @return An object of class \code{\link{lgpmodel}}.
#' @seealso For fitting the model, see \code{\link{lgp_fit}}.
lgp_model <- function(formula,
                      data,
                      likelihood = "gaussian",
                      prior = prior_default(),
                      c_hat = NULL,
                      num_trials = NULL,
                      id_variable = "id",
                      options = NULL,
                      disease_options = NULL) {

  # Parse the formula and options
  model_formula <- parse_formula(formula)
  stan_opts <- parse_options(options)
  stan_dopts <- parse_disease_options(disease_options)
  stan_likelihood <- parse_likelihood(likelihood, c_hat, num_trials)
  parsed_data <- parse_data(data, model_formula)
  stan_prior <- parse_prior(prior)

  # Create the 'lgpmodel' object
  stan_input <- c(parsed_data$stan_data,
                  stan_likelihood,
                  stan_prior,
                  stan_opts, 
                  stan_dopts)
  info <- paste("lgpmodel object, created", date())
  out <- new("lgpmodel",
    model_formula = model_formula,
    stan_input = stan_input,
    scalings = parsed_data$scalings,
    info = info
  )
  return(out)
}
