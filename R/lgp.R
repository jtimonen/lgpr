#' The main function of the `lgpr` package
#'
#' @export
#' @description This is a wrapper for both \code{\link{lgp_model}}
#'  and \code{\link{lgp_fit}}. It first creates an
#'  \linkS4class{lgpmodel} object and then fits the model,
#'  finally returning an \linkS4class{lgpfit} object.
#' @inheritParams lgp_model
#' @inheritParams lgp_fit
#' @return An object of class \code{lgpfit}.
lgp <- function(formula,
                data,
                likelihood = "gaussian",
                prior = prior_default(),
                c_hat = NULL,
                num_trials = NULL,
                id_variable = "id",
                options = NULL,
                disease_options = NULL,
                
                parallel = FALSE,
                threshold = 0.95,
                relevance_method = "f_mean",
                ...) {

  # Create the model
  model <- lgp_model(formula, data, likelihood, prior, c_hat, num_trials,
                     id_variable, options, disease_options)

  # Fit the model
  #fit <- lgp_fit(
  #  model = model,
  #  threshold = threshold,
  #  parallel = parallel,
  #  relevance_method = relevance_method,
  #  ...
  #)
  fit <- model

  # Return the lgpfit object
  return(fit)
}


#' Create an lgpmodel
#'
#' @export
#' @inheritParams parse_formula
#' @inheritParams parse_response
#' @inheritParams parse_likelihood
#' @inheritParams parse_prior
#' @inheritParams parse_covariates
#' @inheritParams parse_options
#' @inheritParams parse_disease_options
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
  list_opts <- parse_options(options)
  list_dopts <- parse_disease_options(disease_options)
  
  # Parse response and likelihood
  parsed <- parse_response(data, likelihood, model_formula)
  list_y <- parsed$y_to_stan
  y_scl <- parsed$y_scaling
  list_lh <- parse_likelihood(likelihood, c_hat, num_trials, list_y)
  
  # Parse covariates and components
  # p_covariates <- parse_covariates(data, model_formula)
  
  # Parse the prior
  # stan_prior <- parse_prior(prior)
  
  # Create the 'lgpmodel' object
  stan_input <- c(
    list_opts,
    list_dopts,
    list_y,
    list_lh
  )
  scalings <- list(response = y_scl)
  
  info <- list(created = date(), pkg_desc = get_pkg_description())
  out <- new("lgpmodel",
             model_formula = model_formula,
             stan_input = stan_input,
             scalings = scalings,
             info = info
  )
  return(out)
} 
