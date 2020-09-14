#' Posterior predictive distribution at test points
#'
#' @export
#' @param fit An object of class \linkS4class{lgpfit}.
#' @param X_pred A data frame of points where predictions are computed.
#' The function \code{\link{new_x}} can help in creating it.
#' @param STREAM external pointer
#' @return A list.
#' @family GP posterior computation functions
posterior_predict <- function(fit, X_pred, STREAM = get_stream()) {
  model <- object_to_model(fit)
  check_type(X_pred, "data.frame")
  obs_model <- get_obs_model(fit)
  if (obs_model != "gaussian") stop("observation model must be gaussian!")

  pred <- format_X_pred(model, X_pred)
  theta <- get_draws_kernel_pars(fit)
  sigma <- get_draws(fit, pars = "sigma", stack_chains = TRUE)
  kernels <- compute_kernel_matrices(model, pred, theta, STREAM)

  si <- get_stan_input(model)
  delta <- dollar(si, "delta")
  compute_gp_posteriors(kernels, y, delta, sigma, STREAM)
}

#' Prediction points to Stan format
#'
#' @param model An object of class \linkS4class{lgpmodel}.
#' @inheritParams posterior_predict
#' @return A list.
#' @family GP posterior computation functions
format_X_pred <- function(model, X_pred) {
  x_names <- rhs_variables(model@model_formula@terms)
  x_names <- unique(x_names)
  for (name in x_names) check_in_data(name, X_pred)
  stan_data_covariates(X_pred, x_names)
}

#' Compute kernel matrices
#'
#' @param model An object of class \linkS4class{lgpmodel}.
#' @param pred object returned by \code{\link{format_X_pred}}
#' @param theta a list with length \code{num_draws}
#' @inheritParams posterior_predict
#' @return A list.
#' @family GP posterior computation functions
compute_kernel_matrices <- function(model, pred, theta, STREAM) {
  "TODO"
}

#' Compute GP posteriors
#'
#' @param kernels a list with length \code{num_draws}
#' @param y response variable vector of length \code{num_obs}
#' @param delta jitter to ensure positive definite matrices
#' @param sigma a vector with length \code{num_draws}
#' @inheritParams posterior_predict
#' @return A list.
#' @family GP posterior computation functions
compute_gp_posteriors <- function(kernels, y, delta, sigma, STREAM) {
  # gp_posterior(KX, KX_s, KX_ss, y, delta, sigma, STREAM)
  "TODO"
}

#' Extract kernel parameter draws
#'
#' @inheritParams posterior_predict
#' @family GP posterior computation functions
#' @return a list with names  "alpha", "ell", "wrp", "beta" and "teff", where
#' each element is an array with number of rows equal to number of posterior
#' samples
get_draws_kernel_pars <- function(fit) {
  num_draws <- get_num_draws(fit)
  na_array <- array(NA, dim = c(num_draws, 1))
  par_names <- c("alpha", "ell", "wrp", "beta", "teff")
  out <- list()
  for (pn in par_names) {
    draws <- get_draws_catch(fit, stack_chains = TRUE, pars = pn)
    out[[pn]] <- if (is.null(draws)) na_array else draws
  }
  return(out)
}
