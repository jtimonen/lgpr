#' Main functions of the `lgpr` package
#'
#' @name lgp
#' @description
#' \itemize{
#'   \item The \code{lgp_model} function creates an object of class
#'   \linkS4class{lgpmodel}.
#'   \item The \code{lgp_fit} function takes an \linkS4class{lgpmodel} object
#'   and fits the model using \code{rstan::sampling}.
#'   \item The \code{lgp} function is a wrapper for both \code{lgp_model}
#'   and \code{lgp_fit}.
#' }
#' @aliases lgp_model, lgp_fit
NULL

#' @export
#' @rdname lgp
lgp <- function(formula,
                data,
                likelihood = "gaussian",
                prior = NULL,
                c_hat = NULL,
                num_trials = NULL,
                options = NULL,
                threshold = 0.95,
                relevance_method = "f_mean",
                ...) {

  # Create and fit the model
  model <- lgp_model(
    formula, data, likelihood, prior, c_hat, num_trials, options
  )
  fit <- lgp_fit(model = model, ...)
  list(model = model, fit = fit)
}

#' @export
#' @rdname lgp
#' @inheritParams parse_formula
#' @inheritParams parse_response
#' @inheritParams parse_likelihood
#' @inheritParams parse_prior
#' @inheritParams parse_covs_and_comps
#' @inheritParams parse_options
lgp_model <- function(formula,
                      data,
                      likelihood = "gaussian",
                      prior = NULL,
                      c_hat = NULL,
                      num_trials = NULL,
                      options = NULL) {

  # Parse the formula and options
  model_formula <- parse_formula(formula)
  list_opts <- parse_options(options)
  list_dopts <- list()

  # Parse response
  parsed <- parse_response(data, likelihood, model_formula)
  list_y <- parsed$to_stan
  y_scaling <- parsed$scaling

  # Parse likelihood
  list_lh <- parse_likelihood(likelihood, c_hat, num_trials, list_y)

  # Parse covariates and components
  parsed <- parse_covs_and_comps(data, model_formula)
  list_x <- parsed$to_stan
  x_cont_scalings <- parsed$x_cont_scalings
  x_cat_levels <- parsed$x_cat_levels

  # Group variable names
  var_names <- list(
    y = model_formula@y_name,
    x_cont = names(x_cont_scalings),
    x_cat = names(x_cat_levels)
  )

  # Parse the prior
  list_prior <- parse_prior(prior, list_x, list_lh$obs_model)

  # Create slots of the 'lgpmodel' object
  stan_input <- c(
    list_opts,
    list_dopts,
    list_y,
    list_lh,
    list_x,
    list_prior
  )
  info <- list(created = date(), pkg_desc = get_pkg_description())
  var_scalings <- list(y = y_scaling, x_cont = x_cont_scalings)
  var_info <- list(x_cat_levels = x_cat_levels)

  # Create the 'lgpmodel' object
  out <- new("lgpmodel",
    model_formula = model_formula,
    var_names = var_names,
    var_info = var_info,
    var_scalings = var_scalings,
    stan_input = stan_input,
    info = info
  )
  return(out)
}

#' @export
#' @rdname lgp
#' @param model An object of class \linkS4class{lgpmodel}.
#' @param ... Optional arguments passed to \code{rstan::sampling}, for example
#' \code{iter}, \code{chains} or \code{control}. See
#' \code{\link[rstan]{sampling}} for the possible arguments.
#' @param ... Additional arguments that are passed to \code{rstan::sampling}.
lgp_fit <- function(model, ...) {
  rstan::sampling(object = stanmodels[["lgp"]], data = model@stan_input, ...)
}
