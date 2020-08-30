#' Main functions of the `lgpr` package
#'
#' @export
#' @description
#' \itemize{
#'   \item The \code{lgp} function is a wrapper for both \code{create_model}
#'   and \code{sample_model}, which returns an object of class
#'   \linkS4class{lgpfit}.
#'   \item The \code{\link{create_model}} function creates an object of class
#'   \linkS4class{lgpmodel}.
#'   \item The \code{sample_model} function takes an \linkS4class{lgpmodel}
#'   object and fits the model using \code{rstan::sampling}.
#'   \item The \code{optimize_model} function takes an \linkS4class{lgpmodel}
#'   object and fits the model using \code{rstan::optimizing}.
#' }
#' @inheritParams create_model
#' @inheritParams sample_model
#' @inheritParams optimize_model
#' @name lgp
NULL

#' @rdname lgp
lgp <- function(formula,
                data,
                likelihood = "gaussian",
                prior = NULL,
                c_hat = NULL,
                num_trials = NULL,
                options = NULL,
                verbose = FALSE,
                ...) {

  # Create and fit the model
  model <- create_model(
    formula, data, likelihood, prior, c_hat, num_trials, options,
    verbose
  )
  stan_fit <- sample_model(model = model, ...)
  new("lgpfit", model = model, stan_fit = stan_fit)
}

#' @rdname lgp
#' @export
#' @param model An object of class \linkS4class{lgpmodel}.
#' @param ... Optional arguments passed to \code{rstan::sampling}, for example
#' \code{iter}, \code{chains} or \code{control}. See
#' \code{\link[rstan]{sampling}} for the possible arguments.
sample_model <- function(model, ...) {
  object <- stanmodels[[model@stan_model_name]]
  data <- model@stan_input
  rstan::sampling(object = object, data = data, check_data = TRUE, ...)
}

#' @rdname lgp
#' @export
#' @param model An object of class \linkS4class{lgpmodel}.
#' @param ... Optional arguments passed to \code{rstan::optimizing}, such as
#' \code{algorithm}. See \code{\link[rstan]{optimizing}} for the possible
#' arguments.
optimize_model <- function(model, ...) {
  object <- stanmodels[[model@stan_model_name]]
  data <- model@stan_input
  rstan::optimizing(object = object, data = data, check_data = TRUE, ...)
}

#' Create a model
#'
#' @inheritParams parse_formula
#' @inheritParams parse_response
#' @inheritParams parse_likelihood
#' @inheritParams parse_prior
#' @inheritParams parse_covs_and_comps
#' @inheritParams parse_options
#' @param verbose Should more verbose output be printed?
create_model <- function(formula,
                         data,
                         likelihood = "gaussian",
                         prior = NULL,
                         c_hat = NULL,
                         num_trials = NULL,
                         options = NULL,
                         verbose = FALSE) {

  # Parse the formula and options
  model_formula <- parse_formula(formula)
  list_opts <- parse_options(options)
  list_dopts <- list()

  # Parse response and likelihood
  parsed <- parse_response(data, likelihood, model_formula)
  list_y <- parsed$to_stan
  y_scaling <- parsed$scaling
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
  parsed <- parse_prior(prior, list_x, list_lh$obs_model)
  list_prior <- parsed$to_stan
  if (verbose) {
    cat(parsed$info)
  }

  # Other
  list_other <- list(is_verbose = as.numeric(verbose))

  # Create slots of the 'lgpmodel' object
  stan_input <- c(
    list_opts,
    list_dopts,
    list_y,
    list_lh,
    list_x,
    list_prior,
    list_other
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
    info = info,
    stan_model_name = "lgp"
  )
  return(out)
}
