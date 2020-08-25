#' Create a model
#'
#' @inheritParams parse_formula
#' @inheritParams parse_response
#' @inheritParams parse_likelihood
#' @inheritParams parse_prior
#' @inheritParams parse_covs_and_comps
#' @inheritParams parse_options
#' @param verbose Should more verbose output be printed?
lgp_model <- function(formula,
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
    info = info
  )
  return(out)
}


#' Get covariate names for a model
#'
#' @param model an object of class \linkS4class{lgpmodel}
#' @param type what types of covariates to include. Must be either
#' "all", "categorical" or "continuous".
#' @return A string where names are separated by a comma.
get_covariate_names <- function(model, type = "all") {
  allowed <- c("continuous", "categorical", "all")
  nam1 <- model@var_names$x_cont
  nam2 <- model@var_names$x_cat
  nam3 <- c(nam1, nam2)
  names <- list(nam1, nam2, nam3)
  idx <- argument_check(arg = type, allowed = allowed)
  out <- names[[idx]]
  paste(out, collapse = ", ")
}

#' Get integer encoding of model components
#'
#' @param model an object of class \linkS4class{lgpmodel}
#' @return A data frame.
get_integer_encoding <- function(model) {
  data.frame(model@stan_input$components)
}
