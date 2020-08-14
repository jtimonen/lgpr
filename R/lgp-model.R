#' Create an lgpmodel
#'
#' @export
#' @inheritParams parse_formula
#' @inheritParams parse_response
#' @inheritParams parse_likelihood
#' @inheritParams parse_prior
#' @inheritParams parse_data
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

  # Parse response
  parsed <- parse_response(data, likelihood, model_formula)
  list_y <- parsed$to_stan
  y_scaling <- parsed$scaling

  # Parse likelihood
  list_lh <- parse_likelihood(likelihood, c_hat, num_trials, list_y)

  # Parse covariates
  parsed <- parse_data(data, model_formula, id_variable)
  list_x <- parsed$to_stan
  x_cont_scalings <- parsed$x_cont_scalings
  x_cat_levels <- parsed$x_cat_levels

  # Group variable names
  var_names <- list(
    y = model_formula@y_name,
    x_cont = names(x_cont_scalings),
    x_cat = names(x_cat_levels)
  )

  # Parse components
  # parsed <- parse_components(model_formula, x_names)

  # Parse the prior
  # stan_prior <- parse_prior(prior)


  # Create slots of the 'lgpmodel' object
  stan_input <- c(
    list_opts,
    list_dopts,
    list_y,
    list_lh,
    list_x
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
#' @export
#' @param model bject of class \code{\link{lgpmodel}}.
#' @param type what types of covariates to include. Must be one of
#' {"all", "categorical", "continuous"}.
#' @return A string where names are separated by a comma.
covariate_names <- function(model, type = "all") {
  allowed <- c("continuous", "categorical", "all")
  nam1 <- model@var_names$x_cont
  nam2 <- model@var_names$x_cat
  nam3 <- c(nam1, nam2)
  names <- list(nam1, nam2, nam3)
  idx <- argument_check(arg = type, allowed = allowed)
  out <- names[[idx]]
  paste(out, collapse = ", ")
}
