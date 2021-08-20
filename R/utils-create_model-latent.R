#' Create common Stan input needed for all models
#'
#' @param formula The model formula, where
#' \itemize{
#'   \item it must contain exactly one tilde (\code{~}), with response
#'   variable on the left-hand side and model terms on the right-hand side
#'   \item terms are be separated by a plus (\code{+}) sign
#'   \item all variables appearing in \code{formula} must be
#'   found in \code{data}
#' }
#' See the "Model formula syntax" section below (\code{\link{lgp}}) for
#' instructions on how to specify the model terms.
#' @param data A data frame.
#' @param options A named list with the following possible fields:
#' \itemize{
#'   \item \code{delta} Amount of added jitter to ensure positive definite
#'   covariance matrices.
#'   \item \code{vm_params} Variance mask function parameters (numeric
#'   vector of length 2).
#' }
#' If \code{options} is \code{NULL}, default options are used. The defaults
#' are equivalent to
#' \code{options = list(delta = 1e-8, vm_params = c(0.025, 1))}.
#' @param prior A named list, defining the prior distribution of model
#' (hyper)parameters. See the "Defining priors" section below
#' (\code{\link{lgp}}).
#' @param prior_only Should likelihood be ignored? See also
#' \code{\link{sample_param_prior}} which can be used for any
#' \linkS4class{lgpmodel}, and whose runtime is independent of the number of
#' observations.
#' @return An object of class \linkS4class{lgpmodel}.
#' @family internal model creation functions
create_model_base <- function(formula, data, prior, options, prior_only,
                                verbose) {

  # Data, formula and common Stan inputs
  data <- convert_to_data_frame(data)
  lgp_formula <- create_lgpformula(formula, data, verbose)
  parsed_input <- standata_common(
    data, lgp_formula, options, prior, prior_only, verbose
  )

  # Variable names
  var_names <- list(
    y = lgp_formula@y_name,
    x = rownames(dollar(parsed_input, "X")),
    z = rownames(dollar(parsed_input, "Z"))
  )

  # Create the 'lgpmodel' object
  new("lgpmodel",
    model_formula = lgp_formula,
    data = data,
    parsed_input = parsed_input,
    var_names = var_names,
    info = creation_info()
  )
}

# Checks if formula is in advanced format and translates if not
create_lgpformula <- function(formula, data, verbose = FALSE) {
  advanced <- is_advanced_formula(formula)
  if (!advanced) formula <- formula_to_advanced(formula, data, verbose)
  fp <- as.character(formula)
  formula_str <- paste(fp[2], fp[1], fp[3])
  log_info(paste0("Formula interpreted as: ", formula_str), verbose)
  parse_formula_advanced(formula)
}

# Create common Stan input needed for all models
standata_common <- function(data, lgp_formula, opts, prior, prior_only, vrb) {
  opts <- standata_common_options(opts, prior_only)
  covs <- standata_covariates(data, lgp_formula)
  comps <- standata_components(lgp_formula, covs)
  expanding <- standata_expanding(covs, comps)
  si <- c(opts, covs, comps, expanding)
  pri <- standata_common_prior(prior, si, vrb)
  lst <- c(si, pri)
  return(lst)
}

# Misc info for created objects
creation_info <- function() {
  list(
    created = date(),
    lgpr_version = utils::packageVersion("lgpr")
  )
}
