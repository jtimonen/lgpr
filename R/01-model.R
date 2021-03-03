#' Create a model
#'
#' @export
#' @inheritParams parse_formula
#' @inheritParams parse_prior
#' @inheritParams parse_covs_and_comps
#' @inheritParams parse_y
#' @inheritParams parse_options
#' @param prior_only Should sampling be done in prior mode,
#' where likelihood is ignored?.
#' @family main functions
create_model <- function(formula,
                         data,
                         likelihood = "gaussian",
                         prior = NULL,
                         c_hat = NULL,
                         num_trials = NULL,
                         options = NULL,
                         prior_only = FALSE,
                         verbose = FALSE,
                         sample_f = !(likelihood == "gaussian")) {

  # Parse common parts (formula, covariates, components, options)
  lgp_formula <- parse_formula(formula, data, verbose)
  cc_info <- parse_covs_and_comps(data, lgp_formula, NA, verbose)
  stan_x <- dollar(cc_info, "to_stan")
  stan_opts <- parse_options(options, verbose)

  # Parse response and likelihood
  y_info <- parse_y(
    data, likelihood, c_hat, num_trials, lgp_formula@y_name, sample_f, verbose
  )
  stan_y <- dollar(y_info, "to_stan")
  stan_input <- c(stan_x, stan_opts, stan_y)

  # Parse the prior
  prior_info <- parse_prior(prior, stan_input, verbose)
  full_prior <- dollar(prior_info, "raw")
  stan_input <- c(stan_input, dollar(prior_info, "to_stan"))

  # Binary option switches
  stan_switches <- list(
    is_verbose = as.numeric(verbose),
    is_likelihood_skipped = as.numeric(prior_only)
  )
  stan_input <- c(stan_input, stan_switches)

  # Variable info
  var_names <- dollar(cc_info, "var_names")
  var_scalings <- list(
    y = dollar(y_info, "scaling"),
    x_cont = dollar(cc_info, "x_cont_scalings")
  )
  var_info <- list(x_cat_levels = dollar(cc_info, "x_cat_levels"))

  # Misc model info
  info <- list(
    created = date(),
    lgpr_version = utils::packageVersion("lgpr"),
    caseid_map = dollar(cc_info, "caseid_map")
  )

  # Create the 'lgpmodel' object
  if (verbose) cat("Creating lgpmodel object...\n")
  out <- new("lgpmodel",
    model_formula = lgp_formula,
    data = data,
    var_names = var_names,
    var_info = var_info,
    var_scalings = var_scalings,
    stan_input = stan_input,
    info = info,
    sample_f = sample_f,
    full_prior = full_prior
  )
  if (verbose) cat("Done.\n")
  return(out)
}

#' Parse the given modeling options
#'
#' @inheritParams parse_y
#' @param options A named list with the following possible fields:
#' \itemize{
#'   \item \code{delta} Amount of added jitter to ensure positive definite
#'   covariance matrices.
#'   \item \code{vm_params} Variance mask function parameters (numeric
#'   vector of length 2).
#' }
#' If \code{options} is \code{NULL}, default options are used. The defaults
#' are equivalent to
#' \code{options = list(delta = 1e-8,  vm_params = c(0.025, 1))}.
#' @return a named list of parsed options
parse_options <- function(options, verbose) {

  # Default options
  if (verbose) cat("Parsing options...\n")
  input <- options
  opts <- list(delta = 1e-8, vm_params = c(0.025, 1))

  # Replace defaults if found from input
  for (opt_name in names(opts)) {
    if (opt_name %in% names(input)) {
      opts[[opt_name]] <- input[[opt_name]]
    }
  }

  # Validate and format for Stan input
  delta <- dollar(opts, "delta")
  check_positive(delta)
  vm_params <- dollar(opts, "vm_params")
  check_length(vm_params, 2)
  check_positive_all(vm_params)
  check_all_leq(vm_params, c(1, 1))
  list(delta = delta, vm_params = vm_params)
}
