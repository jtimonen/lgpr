#' Create a model
#'
#' @description See the
#' \href{https://jtimonen.github.io/lgpr-usage/articles/math.html}{Mathematical description of lgpr models}
#' vignette for more information about the connection between different options
#' and the created statistical model.
#' @export
#' @inheritParams create_model.formula
#' @inheritParams create_model.prior
#' @inheritParams create_model.covs_and_comps
#' @inheritParams create_model.likelihood
#' @inheritParams create_model.options
#' @param prior_only Should likelihood be ignored? See also
#' \code{\link{sample_param_prior}} which can be used for any
#' \linkS4class{lgpmodel}, and whose runtime is independent of the number of
#' observations.
#' @family main functions
#' @return An object of class \linkS4class{lgpmodel}, containing the
#' Stan input created based on parsing the specified \code{formula},
#' \code{prior}, and other options.
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
  data <- convert_to_data_frame(data)
  lgp_formula <- create_model.formula(formula, data, verbose)
  cc_info <- create_model.covs_and_comps(data, lgp_formula, NA, verbose)
  stan_x <- dollar(cc_info, "to_stan")
  stan_opts <- create_model.options(options, verbose)

  # Parse response and likelihood
  y_info <- create_model.likelihood(
    data, likelihood, c_hat, num_trials, lgp_formula@y_name, sample_f, verbose
  )
  stan_y <- dollar(y_info, "to_stan")
  stan_input <- c(stan_x, stan_opts, stan_y)

  # Parse the prior
  prior_info <- create_model.prior(prior, stan_input, verbose)
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
  return(out)
}

#' Parse the given modeling options
#'
#' @inheritParams create_model.likelihood
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
create_model.options <- function(options, verbose) {
  # Default options
  log_progress("Parsing options...", verbose)
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


#' Prior definitions
#'
#' @param square is prior for a square-transformed parameter?
#' @name priors
#' @aliases normal, log_normal, gam, igam, uniform, student_t, bet
#' @return a named list
#' @description These use the same parametrizations as defined in the 'Stan'
#' documentation. See the docs for
#' \href{https://mc-stan.org/docs/2_24/functions-reference/gamma-distribution.html}{gamma} and
#' \href{https://mc-stan.org/docs/2_24/functions-reference/inverse-gamma-distribution.html}{inverse gamma} distributions.
#' @examples
#' # Log-normal prior
#' log_normal(mu = 1, sigma = 1)
#'
#' # Cauchy prior
#' student_t(nu = 1)
#'
#' # Exponential prior with rate = 0.1
#' gam(shape = 1, inv_scale = 0.1)
#'
#' # Create a similar priors as in LonGP (Cheng et al., 2019)
#' # Not recommended, because a lengthscale close to 0 is possible.
#' a <- log(1) - log(0.1)
#' log_normal(mu = 0, sigma = a / 2) # for continuous lengthscale
#' student_t(nu = 4) # for interaction lengthscale
#' igam(shape = 0.5, scale = 0.005, square = TRUE) # for sigma
NULL

#' @export
#' @rdname priors
uniform <- function(square = FALSE) {
  list(
    dist = "uniform",
    square = square
  )
}

#' @export
#' @rdname priors
#' @param mu mean
#' @param sigma standard deviation
normal <- function(mu, sigma, square = FALSE) {
  check_numeric(mu)
  check_positive(sigma)
  list(
    dist = "normal",
    square = square,
    mu = mu,
    sigma = sigma
  )
}

#' @export
#' @rdname priors
#' @param nu degrees of freedom
student_t <- function(nu, square = FALSE) {
  check_positive(nu)
  list(
    dist = "student-t",
    square = square,
    nu = nu
  )
}

#' @export
#' @rdname priors
#' @param shape shape parameter (alpha)
#' @param inv_scale inverse scale parameter (beta)
gam <- function(shape, inv_scale, square = FALSE) {
  check_positive(shape)
  check_positive(inv_scale)
  list(
    dist = "gamma",
    alpha = shape,
    beta = inv_scale,
    square = square
  )
}

#' @export
#' @rdname priors
#' @param shape shape parameter (alpha)
#' @param scale scale parameter (beta)
#' @family functions related to the inverse-gamma distribution
igam <- function(shape, scale, square = FALSE) {
  check_positive(shape)
  check_positive(scale)
  list(
    dist = "inv-gamma",
    alpha = shape,
    beta = scale,
    square = square
  )
}

#' @export
#' @rdname priors
#' @param mu mean
#' @param sigma standard deviation
log_normal <- function(mu, sigma, square = FALSE) {
  check_numeric(mu)
  check_positive(sigma)
  list(
    dist = "log-normal",
    square = square,
    mu = mu,
    sigma = sigma
  )
}

#' @export
#' @rdname priors
#' @param a shape parameter
#' @param b shape parameter
bet <- function(a, b) {
  check_positive(a)
  check_positive(b)
  list(
    dist = "beta",
    square = FALSE,
    alpha = a,
    beta = b
  )
}


#' Create a standardizing transform
#'
#' @param x variable measurements (might contain \code{NA} or \code{NaN})
#' @param name variable name
#' @return an object of class \linkS4class{lgpscaling}
#' @family variable scaling functions
create_scaling <- function(x, name) {
  check_length_geq(x, 2)
  loc <- mean(x, na.rm = TRUE)
  scale <- stats::sd(x, na.rm = TRUE)
  if (scale == 0) {
    msg <- paste0("the variable <", name, "> has zero variance!")
    stop(msg)
  }
  new("lgpscaling", loc = loc, scale = scale, var_name = name)
}


#' Apply variable scaling
#'
#' @param scaling an object of class \linkS4class{lgpscaling}
#' @param x object to which apply the scaling (numeric)
#' @param inverse whether scaling should be done in inverse direction
#' @return a similar object as \code{x}
#' @family variable scaling functions
apply_scaling <- function(scaling, x, inverse = FALSE) {
  loc <- scaling@loc
  scale <- scaling@scale
  if (inverse) {
    x <- scale * x + loc
  } else {
    x <- (x - loc) / scale
  }
  return(x)
}
