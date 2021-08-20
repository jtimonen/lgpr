# MAIN LEVEL ----------------------------------------------------------

#' Create a model
#'
#' @description See the
#' \href{https://jtimonen.github.io/lgpr-usage/articles/math.html}{Mathematical description of lgpr models}
#' vignette for more information about the connection between different options
#' and the created statistical model.
#' @export
#' @inheritParams create_model.common
#' @inheritParams create_model.likelihood
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
                         sample_f = "auto") {
  m <- create_model.common(formula, data, options, prior, prior_only, verbose)
  return(m)
}




# COMMON --------------------------------------------------------------

#' Create common Stan input needed for all models
#'
#' @inheritParams create_model.formula
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
#' @return An object of class \linkS4class{lgpmodel}
#' @family internal model creation functions
create_model.common <- function(formula, data, prior, options, prior_only,
                                verbose) {

  # Data, formula and common Stan inputs
  data <- convert_to_data_frame(data)
  lgp_formula <- create_model.formula(formula, data, verbose)
  standata_common <- standata_common(
    data, lgp_formula, options, prior, prior_only, verbose
  )

  # Variable names
  var_names <- list(
    y = lgp_formula@y_name,
    x = rownames(dollar(standata_common, "X")),
    z = rownames(dollar(standata_common, "Z"))
  )

  # Create the 'lgpmodel' object
  new("lgpmodel",
    model_formula = lgp_formula,
    data = data,
    parsed_input = standata_common,
    var_names = var_names,
    info = creation_info()
  )
}

# Misc info for created objects
creation_info <- function() {
  list(
    created = date(),
    lgpr_version = utils::packageVersion("lgpr")
  )
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


# APPROXIMATION -----------------------------------------------------------


#' Parse the given approximation options
#'
#' @param approx A named list with the following possible fields:
#' \itemize{
#'   \item \code{num_bf} Number of basis functions (0 = no approximation).
#'   \item \code{scale_bf} Scale of the domain to be used in basis
#'   function approximation. Has no effect if \code{num_bf = 0}.
#' }
#' If \code{approx} is \code{NULL}, default options are used. The defaults
#' are equivalent to
#' \code{options = list(
#'   num_bf = 0,
#'   scale_bf = 1.5
#' )
#' }.
#' @return a named list of parsed options
create_model.approx_options <- function(approx) {

  # Default options
  input <- approx
  opts <- list(num_bf = 0, scale_bf = 1.5)

  # Replace defaults if found from input
  for (opt_name in names(opts)) {
    if (opt_name %in% names(input)) {
      opts[[opt_name]] <- input[[opt_name]]
    }
  }

  # Validate
  num_bf <- dollar(opts, "num_bf")
  scale_bf <- dollar(opts, "scale_bf")
  check_non_negative_all(num_bf)
  check_non_negative_all(scale_bf)
  return(opts)
}

# PRIOR DEFINITIONS -------------------------------------------------------

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
