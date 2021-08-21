#' Create a model
#'
#' @description See the
#' \href{https://jtimonen.github.io/lgpr-usage/articles/math.html}{Mathematical description of lgpr models}
#' vignette for more information about the connection between different options
#' and the created statistical model.
#' @export
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
#' @param data A \code{data.frame} where each column corresponds to one
#' variable, and each row is one observation. Continuous covariates and the
#' response variable must have type \code{"numeric"} and categorical covariates
#' must have type \code{"factor"}. Missing values should be indicated with
#' \code{NaN} or \code{NA}. The response variable cannot contain missing
#' values. Column names should not contain trailing or leading underscores.
#' @param likelihood Determines the observation model. Must be either
#' \code{"gaussian"} (default), \code{"poisson"}, \code{"nb"} (negative
#' binomial), \code{"binomial"} or \code{"bb"} (beta binomial).
#' @param prior A named list, defining the prior distribution of model
#' (hyper)parameters. See the "Defining priors" section below
#' (\code{\link{lgp}}).
#' @param c_hat Constant added to the zero-mean GP before mapping through
#' inverse link function. This should only be given if \code{sample_f} is
#' \code{TRUE}, otherwise this is zero (as response is normalized to have
#' zero mean and unit variance). If \code{sample_f}
#' is \code{TRUE}, the given \code{c_hat} can be a vector of length
#' \code{dim(data)[1]}, or a real number defining a constant GP mean. If not
#' specified and \code{sample_f} is \code{TRUE}, \code{c_hat} is set to
#'  \itemize{
#'    \item \code{c_hat = mean(y)}, if \code{likelihood} is \code{"gaussian"},
#'    \item \code{c_hat = } \code{log(mean(y))} if \code{likelihood} is
#'    \code{"poisson"} or \code{"nb"},
#'    \item \code{c_hat = } \code{log(p/(1-p))}, where
#'    \code{p = mean(y/num_trials)} if \code{likelihood} is \code{"binomial"}
#'    or \code{"bb"},
#'  }
#' where \code{y} denotes the response variable measurements.
#' @param num_trials This argument (number of trials) is only needed when
#' likelihood is \code{"binomial"} or \code{"bb"}. Must have length one or
#' equal to the number of data points. Setting \code{num_trials=1} and
#' \code{likelihood="binomial"} corresponds to Bernoulli observation model.
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
#' @param prior_only Should likelihood be ignored? See also
#' \code{\link{sample_param_prior}} which can be used for any
#' \linkS4class{lgpmodel}, and whose runtime is independent of the number of
#' observations.
#' @param verbose Should more informational messages be printed?
#' @param sample_f Determines if the latent function values are sampled.
#' Can be either \code{TRUE}, \code{FALSE} or \code{"auto"} (default). In
#' the latter case, this is set to \code{TRUE} if likelihood is not
#' \code{"gaussian"}, or if a basis function approximation is specified.
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
                         sample_f = "auto",
                         approx = NULL) {
  bm <- create_model.base(formula, data, options, prior, prior_only, verbose)
  sample_f <- determine_sample_f(approx, likelihood, sample_f)
  if (!sample_f) {
    m <- create_model.marginal(bm, prior, verbose)
  } else {
    m <- create_model.latent(
      bm, likelihood, prior, c_hat, num_trials,
      approx, verbose
    )
  }
  return(m)
}

# Convert sample_f input to TRUE or FALSE
determine_sample_f <- function(approx, likelihood, sample_f_input) {
  if (is.logical(sample_f_input)) {
    val <- sample_f_input
  } else {
    if (sample_f_input == "auto") {
      val <- (likelihood != "gaussian") || !is.null(approx)
    } else {
      stop("unrecognized argument sample_f=", sample_f_input)
    }
  }
  return(val)
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
#'
#' @section Prior distributions:
#' \itemize{
#'    \item \code{normal} - normal distribution
#'    \item \code{log_normal} - log-normal distribution
#'    \item \code{student_t} - Student-\emph{t} distribution
#'    \item \code{gam} - gamma distribution
#'    \item \code{igam} - inverse gamma distribution
#'    \item \code{bet} - beta distribution
#' }
#' The beta prior can be specified only for parameters that are restricted
#' on the unit interval (i.e for the \code{beta} and \code{gamma} parameters).
#' Other priors can be used for any other parameters. Note that for example
#' if the parameter is restricted to  be positive, then for example the
#' normal prior is actually a half-normal prior.
#'
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
