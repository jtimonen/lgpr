#' Create a model
#'
#' @description See the
#' \href{https://jtimonen.github.io/lgpr-usage/articles/math.html}{Mathematical description of lgpr models}
#' vignette for more information about the connection between different options
#' and the created statistical model.
#' @export
#' @inheritParams create_model.base
#' @inheritParams create_model.base
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
                         sample_f = "auto",
                         approx = NULL) {
  bm <- create_model.base(formula, data, options, prior, prior_only, verbose)
  sample_f <- determine_sample_f(approx, likelihood, sample_f)
  if (!sample_f) {
    m <- create_model.marginal(bm, prior, verbose)
  } else {
    m <- create_model.latent(bm, likelihood, prior, approx, verbose)
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
