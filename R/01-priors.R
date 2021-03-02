#' Prior definitions
#'
#' @param square is prior for a square-transformed parameter?
#' @name priors
#' @aliases normal, log_normal, gam, igam, uniform, student_t, bet
#' @return a named list
#' @description These use the same parametrizations as defined in the Stan
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
