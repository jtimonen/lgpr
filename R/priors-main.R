#' Prior definitions
#'
#' @param square is prior for a square-transformed parameter?
#' @name priors
#' @aliases normal, log_normal, gamma, inv_gamma, uniform, student_t
#' @return a named list
#' @description These use the same parametrizations as defined in the Stan
#' documentation. See the docs for
#' \href{https://mc-stan.org/docs/2_24/functions-reference/
#' gamma-distribution.html}{gamma} and
#' \href{https://mc-stan.org/docs/2_24/functions-reference/
#' inverse-gamma-distribution.html}{inverse gamma} distributions.
#' @examples
#' # Log-normal prior
#' log_normal(mu = 1, sigma = 1)
#'
#' # Cauchy prior
#' student_t(nu = 1)
#'
#' # Exponential prior with rate = 0.1
#' gamma(shape = 1, inv_scale = 0.1)
NULL

#' @export
#' @rdname priors
uniform <- function(square = FALSE) {
  warning("specifying uniform prior")
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
  check_numeric(sigma, require_positive = TRUE)
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
  check_numeric(nu, require_positive = TRUE)
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
gamma <- function(shape, inv_scale, square = FALSE) {
  check_numeric(shape, require_positive = TRUE)
  check_numeric(inv_scale, require_positive = TRUE)
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
inv_gamma <- function(shape, scale, square = FALSE) {
  check_numeric(shape, require_positive = TRUE)
  check_numeric(scale, require_positive = TRUE)
  list(
    dist = "inverse-gamma",
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
  check_numeric(sigma, require_positive = TRUE)
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
beta <- function(a, b) {
  check_numeric(a, require_positive = TRUE)
  check_numeric(b, require_positive = TRUE)
  list(
    dist = "beta",
    square = FALSE,
    alpha = a,
    beta = b
  )
}
