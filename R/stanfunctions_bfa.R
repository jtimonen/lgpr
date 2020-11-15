
#' The m'th eigenfunction of the Dirichlet boundary value problem
#'
#' @description This is a wrapper for \code{STAN_phi}, with some input
#' validation added.
#' @param x a numeric vector
#' @param m index of the eigenfunction (non-negative integer)
#' @param L domain width (positive real number)
#' @inheritParams cpp_kernel_const_all
#' @family basis function approximation functions
#' @return with length equal to that of \code{x}
cpp_phi <- function(x, m, L, STREAM = get_stream()) {
  check_positive(m)
  check_positive(L)
  STAN_phi(x, m, L, STREAM)
}

#' The m'th eigenvalue of the Dirichlet boundary value problem
#'
#' @description This is a wrapper for \code{STAN_lambda}, with some input
#' validation added.
#' @param m index of the eigenfunction (non-negative integer)
#' @param L domain width (positive real number)
#' @inheritParams cpp_kernel_const_all
#' @family basis function approximation functions
#' @return a number
cpp_lambda <- function(m, L, STREAM = get_stream()) {
  check_positive(m)
  check_positive(L)
  STAN_lambda(m, L, STREAM)
}

#' Spectral density function of the exponentiated quadratic kernel
#'
#' @description This is a wrapper for \code{STAN_spd_eq}, with some input
#' validation added.
#' @param w the frequency, real number
#' @param ell lengthscale of the kernel
#' @inheritParams cpp_kernel_const_all
#' @family basis function approximation functions
#' @return a real number
cpp_spd_eq <- function(w, ell, STREAM = get_stream()) {
  check_positive(ell)
  STAN_spd_eq(w, ell, STREAM)
}

#' Multivariate normal density using basis function approximation
#'
#' @description This is a wrapper for \code{STAN_multi_normal_bfa_logpdf},
#' with some input validation added.
#' @param y a vector of length \code{n}
#' @param V a matrix of size \code{n} x \code{R*M}
#' @param D_diag a vector of length \code{R*M}
#' @param sigma noise standard deviation, a positive real number
#' @inheritParams cpp_kernel_const_all
#' @family basis function approximation functions
#' @return a real number
cpp_multi_normal_bfa_logpdf <- function(y, V, D_diag, sigma,
                                        STREAM = get_stream()) {
  check_positive(sigma)
  STAN_multi_normal_bfa_logpdf(y, V, D_diag, sigma, STREAM)
}
