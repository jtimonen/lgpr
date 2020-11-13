
#' The m'th eigenfunction of the Dirichlet boundary value problem
#'
#' @description This is a wrapper for \code{STAN_bfa_phi}, with some input
#' validation added.
#' @param x a numeric vector
#' @param m index of the eigenfunction (non-negative integer)
#' @param L domain width (positive real number)
#' @inheritParams cpp_kernel_const_all
#' @family basisfunction approximation functions
#' @return with length equal to that of \code{x}
cpp_bfa_phi <- function(x, m, L, STREAM = get_stream()) {
  check_positive(m)
  check_positive(L)
  STAN_bfa_phi(x, m, L, STREAM)
}

#' The m'th eigenvalue of the Dirichlet boundary value problem
#'
#' @description This is a wrapper for \code{STAN_bfa_lambda}, with some input
#' validation added.
#' @param m index of the eigenfunction (non-negative integer)
#' @param L domain width (positive real number)
#' @inheritParams cpp_kernel_const_all
#' @family basisfunction approximation functions
#' @return a number
cpp_bfa_lambda <- function(m, L, STREAM = get_stream()) {
  check_positive(m)
  check_positive(L)
  STAN_bfa_lambda(m, L, STREAM)
}

#' Spectral density function of the exponentiated quadratic kernel
#'
#' @description This is a wrapper for \code{STAN_spd_eq}, with some input
#' validation added.
#' @param w the frequency, real number
#' @param alpha marginal std of the kernel
#' @param ell lengthscale of the kernel
#' @inheritParams cpp_kernel_const_all
#' @family basisfunction approximation functions
#' @return a real number
cpp_spd_eq <- function(w, alpha, ell, STREAM = get_stream()) {
  check_positive(alpha)
  check_positive(ell)
  STAN_spd_eq(w, alpha, ell, STREAM)
}

#' Multivariate normal density using basis function approximation
#'
#' @description This is a wrapper for \code{STAN_bfa_multi_normal_lpdf},
#' with some input validation added.
#' @param y a vector of length \code{n}
#' @param V a matrix of size \code{n} x \code{R*M}
#' @param D_diag a vector of length \code{R*M}
#' @param sigma noise standard deviation, a positive real number
#' @inheritParams cpp_kernel_const_all
#' @family basisfunction approximation functions
#' @return a real number
cpp_bfa_multi_normal_lpdf <- function(y, V, D_diag, sigma,
                                      STREAM = get_stream()) {
  check_positive(sigma)
  STAN_bfa_multi_normal_lpdf(y, V, D_diag, sigma, STREAM)
}
