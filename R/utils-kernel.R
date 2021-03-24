#' Compute a kernel matrix (covariance matrix)
#'
#' @param x1 vector of length \eqn{n}
#' @param x2 vector of length \eqn{m}
#' @param alpha marginal std (default = 1)
#' @param ell lengthscale
#' @return A matrix of size \eqn{n} x \eqn{m}.
#' @name kernel
NULL

#' @describeIn kernel Use the squared exponential kernel.
kernel_se <- function(x1, x2, alpha = 1.0, ell) {
  check_positive(ell)
  check_non_negative(alpha)
  n1 <- length(x1)
  n2 <- length(x2)
  X1 <- matrix(rep(x1, each = n2), n1, n2, byrow = T)
  X2 <- matrix(rep(x2, n1), n1, n2, byrow = T)
  K <- alpha^2 * exp(-0.5 * (X1 - X2)^2 / ell^2)
  return(K)
}

#' @describeIn kernel Use the nonstationary kernel (input warping + squared
#' exponential).
#' @param a steepness of the warping function rise
kernel_ns <- function(x1, x2, alpha = 1.0, ell, a) {
  nan_replace <- 0
  check_positive(a)
  x1[is.nan(x1)] <- nan_replace
  x2[is.nan(x2)] <- nan_replace
  w1 <- warp_input(x1, a)
  w2 <- warp_input(x2, a)
  K <- kernel_se(w1, w2, alpha, ell)
  return(K)
}


#' @describeIn kernel Uses the zero-sum kernel. Here, \code{x1} and
#' \code{x2} must be integer vectors (integers denoting different categories).
#' Returns a binary matrix.
#' @param M number of categories
kernel_zerosum <- function(x1, x2, M, alpha = 1.0) {
  check_non_negative(alpha)
  n1 <- length(x1)
  n2 <- length(x2)
  K <- matrix(0, n1, n2)
  for (i in 1:n1) {
    for (j in 1:n2) {
      if (x1[i] == x2[j]) {
        K[i, j] <- 1.0
      } else {
        K[i, j] <- -1.0 / (M - 1)
      }
    }
  }
  return(alpha^2 * K)
}

#' @describeIn kernel Uses the binary (mask) kernel. Here, \code{x1} and
#' \code{x2} must be integer vectors (integers denoting different categories).
#' Returns a binary matrix.
#' @param pos_class the positive class label (integer)
kernel_bin <- function(x1, x2 = NULL, alpha = 1.0, pos_class = 1) {
  check_non_negative(alpha)
  n1 <- length(x1)
  n2 <- length(x2)
  X1 <- matrix(rep(x1, each = n2), n1, n2, byrow = T)
  X2 <- matrix(rep(x2, n1), n1, n2, byrow = T)
  K1 <- matrix(as.numeric(X1 == pos_class), n1, n2)
  K2 <- matrix(as.numeric(X2 == pos_class), n1, n2)
  return(alpha^2 * K1 * K2)
}


#' @describeIn kernel Computes variance mask multiplier matrix.
#'
#' @param vm_params vector of two mask function parameters
#' @param nan_replace value to replace \code{NaN}s in \code{x1} and \code{x2}
kernel_var_mask <- function(x1, x2, vm_params,
                            a, nan_replace = 0) {
  x1[is.nan(x1)] <- nan_replace
  x2[is.nan(x2)] <- nan_replace
  check_interval(vm_params[1], 0, 1)
  check_positive(vm_params[1])
  check_positive(vm_params[2])
  stp <- a * vm_params[2] # steepness of mask function
  r <- 1 / stp * log(vm_params[1] / (1 - vm_params[1]))
  s1 <- var_mask(x1 - r, stp)
  s2 <- var_mask(x2 - r, stp)
  M <- tcrossprod(s1, s2)
  return(M)
}


#' @describeIn kernel Computes the heterogeneity multiplier matrix.
#'
#' @param beta a parameter vector (row vector) of length \code{N_cases}
#' @param row_to_caseID_1 mapping from row index to case ID (integer vector of
#' length \eqn{n})
#' @param row_to_caseID_2 mapping from row index to case ID (integer vector of
#' length \eqn{m})
kernel_beta <- function(beta, row_to_caseID_1, row_to_caseID_2) {
  n1 <- length(row_to_caseID_1)
  n2 <- length(row_to_caseID_2)
  BETA <- matrix(0, n1, n2)
  for (i in 1:n1) {
    i_case <- row_to_caseID_1[i]
    if (i_case > 0) {
      b1 <- beta[i_case]
    } else {
      b1 <- 0
    }
    for (j in 1:n2) {
      j_case <- row_to_caseID_2[j]
      if (j_case > 0) {
        b2 <- beta[j_case]
      } else {
        b2 <- 0
      }
      BETA[i, j] <- sqrt(b1 * b2)
    }
  }
  return(BETA)
}

#' Input warping function
#'
#' @param x a vector of length \eqn{n}
#' @inheritParams kernel
#' @return a vector of warped inputs \eqn{w(x)}, length \eqn{n}
#' @family kernel utility functions
warp_input <- function(x, a) {
  b <- 0.0
  c <- 1.0
  w <- 2 * c * (-0.5 + 1 / (1 + exp(-a * (x - b))))
  return(w)
}

#' Variance masking function
#'
#' @param x a vector of length \eqn{n}
#' @param stp a positive real number (steepness of mask function)
#' @return a vector of length \eqn{n}
#' @family kernel utility functions
var_mask <- function(x, stp) {
  check_positive(stp)
  y <- 1 / (1 + exp(-stp * x))
  return(y)
}
