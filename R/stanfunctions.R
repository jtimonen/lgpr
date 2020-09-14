#' Compute all constant kernel matrices of a model
#'
#' @description This is a wrapper for \code{STAN_kernel_const_all}
#' @param x1_cat categorical covariates
#' @param x2_cat categorical covariates
#' @param x1_cont_mask continuous covariate masks
#' @param x2_cont_mask continuous covariate masks
#' @param x_cat_num_levels number of levels for each categorical covariate
#' @param components an integer array
#' @param STREAM an external pointer
#' @return an array of kernel matrices
kernel_const_all <- function(x1_cat, x2_cat,
                             x1_cont_mask, x2_cont_mask,
                             x_cat_num_levels, components,
                             STREAM = get_stream()) {
  n1 <- length(x1_cat[1, ])
  n2 <- length(x2_cat[1, ])
  x1_cat <- matrix_to_list(x1_cat)
  x2_cat <- matrix_to_list(x2_cat)
  x1_cont_mask <- matrix_to_list(x1_cont_mask)
  x2_cont_mask <- matrix_to_list(x2_cont_mask)
  components <- matrix_to_list(components)

  STAN_kernel_const_all(
    n1, n2, x1_cat, x2_cat, x1_cont_mask, x2_cont_mask,
    x_cat_num_levels, components, STREAM
  )
}


#' Compute all kernel matrices of a model
#'
#' @description This is a wrapper for \code{STAN_kernel_all}
#' @param K_const list of constant matrices
#' @param components an integer array
#' @param x1 continuous covariates
#' @param x2 continuous covariates
#' @param x1_unnorm unnormalized versions of the continuous covariates
#' @param x2_unnorm unnormalized versions of the continuous covariates
#' @param alpha magnitude parameters
#' @param ell lengthscale parameters
#' @param wrp input warping parameters
#' @param beta heterogeneity parameters
#' @param teff effect time parameters
#' @param vm_params variance mask parameters
#' @param idx1_expand expansion of beta and/or teff
#' @param idx2_expand expansion of beta and/or teff
#' @param teff_obs observed effect times
#' @param STREAM an external pointer
#' @return an array of kernel matrices
kernel_all <- function(K_const, components, x1, x2, x1_unnorm, x2_unnorm,
                       alpha, ell, wrp, beta, teff, vm_params,
                       idx1_expand, idx2_expand, teff_obs,
                       STREAM = get_stream()) {
  n1 <- length(x1[1, ])
  n2 <- length(x2[1, ])

  # List format
  x1 <- matrix_to_list(x1)
  x2 <- matrix_to_list(x2)
  x1_unnorm <- matrix_to_list(x1_unnorm)
  x2_unnorm <- matrix_to_list(x2_unnorm)
  components <- matrix_to_list(components)
  vm_params <- matrix_to_list(vm_params)
  teff_obs <- matrix_to_list(teff_obs)

  # Return
  STAN_kernel_all(
    n1, n2, K_const, components,
    x1, x2, x1_unnorm, x2_unnorm,
    alpha, ell, wrp, beta, teff, vm_params,
    idx1_expand, idx2_expand, teff_obs, STREAM
  )
}

#' Compute all componentwise GP posteriors of a model
#'
#' @description This is a wrapper for \code{STAN_gp_posterior}
#' @param KX list of kernel matrices fo size \code{num_obs} x \code{num_obs}
#' @param KX_s list of kernel matrices fo size \code{num_out} x \code{num_obs}
#' @param KX_ss list of kernel matrices fo size \code{num_out} x \code{num_out}
#' @param y response vector
#' @param delta jitter
#' @param sigma Gaussian noise magnitude parameter
#' @param STREAM an external pointer
#' @return an array of kernel matrices
gp_posterior <- function(KX, KX_s, KX_ss, y, delta, sigma,
                         STREAM = get_stream()) {
  STAN_gp_posterior(KX, KX_s, KX_ss, y, delta, sigma, STREAM)
}

#' Input warping function
#'
#' @description This is a wrapper for \code{STAN_warp_input}
#' @param x a numeric vector
#' @param a steepness parameter
#' @param STREAM an external pointer
#' @return a vector
warp_input <- function(x, a, STREAM = get_stream()) {
  STAN_warp_input(x, a, STREAM)
}

#' Variance masking function
#'
#' @description This is a wrapper for \code{STAN_var_mask}
#' @param x a numeric vector
#' @param a steepness parameter
#' @param STREAM an external pointer
#' @return a vector
var_mask <- function(x, a, STREAM = get_stream()) {
  STAN_var_mask(x, a, STREAM)
}
