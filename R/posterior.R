#' Compute total and componentwise GP posteriors at given points
#'
#' @description
#' \itemize{
#'   \item Total and componentwise GP posteriors at test points can be
#'   computed using \code{gp_posteriors}. This function requires
#'   the observation model to be Gaussian. The actual work is done
#'   in \code{\link{gp_posteriors_compute}}. Return?
#'   \item All kernel matrices required for computing posterior distributions
#'   at test points are computed by \code{gp_kernels}. It returns a list
#'   with fields named \code{data_vs_data}, \code{pred_vs_data} and
#'   \code{pred_vs_pred}. Each of them is a list with length equal to the
#'   number of posterior draws. The actual work is done
#'   in \code{\link{kernel_matrices}}
#' }
#' @param fit An object of class \linkS4class{lgpfit}.
#' @param x A data frame of points where the posteriors are computed.
#' The function \code{\link{new_x}} provides an easy way to create it.
#' @param STREAM External pointer. By default obtained with
#' \code{\link[rstan]{get_stream}}.
#' @return See description.
#' @name gp_posteriors
#' @family GP posterior computation functions
NULL

#' @export
#' @rdname gp_posteriors
gp_posteriors <- function(fit, x, STREAM = get_stream()) {
  check_type(x, "data.frame")
  obs_model <- get_obs_model(fit)
  if (obs_model != "gaussian") stop("observation model must be gaussian!")
  sigma <- get_draws(fit, pars = "sigma", stack_chains = TRUE)
  kernels <- gp_kernels(fit, x, STREAM)
  si <- get_stan_input(fit)
  delta <- dollar(si, "delta")
  y <- get_y(fit)
  y <- as.vector(get_y(fit))
  gp_posteriors_compute(kernels, y, sigma, delta, STREAM)
}

#' @export
#' @rdname gp_posteriors
gp_kernels <- function(fit, x, STREAM = get_stream()) {
  model <- object_to_model(fit)
  x_data <- model@stan_input
  x_pred <- kernel_matrices_create_input(model, x)
  theta <- get_draws_kernel_pars(fit)

  # Computation
  K <- kernel_matrices(model, x_data, x_data, theta, STREAM)
  Ks <- kernel_matrices(model, x_pred, x_data, theta, STREAM)
  Kss <- kernel_matrices(model, x_pred, x_pred, theta, STREAM)

  # Return
  list(
    data_vs_data = K,
    pred_vs_data = Ks,
    pred_vs_pred = Kss
  )
}

#' Compute GP posteriors
#'
#' @param kernels a list of arrays returned by \code{\link{kernel_matrices}}
#' @param y response variable vector of length \code{num_obs}
#' @param sigma a vector with length \code{num_draws}
#' @param delta jitter to ensure positive definite matrices
#' @inheritParams gp_posteriors
#' @return A list.
#' @family GP posterior computation functions
gp_posteriors_compute <- function(kernels, y, sigma, delta,
                                  STREAM = get_stream()) {
  sigma <- as.vector(sigma)
  K <- dollar(kernels, "data_vs_data")
  Ks <- dollar(kernels, "pred_vs_data")
  Kss <- dollar(kernels, "pred_vs_pred")
  L <- length(sigma)
  arr3_to_list <- function(x) {
    out <- list()
    d <- dim(x)[1]
    for (j in seq_len(d)) {
      out[[j]] <- x[j, , ]
    }
    return(out)
  }
  out <- list()
  for (i in seq_len(L)) {
    k <- arr3_to_list(K[i, , , ])
    ks <- arr3_to_list(Ks[i, , , ])
    kss <- arr3_to_list(Kss[i, , , ])
    post <- cpp_gp_posterior(k, ks, kss, y, delta, sigma[i], STREAM)
    out[[i]] <- post
  }
  return(out)
}


#' Kernel matrix computation using the C++ functions of the package namespace
#'
#' @description
#' \itemize{
#'   \item Kernel matrices using parameter values defined by each posterior
#'   draw can be computed using \code{kernel_matrices}. It returns a
#'   4-dimensional array of shape \code{num_draws} x \code{num_comps} x
#'   \code{n1} x \code{n2}.
#'   \item The constant kernel matrices that do not depend on parameters
#'   can be computed using \code{kernel_matrices_const}. It returns a list
#'   of kernel matrices with length equal to number of model components.
#'   \item Input information can be translated from a data frame to the
#'   needed format using \code{kernel_matrices_create_input}.
#'   \item Kernel parameter draws can be extracted using
#'   \code{get_draws_kernel_params}. It returns a list with names
#'   "alpha", "ell", "wrp", "beta" and "teff", where each element is an
#'   array with number of rows equal to number of posterior draws.
#' }
#' @param model An object of class \linkS4class{lgpmodel}.
#' @param x1 A list defining the first input argument. Can be
#' created by \code{\link{kernel_matrices_create_input}}.
#' @param x2 A list defining the second input argument. Can be
#' created by \code{\link{kernel_matrices_create_input}}.
#' @inheritParams gp_posteriors
#' @name kernel_matrices
NULL

#' @export
#' @rdname kernel_matrices
#' @param theta A list of kernel parameter draws with length \code{num_draws}.
kernel_matrices <- function(model, x1, x2, theta,
                            STREAM = get_stream()) {

  # Compute constant matrices
  K_const <- kernel_matrices_const(model, x1, x2, STREAM)

  # Get variables related to model
  si <- model@stan_input
  components <- dollar(si, "components")
  vm_params <- dollar(si, "vm_params")
  teff_zero <- dollar(si, "teff_zero")

  # Get variables related to first input
  x1_cont <- dollar(x1, "x_cont")
  x1_unnorm <- dollar(x1, "x_cont_unnorm")
  idx1_expand <- dollar(x1, "idx_expand")
  n1 <- length(idx1_expand)

  # Get variables related to second input
  x2_cont <- dollar(x2, "x_cont")
  x2_unnorm <- dollar(x2, "x_cont_unnorm")
  idx2_expand <- dollar(x2, "idx_expand")
  n2 <- length(idx2_expand)

  # Get parameter draws
  alpha <- dollar(theta, "alpha")
  ell <- dollar(theta, "ell")
  wrp <- dollar(theta, "wrp")
  beta <- dollar(theta, "beta")
  teff <- dollar(theta, "teff")

  # Loop through each posterior draw
  J <- nrow(components)
  L <- nrow(alpha)
  K_out <- array(0.0, dim = c(L, J, n1, n2))
  for (i in seq_len(L)) {
    K_i <- cpp_kernel_all(
      n1, n2, K_const, components,
      x1_cont, x2_cont, x1_unnorm, x2_unnorm,
      alpha[i, ], ell[i, ], wrp[i, ], beta[i, ], teff[i, ],
      vm_params, idx1_expand, idx2_expand, teff_zero,
      STREAM
    )
    for (j in seq_len(J)) {
      K_out[i, j, , ] <- K_i[[j]]
    }
  }
  return(K_out)
}

#' @rdname kernel_matrices
kernel_matrices_const <- function(model, x1, x2, STREAM = get_stream()) {

  # Get variables related to model
  si <- model@stan_input
  num_levels <- dollar(si, "x_cat_num_levels")
  components <- dollar(si, "components")

  # Get variables related to first input
  x1_cat <- dollar(x1, "x_cat")
  x1_cont_mask <- dollar(x1, "x_cont_mask")
  n1 <- length(dollar(x1, "idx_expand"))

  # Get variables related to second input
  x2_cat <- dollar(x2, "x_cat")
  x2_cont_mask <- dollar(x2, "x_cont_mask")
  n2 <- length(dollar(x2, "idx_expand"))

  # Compute matrices and return
  cpp_kernel_const_all(
    n1, n2, x1_cat, x2_cat, x1_cont_mask, x2_cont_mask,
    num_levels, components, STREAM
  )
}

#' @rdname kernel_matrices
kernel_matrices_create_input <- function(model, x) {
  proc <- parse_covs_and_comps(x, model@model_formula)
  dollar(proc, "to_stan")
}

#' @rdname kernel_matrices
get_draws_kernel_pars <- function(fit) {
  num_draws <- get_num_draws(fit)
  na_array <- array(NaN, dim = c(num_draws, 1))
  par_names <- c("alpha", "ell", "wrp", "beta", "teff")
  out <- list()
  for (pn in par_names) {
    draws <- get_draws_catch(fit, stack_chains = TRUE, pars = pn)
    out[[pn]] <- if (is.null(draws)) na_array else draws
  }
  return(out)
}
