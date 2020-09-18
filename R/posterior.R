#' Posterior predictive distribution at test points
#'
#' @export
#' @param fit An object of class \linkS4class{lgpfit}.
#' @param X_pred A data frame of points where predictions are computed.
#' The function \code{\link{new_x}} can help in creating it.
#' @param STREAM External pointer. By default obtained with
#' \code{\link[rstan]{get_stream}}.
#' @return A list.
#' @family GP posterior computation functions
posterior_predict <- function(fit, X_pred, STREAM = get_stream()) {
  check_type(X_pred, "data.frame")
  obs_model <- get_obs_model(fit)
  if (obs_model != "gaussian") stop("observation model must be gaussian!")
  sigma <- get_draws(fit, pars = "sigma", stack_chains = TRUE)
  kernels <- posterior_predict_kernels(fit, X_pred, STREAM)
  si <- get_stan_input(fit)
  delta <- dollar(si, "delta")
  y <- get_y(fit)
  y <- as.vector(get_y(fit))
  compute_gp_posteriors(kernels, y, sigma, delta, STREAM)
}

#' Compute all kernel matrices required for computing the function
#' posterior distributions at test points
#'
#' @inheritParams posterior_predict
#' @return A list.
#' @family GP posterior computation functions
posterior_predict_kernels <- function(fit, X_pred, STREAM = get_stream()) {
  model <- object_to_model(fit)
  x_data <- model@stan_input
  x_pred <- create_kernel_input(model, X_pred)
  theta <- get_draws_kernel_pars(fit)

  # Computation
  K <- compute_kernel_matrices(model, x_data, x_data, theta, STREAM)
  Ks <- compute_kernel_matrices(model, x_pred, x_data, theta, STREAM)
  Kss <- compute_kernel_matrices(model, x_pred, x_pred, theta, STREAM)

  # Return
  list(
    data_vs_data = K,
    pred_vs_data = Ks,
    pred_vs_pred = Kss
  )
}

#' Kernel input information to Stan format
#'
#' @param model An object of class \linkS4class{lgpmodel}.
#' @inheritParams posterior_predict
#' @return A list.
#' @family GP posterior computation functions
create_kernel_input <- function(model, X_pred) {
  proc <- parse_covs_and_comps(X_pred, model@model_formula)
  dollar(proc, "to_stan")
}

#' Compute constant kernel matrices
#'
#' @param model An object of class \linkS4class{lgpmodel}.
#' @param x1 A list defining the first input argument. Can be
#' created by \code{\link{create_kernel_input}}.
#' @param x2 A list defining the second input argument. Can be
#' created by \code{\link{create_kernel_input}}.
#' @inheritParams posterior_predict
#' @return A list of kernel matrices with length equal to number of components.
#' @family GP posterior computation functions
compute_kernel_matrices_const <- function(model, x1, x2,
                                          STREAM = get_stream()) {

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
  kernel_const_all(
    n1, n2, x1_cat, x2_cat, x1_cont_mask, x2_cont_mask,
    num_levels, components, STREAM
  )
}

#' Extract kernel parameter draws
#'
#' @inheritParams posterior_predict
#' @family GP posterior computation functions
#' @return a list with names  "alpha", "ell", "wrp", "beta" and "teff", where
#' each element is an array with number of rows equal to number of posterior
#' samples
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

#' Compute kernel matrices for each posterior draw
#'
#' @inheritParams compute_kernel_matrices_const
#' @param theta A list of parameter draws with length \code{num_draws}.
#' @inheritParams posterior_predict
#' @return Returns a 4-dimensional array of shape \code{num_draws} x
#' \code{num_comps} x \code{n1} x \code{n2}.
#' @family GP posterior computation functions
compute_kernel_matrices <- function(model, x1, x2, theta,
                                    STREAM = get_stream()) {

  # Compute constant matrices
  K_const <- compute_kernel_matrices_const(model, x1, x2, STREAM)

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
    K_i <- kernel_all(
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

#' Compute GP posteriors
#'
#' @export
#' @param kernels a list of arrays returned by
#' \code{\link{compute_kernel_matrices}}
#' @param y response variable vector of length \code{num_obs}
#' @param sigma a vector with length \code{num_draws}
#' @param delta jitter to ensure positive definite matrices
#' @inheritParams posterior_predict
#' @return A list.
#' @family GP posterior computation functions
compute_gp_posteriors <- function(kernels, y, sigma, delta,
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
    post <- gp_posterior(k, ks, kss, y, delta, sigma[i], STREAM)
    out[[i]] <- post
  }
  return(out)
}
