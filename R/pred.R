#' Computing out-of-sample predictions
#'
#' @description
#' \itemize{
#'   \item \code{pred} calls either \code{pred.gaussian} or
#'   \code{pred.kr} depending on whether \code{fit} is for a model
#'   that marginalized or sampled the latent signal \code{f}
#'   \item \code{pred.gaussian} computes analytic posterior or prior
#'   predictive distribution and returns an object of class
#'   \linkS4class{GaussianPrediction}
#'   \item \code{pred.kr} computes predictions using kernel regression for
#'   the samples of \code{f} and returns an object of class
#'   \linkS4class{Prediction}
#' }
#' These also compute the predicted (distribution of the) total signal
#' \code{f} and its additive components. All these are computed for
#' each parameter draw (defined by \code{draws}), or other parameter set
#' (obtained by a reduction defined by \code{reduce}).
#'
#' @param fit An object of class \linkS4class{lgpfit}.
#' @param x A data frame of points where the predictions are computed.
#' The function \code{\link{new_x}} provides an easy way to create it.
#' @param reduce reduction for parameters draws
#' @param draws indices of parameter draws
#' @param STREAM External pointer. By default obtained with
#' \code{\link[rstan]{get_stream}}.
#' @return An object of class \linkS4class{GaussianPrediction} or
#' \linkS4class{Prediction}.
#' @family prediction functions
#' @name pred
NULL

#' @export
#' @rdname pred
pred <- function(fit, x, reduce = mean, draws = NULL,
                 STREAM = get_stream()) {
  check_type(x, "data.frame")
  f_sampled <- is_f_sampled(fit)
  if (f_sampled) {
    out <- pred.kr(fit, x, reduce, draws, STREAM)
  } else {
    out <- pred.gaussian(fit, x, reduce, draws, STREAM)
  }
  return(out)
}

#' @rdname pred
pred.gaussian <- function(fit, x, reduce, draws, STREAM = get_stream()) {
  sigma <- get_draws(fit, draws = draws, reduce = reduce, pars = "sigma")
  kernels <- pred_kernels(fit, x, reduce, draws, STREAM)
  si <- get_stan_input(fit)
  delta <- dollar(si, "delta")
  y_norm <- as.vector(get_y(fit, original = FALSE))
  post <- pred.gaussian_compute(kernels, y_norm, sigma, delta, STREAM)

  f_mean <- dollar(post, "f_mean")
  f_std <- dollar(post, "f_std")
  y_scl <- dollar(fit@model@var_scalings, "y")
  y_pred <- pred.gaussian.f_to_y(f_mean, f_std, sigma, y_scl@fun_inv)

  new("GaussianPrediction",
    f_comp_mean = arr3_to_list(dollar(post, "f_comp_mean")),
    f_comp_std = arr3_to_list(dollar(post, "f_comp_std")),
    f_mean = f_mean,
    f_std = f_std,
    y_mean = dollar(y_pred, "mean"),
    y_std = dollar(y_pred, "std")
  )
}

#' @rdname pred
pred.kr <- function(fit, x, reduce, draws, STREAM = get_stream()) {
  kernels <- pred_kernels(fit, x, reduce, draws, STREAM)
  si <- get_stan_input(fit)
  delta <- dollar(si, "delta")
  pred <- get_pred(fit, draws = draws, reduce = reduce)
  kr <- pred.kr_compute(kernels, pred, delta, STREAM)
  f <- dollar(kr, "f")
  likelihood <- get_obs_model(fit)

  # Return
  new("Prediction",
    f_comp = arr3_to_list(dollar(kr, "f_comp")),
    f = f,
    h = link_inv(f, likelihood)
  )
}

#' Compute all kernel matrices required for computing predictions
#'
#' @inheritParams pred
#' @return Returns a list with fields named \code{data_vs_data},
#' \code{pred_vs_data} and \code{pred_vs_pred}. Each of them is a list with
#' length equal to the number of parameter sets.
#' @family prediction functions
pred_kernels <- function(fit, x, reduce, draws, STREAM = get_stream()) {
  model <- object_to_model(fit)
  x_data <- model@stan_input
  x_pred <- kernel_matrices_create_input(model, x)
  theta <- get_kernel_pars(fit, reduce, draws)

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

#' Compute analytic out-of-sample predictions
#'
#' @param kernels A list returned by \code{\link{pred_kernels}}.
#' @param y response variable vector of length \code{num_obs}
#' @param sigma a vector with length equal to number of parameter sets
#' @param delta jitter to ensure positive definite matrices
#' @inheritParams pred
#' @return A list.
#' @family prediction functions
pred.gaussian_compute <- function(kernels, y, sigma, delta,
                                  STREAM = get_stream()) {
  sigma <- as.vector(sigma)
  K <- dollar(kernels, "data_vs_data")
  Ks <- dollar(kernels, "pred_vs_data")
  Kss <- dollar(kernels, "pred_vs_pred")
  num_draws <- dim(Kss)[1]
  D <- dim(Kss)[2]
  R <- D + 1
  num_pred <- dim(Kss)[3]
  out <- array(0, dim = c(2 * R, num_pred, num_draws))
  for (i in seq_len(num_draws)) {
    k <- arr3_to_list(K[i, , , ])
    ks <- arr3_to_list(Ks[i, , , ])
    kss <- arr3_to_list(Kss[i, , , ])
    post <- cpp_gp_posterior(k, ks, kss, y, delta, sigma[i], STREAM)
    post <- list_to_matrix(post, num_pred)
    out[, , i] <- post
  }
  out <- aperm(out, c(1, 3, 2))
  # Return
  list(
    f_comp_mean = out[1:D, , , drop = FALSE],
    f_comp_std = out[(R + 1):(R + D), , , drop = FALSE],
    f_mean = arr3_select(out, R),
    f_std = arr3_select(out, 2 * R)
  )
}

#' Compute out-of-sample predictions using kernel regression on
#' sampled function values
#'
#' @param pred An object of class \linkS4class{Prediction}.
#' @inheritParams pred.gaussian_compute
#' @return An object of class \linkS4class{Prediction}.
#' @family prediction functions
pred.kr_compute <- function(kernels, pred, delta, STREAM = get_stream()) {
  K <- dollar(kernels, "data_vs_data")
  Ks <- dollar(kernels, "pred_vs_data")
  Kss <- dollar(kernels, "pred_vs_pred")
  f_comp <- pred@f_comp # list, each elem has shape num_draws x num_obs
  num_draws <- dim(Kss)[1]
  D <- dim(Kss)[2]
  num_obs <- dim(K)[3]
  num_pred <- dim(Kss)[3]
  out <- array(0, dim = c(D + 1, num_draws, num_pred))
  DELTA <- DELTA <- delta * diag(num_obs)
  for (i in seq_len(num_draws)) {
    f_sum <- 0
    for (j in seq_len(D)) {
      fj <- f_comp[[j]]
      k <- K[i, j, , ] + DELTA
      ks <- Ks[i, j, , ]
      f <- fj[i, ]
      f_pred <- ks %*% solve(k, f)
      out[j, i, ] <- f_pred
      f_sum <- f_sum + f_pred
    }
    out[D + 1, i, ] <- f_sum
  }

  # Return
  list(
    f_comp = out[1:D, , , drop = FALSE],
    f = arr3_select(out, D + 1)
  )
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
#' }
#' @param model An object of class \linkS4class{lgpmodel}.
#' @param x1 A list defining the first input argument. Can be
#' created by \code{\link{kernel_matrices_create_input}}.
#' @param x2 A list defining the second input argument. Can be
#' created by \code{\link{kernel_matrices_create_input}}.
#' @inheritParams pred
#' @name kernel_matrices
#' @family prediction functions
NULL

#' @rdname kernel_matrices
#' @param theta A list of kernel parameter sets.
kernel_matrices <- function(model, x1, x2, theta, STREAM = get_stream()) {

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

#' Extract kernel parameter draws
#'
#' @return a list with names "alpha", "ell", "wrp", "beta" and "teff",
#' where each element is an array with number of rows equal to number of
#' parameter sets
#' @inheritParams pred
get_kernel_pars <- function(fit, reduce = NULL, draws = NULL) {
  num_draws <- get_num_draws(fit)
  na_array <- array(NaN, dim = c(num_draws, 1))
  par_names <- c("alpha", "ell", "wrp", "beta", "teff")
  out <- list()
  for (pn in par_names) {
    s <- get_draws.catch(fit, draws = draws, reduce = reduce, pars = pn)
    out[[pn]] <- if (is.null(s)) na_array else s
  }
  return(out)
}
