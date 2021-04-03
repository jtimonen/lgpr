#' Posterior predictions and function posteriors
#'
#' @description
#' \itemize{
#'   \item If \code{fit} is for a model that marginalizes the latent
#'   signal \code{f}, this computes the analytical conditional posterior
#'   distributions of each model component, their sum, and the conditional
#'   predictive distribution. All these are computed for
#'   each (hyper)parameter draw (defined by \code{draws}), or other parameter
#'   set (obtained by a reduction defined by \code{reduce}). Results are stored
#'   in a \linkS4class{GaussianPrediction} object which is then returned.
#'
#'   \item If \code{fit} is for a model that samples the latent
#'   signal components (and their sum \code{f}), this will in addition to
#'   these function samples, compute versions of the sum \code{f} samples
#'   which are transformed through the inverse link function.
#'   These are stored in a \linkS4class{Prediction}
#'   object which is then returned.
#' }
#'
#' @export
#' @param fit An object of class \linkS4class{lgpfit}.
#' @param x A data frame of points where function posterior distributions
#' and predictions should be computed or sampled.
#' The function \code{\link{new_x}} provides an easy way to create it.
#' If this is NULL, the data points are used.
#' @param c_hat_pred This is only used if the latent signal \code{f} was
#' sampled. This input contains the values added to the sum \code{f} before
#' passing through inverse link function. Must be a vector with length equal to
#' the number of prediction points. If original \code{c_hat} was constant,
#' then \code{c_hat_pred} can be ignored, in which case this will by default
#' use the same constant.
#' @param reduce Reduction for parameters draws. Can be a function that
#' is applied to reduce all parameter draws into one parameter set, or
#' NULL (no reduction). Has no effect if \code{draws} is specified.
#' @param draws Indices of parameter draws to use, or \code{NULL} to use all
#' draws.
#' @param STREAM External pointer. By default obtained with
#' \code{\link[rstan]{get_stream}}.
#' @param verbose Should more information and a possible progress bar be
#' printed?
#' @param ... optional arguments passed to \code{\link{posterior_f}}
#' @return An object of class \linkS4class{GaussianPrediction} or
#' \linkS4class{Prediction}.
pred <- function(fit,
                 x = NULL,
                 reduce = function(x) base::mean(x),
                 draws = NULL,
                 verbose = TRUE,
                 STREAM = get_stream(),
                 c_hat_pred = NULL,
                 ...) {
  f_sampled <- is_f_sampled(fit)
  if (!is.null(draws)) reduce <- NULL
  fp <- posterior_f(
    fit = fit, x = x, reduce = reduce, draws = draws,
    verbose = verbose, STREAM = STREAM, ...
  )
  if (f_sampled) {
    out <- pred_kr(fit, fp, verbose, STREAM, c_hat_pred)
  } else {
    out <- pred_gaussian(fit, fp, verbose)
  }
  return(out)
}

# pred when sample_f = FALSE
pred_gaussian <- function(fit, fp, verbose) {
  if (verbose) {
    cat("Computing preditive distribution on original data scale...\n")
  }
  f_mean <- dollar(fp, "f_mean")
  f_std <- dollar(fp, "f_std")
  sigma2 <- dollar(fp, "sigma2")
  y_scl <- dollar(fit@model@var_scalings, "y")
  y_pred <- pred_gaussian.f_to_y(f_mean, f_std, sigma2, y_scl)

  if (verbose) cat("Done.\n")
  new("GaussianPrediction",
    f_comp_mean = dollar(fp, "f_comp_mean"),
    f_comp_std = dollar(fp, "f_comp_std"),
    f_mean = f_mean,
    f_std = f_std,
    y_mean = dollar(y_pred, "mean"),
    y_std = dollar(y_pred, "std"),
    x = dollar(fp, "x")
  )
}

# Tranform distribution of f to distribution of y
pred_gaussian.f_to_y <- function(f_mean, f_std, sigma2, y_scl) {

  # Compute y_mean and y_std on normalized scale
  y_mean <- f_mean
  y_var <- add_to_columns(f_std^2, sigma2)
  y_std <- sqrt(y_var)

  # Scale y_mean and y_std to original scale and return
  list(
    mean = apply_scaling(y_scl, y_mean, inverse = TRUE),
    std = y_scl@scale * y_std
  )
}

# pred when sample_f = TRUE
pred_kr <- function(fit, fp, verbose, STREAM, c_hat_pred) {
  kernels <- pred_kernels(fit, x, reduce, draws, verbose, STREAM)
  si <- get_stan_input(fit)
  delta <- dollar(si, "delta")
  if (verbose) cat("Extracting sampled function components...\n")
  pred <- get_pred(fit, draws = draws, reduce = reduce)
  if (verbose) cat("Computing kernel regression...\n")
  kr <- pred.kr_compute(kernels, pred, delta, verbose, STREAM)
  f <- dollar(kr, "f")
  if (verbose) cat("Done.\n")

  h <- pred.kr_h(fit, f, c_hat_pred, verbose)

  # Return
  new("Prediction",
    f_comp = arr3_to_list(dollar(kr, "f_comp")),
    f = f,
    h = h
  )
}

# Map the sum f from pred.kr_compute to h
#
# @inheritParams pred
# @param f an array of shape (num_draws, num_pred_points)
# @return an array with same shape as \code{f}
pred.kr_h <- function(fit, f, c_hat_pred, verbose) {

  # helper function
  is_constant <- function(x) {
    s <- sum(x == x[1])
    s == length(x)
  }

  num_draws <- dim(f)[1]
  num_pred <- dim(f)[2]
  if (is.null(c_hat_pred)) {
    c_hat_data <- get_chat(fit)
    if (is_constant(c_hat_data)) {
      msg <- paste0(
        "c_hat_pred not given,",
        " using constant c_hat_pred = ", c_hat_data[1], "\n"
      )
      if (verbose) cat(msg)
      c_hat_pred <- rep(c_hat_data[1], num_pred)
    } else {
      msg <- paste0(
        "c_hat (at data points) is not constant! ",
        "you must give c_hat_pred (at prediction points) as input!"
      )
      stop(msg)
    }
  }

  f <- f + repvec(c_hat_pred, num_draws)
  h <- link_inv(f, get_obs_model(fit))
  return(h)
}

# Compute out-of-sample predictions using kernel regression on
# sampled function values
#
# @param pred An object of class \linkS4class{Prediction}.
# @inheritParams pred.gaussian_compute
# @return An object of class \linkS4class{Prediction}.
pred.kr_compute <- function(kernels, pred, delta, verbose,
                            STREAM = get_stream()) {
  K <- dollar(kernels, "data_vs_data")
  Ks <- dollar(kernels, "pred_vs_data")
  Kss <- dollar(kernels, "pred_vs_pred")
  f_comp <- pred@f_comp # list, each elem has shape num_draws x num_obs
  num_draws <- dim(Kss)[1]
  D <- dim(Kss)[2]
  num_obs <- dim(K)[3]
  num_pred <- dim(Kss)[3]
  out <- array(0, dim = c(D + 1, num_draws, num_pred))
  DELTA <- delta * diag(num_obs)
  hdr <- progbar_header(num_draws)
  idx_print <- dollar(hdr, "idx_print")
  if (verbose) cat(dollar(hdr, "header"), "\n")
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
    if (verbose) progbar_print(i, idx_print)
  }
  if (verbose) cat("\n")

  # Return
  list(
    f_comp = out[1:D, , , drop = FALSE],
    f = arr3_select(out, D + 1)
  )
}
