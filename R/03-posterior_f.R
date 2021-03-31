#' Evaluate or extract function posterior distributions
#'
#' @description Evaluates (or extracts) the posterior distribution of the
#' total signal \code{f} and its additive components (or draws from these
#' distributions). All these are computed for
#' each parameter draw (defined by \code{draws}), or other parameter set
#' (obtained by a reduction defined by \code{reduce}).
#'
#' @param fit An object of class \linkS4class{lgpfit}.
#' @param x A data frame of points where the predictions are computed.
#' The function \code{\link{new_x}} can help in creating it.
#' @param c_hat_pred This is only used if the latent signal \code{f} was
#' sampled. This input contains the values added to the sum \code{f} before
#' passing through inverse link function. Must be a vector with length equal to
#' the number of prediction points. If original \code{c_hat} was constant,
#' then \code{c_hat_pred} can be ignored, in which case this will by default
#' use the same constant.
#' @param reduce Reduction for parameters draws. Can be a function that
#' is applied to reduce all parameter draws into one parameter set, or
#' \code{NULL} (no reduction). Has no effect if \code{draws} is specified.
#' @param draws Indices of parameter draws to use, or \code{NULL} to use all
#' draws.
#' @return An object of class \linkS4class{FunctionPosterior} or
#' \linkS4class{FunctionDraws}.
#' @param verbose Should more some informational messages be printed?
#' @param debug_km Should this only return the required kernel matrices.
#' Can be used for debugging or testing \code{\link{kernels_fpost}}.
#' @param debug_dims Should this print dimensions of some variables.
#' Can be used for debugging \code{\link{kernels_fpost}}.
#' @param force This is by default \code{FALSE} to prevent unnecessarily
#' large computations that might crash R or take forever.
#' @param STREAM an external pointer
posterior_f <- function(fit,
                        x = NULL,
                        c_hat_pred = NULL,
                        reduce = function(x) base::mean(x),
                        draws = NULL,
                        verbose = TRUE,
                        debug_km = FALSE,
                        debug_dims = FALSE,
                        force = FALSE,
                        STREAM = get_stream()) {

  # Settings
  if (is.null(x)) x <- get_data(fit)
  if (!is.null(draws)) reduce <- NULL
  f_sampled <- is_f_sampled(fit)

  # Stop if potentially going to crash the computer
  prevent_too_large_mats(fit, x, reduce, draws, verbose, force)

  # Compute all required kernel matrices
  km <- kernels_fpost(fit, x, reduce, draws, verbose, debug_dims, STREAM)
  if (debug_km) {
    return(km)
  }

  # Compute the function posteriors
  if (f_sampled) {
    chp <- c_hat_pred
    out <- fp_latent(km, fit, x, chp, reduce, draws, verbose, STREAM)
  } else {
    out <- fp_marginal(km, fit, x, reduce, draws, verbose, STREAM)
  }
  return(out)
}

# Analytic function posteriors
fp_marginal <- function(km, fit, x, reduce, draws, verbose, STREAM) {

  # Fetch sigma2, delta and normalized y
  y <- get_y(fit, original = FALSE)
  d_sigma <- get_draws(fit, pars = "sigma[1]", reduce = reduce, draws = draws)
  sigma2 <- as.vector(d_sigma^2)
  delta <- dollar(get_stan_input(fit), "delta")

  # Get kernel matrices
  K <- dollar(km, "K")
  Ks <- dollar(km, "Ks")
  Kss <- dollar(km, "Kss")
  S <- dim(Ks)[1] # number of parameter sets

  # Setup
  fp <- list()
  progbar <- verbose && S > 1
  pb <- progbar_setup(L = S)
  hdr <- dollar(pb, "header")
  idx_print <- dollar(pb, "idx_print")
  if (verbose) cat("Computing analytic function posteriors...\n")
  if (progbar) cat(hdr)

  # Loop through parameter sets
  for (idx in seq_len(S)) {

    # Perform computations for one parameter set
    K_i <- K[idx, , , ]
    Ks_i <- Ks[idx, , , ]
    Kss_i <- Kss[idx, , , ]
    fp_i <- fp_marginal.compute(K_i, Ks_i, Kss_i, sigma2[idx], delta, y)

    # Format as data.frame and update progress
    fp[[idx]] <- fp_marginal.format(fp_i, dollar(km, "comp_names"))
    if (progbar) progbar_print(idx, idx_print)
  }
  if (progbar) cat("\n")
  if (verbose) cat("\n")

  # Return
  new("FunctionPosterior",
    f = fp,
    x = x,
    model = fit@model,
    num_paramsets = length(fp),
    sigma2 = sigma2
  )
}

# Compute componentwise and total function posteriors
fp_marginal.compute <- function(K, Ks, Kss, sigma2, delta, y) {

  # Helper function for linear algebra
  gp_posterior_helper <- function(Ly, Ks, Kss_diag, v) {
    P <- length(Kss_diag)
    A <- t(forwardsolve(Ly, t(Ks)))
    f_post <- matrix(0, P, 2)
    f_post[, 1] <- A %*% v # mean
    f_post[, 2] <- sqrt(Kss_diag - rowSums(A * A)) # sd
    return(f_post)
  }

  # Function that sums along the first dimension of a 3D array
  matsum1 <- function(A) {
    J <- dim(A)[1]
    A_sum <- A[1, , ]
    for (j in 2:J) A_sum <- A_sum + A[j, , ]
    return(A_sum)
  }
  Kss_diag <- t(apply(Kss, 1, diag)) # has the diagonals as rows

  # Setup output arrays
  J <- dim(Ks)[1] # number of components
  P <- dim(Ks)[2] # number of output points
  N <- dim(Ks)[3] # number of data points
  F_MU <- matrix(0, P, J + 1)
  F_SD <- matrix(0, P, J + 1)

  # Compute Ky and Cholesky decompose it
  Ky <- matsum1(K) + (delta + sigma2) * diag(N)
  Ly <- t(chol(Ky))
  v <- forwardsolve(Ly, y)

  # Component-wise means and sds
  for (j in seq_len(J)) {
    fp_j <- gp_posterior_helper(Ly, Ks[j, , ], Kss_diag[j, ], v)
    F_MU[, j] <- fp_j[, 1]
    F_SD[, j] <- fp_j[, 2]
  }

  # Total mean and sd
  Ks_sum <- matsum1(Ks)
  fp_sum <- gp_posterior_helper(Ly, Ks_sum, colSums(Kss_diag), v)
  F_MU[, J + 1] <- fp_sum[, 1]
  F_SD[, J + 1] <- fp_sum[, 2]

  # Return
  list(mean = F_MU, sd = F_SD)
}

# Format fp_marginal computation results as a data frame
fp_marginal.format <- function(fp, comp_names) {
  m <- dollar(fp, "mean") # shape (P, J+1)
  s <- dollar(fp, "sd") # shape (P, J+1)
  P <- dim(m)[1]
  J <- dim(m)[2] - 1
  comp_names <- c(comp_names, "f_sum")
  component <- as.factor(rep(comp_names, each = P))
  eval_point_idx <- as.factor(rep(1:P, times = J + 1))
  m <- as.vector(m)
  s <- as.vector(s)
  check_lengths(m, s)
  check_lengths(m, component)
  check_lengths(m, eval_point_idx)
  df <- data.frame(eval_point_idx, component, m, s)
  colnames(df) <- c("eval_point_idx", "component", "mean", "sd")
  return(df)
}

# Compute total posterior
fp_marginal.total <- function(km, sigma2, y_norm, delta) {
  stop("NOT IMPLEMENTED!") # TODO
  S <- dim(m)[2] # number of param sets
  P <- dim(m)[3] # number of prediction points
  paramset <- as.factor(rep(1:S, P))
  eval_point <- as.factor(rep(1:P, each = S))
  sigma <- rep(sigma, P)
  m <- as.vector(m)
  s <- as.vector(s)
  check_lengths(m, s)
  check_lengths(m, paramset)
  check_lengths(m, eval_point)
  check_lengths(m, sigma)
  df <- data.frame(paramset, eval_point, m, s, sigma)
  colnames(df) <- c("paramset", "eval_point", "mean", "std", "sigma")
  return(df)
}

# Function posterior draws
fp_latent <- function(km, fit, x, c_hat_pred, reduce, draws, verbose, STREAM) {
  fpred <- NULL
  stop("NOT IMPLEMENTED!") # TODO

  # Create the prediction
  # h <- pred.latent_h(fit, f_pred, c_hat_pred, verbose)
  #
  ## Return
  # new("Prediction",
  #    f_comp = arr3_to_list(dollar(kr, "f_comp")),
  #    f = f,
  #    h = h
  # )
  return(f_pred)
}

# Transform distribution of f to distribution of y
pred_marginal.f_to_y <- function(f_mean, f_std, sigma, y_norm_inv) {
  y_mean <- f_mean
  y_var <- add_to_columns(f_std^2, sigma^2)
  y_std <- sqrt(y_var)
  y_upper <- y_mean + y_std

  # Scale y_pred to original scale
  y_mean <- call_fun(y_norm_inv, y_mean)
  y_upper <- call_fun(y_norm_inv, y_upper)
  y_std <- y_upper - y_mean

  # Return
  list(mean = y_mean, std = y_std)
}

# Map the sum f from pred.latent_compute to h
pred.latent_h <- function(fit, f, c_hat_pred, verbose) {

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

# Safeguard
prevent_too_large_mats <- function(fit, x, reduce, draws, verbose, force) {
  S <- determine_num_paramsets(fit, draws, reduce)
  J <- length(component_names(fit@model))
  N <- get_num_obs(fit)
  P <- nrow(x)
  msg <- paste0(
    "Computations will require creating and handling ",
    P, " x ", N, " matrices ", S, " x ", J + 1,
    " times. \n"
  )
  P_LIMIT <- 6000
  if (verbose) cat(msg)
  if (P > P_LIMIT) {
    if (!force) {
      msg <- paste0(
        msg, "Computations might take very long. Give a smaller number of ",
        "output points (now ", P, "), or set force = TRUE if you are sure ",
        "you want to do this."
      )
      stop(msg)
    }
  }
  TRUE
}
