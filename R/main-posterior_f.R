#' Evaluate or extract function posterior distributions
#'
#' @description Evaluates (or extracts) the posterior distribution of the
#' total signal \code{f} and its additive components (or draws from these
#' distributions). All these are computed for
#' each parameter draw (defined by \code{draws}), or other parameter set
#' (obtained by a reduction defined by \code{reduce}).
#'
#' @inheritParams pred
#' @param debug_km Should this only return the required kernel matrices.
#' Can be used for debugging or testing \code{\link{kernels_fpost}}.
#' @param debug_dims Should this print dimensions of some variables.
#' Can be used for debugging \code{\link{kernels_fpost}}.
#' @param force This is by default \code{FALSE} to prevent unnecessarily
#' large computations that might crash R or take forever.
#' @return a named list
posterior_f <- function(fit,
                        x = NULL,
                        reduce = function(x) base::mean(x),
                        draws = NULL,
                        verbose = TRUE,
                        STREAM = get_stream(),
                        debug_km = FALSE,
                        debug_dims = FALSE,
                        force = FALSE) {

  # Settings
  if (!is.null(draws)) reduce <- NULL
  f_sampled <- is_f_sampled(fit)

  # Stop if potentially going to crash the computer
  prevent_too_large_mats(fit, x, reduce, draws, verbose, force)

  # Compute all required kernel matrices
  km <- kernels_fpost(fit, x, reduce, draws, verbose, debug_dims, STREAM)
  if (debug_km) {
    return(km)
  }
  if (is.null(x)) x <- get_data(fit@model)

  # Compute the function posteriors
  if (f_sampled) {
    out <- fp_extrapolate(km, fit, x, reduce, draws, verbose, STREAM)
  } else {
    out <- fp_gaussian(km, fit, x, reduce, draws, verbose, STREAM)
  }
  return(out)
}

# Analytic function posteriors
fp_gaussian <- function(km, fit, x, reduce, draws, verbose, STREAM) {

  # Fetch sigma2, delta and normalized y
  y <- get_y(fit, original = FALSE)
  d_sigma <- get_draws(fit, pars = "sigma[1]", reduce = reduce, draws = draws)
  sigma2 <- as.vector(d_sigma^2)
  delta <- dollar(get_stan_input(fit), "delta")

  # Get kernel matrices
  K <- dollar(km, "K")
  Ks <- dollar(km, "Ks")
  Kss <- dollar(km, "Kss")

  # Create output arrays
  S <- length(Ks) # number of parameter sets
  P <- nrow(x) # number of output points
  J <- get_num_comps(fit) # number of components
  f_comp_mean <- array(0, c(S, P, J))
  f_comp_std <- array(0, c(S, P, J))
  f_mean <- array(0, c(S, P))
  f_std <- array(0, c(S, P))

  # Setup
  progbar <- verbose && S > 1
  pb <- progbar_setup(L = S)
  hdr <- dollar(pb, "header")
  idx_print <- dollar(pb, "idx_print")
  log_progress("Computing analytic function posteriors...", verbose)
  log_progress(hdr, progbar)

  # Loop through parameter sets
  for (idx in seq_len(S)) {

    # Perform computations for one parameter set
    K_i <- K[[idx]]
    Ks_i <- Ks[[idx]]
    Kss_i <- Kss[[idx]]
    fp_i <- fp_gaussian.compute(K_i, Ks_i, Kss_i, sigma2[idx], delta, y)
    mean_i <- dollar(fp_i, "mean") # matrix with shape (P, J + 1)
    std_i <- dollar(fp_i, "sd") # matrix with shape (P, J + 1)

    # Store result and update progress
    f_comp_mean[idx, , ] <- mean_i[, 1:J]
    f_comp_std[idx, , ] <- std_i[, 1:J]
    f_mean[idx, ] <- mean_i[, J + 1]
    f_std[idx, ] <- std_i[, J + 1]

    if (progbar) progbar_print(idx, idx_print)
  }
  log_progress(" ", progbar)

  # Return
  comp_names <- dollar(km, "comp_names")
  f_comp_mean <- aperm(f_comp_mean, c(3, 1, 2)) # dim (S, P, J) -> (J, S, P)
  f_comp_std <- aperm(f_comp_std, c(3, 1, 2)) # dim (S, P, J) -> (J, S, P)
  list(
    f_comp_mean = arr3_to_list(f_comp_mean, comp_names), # list with len J
    f_comp_std = arr3_to_list(f_comp_std, comp_names), # list with len J
    f_mean = f_mean, # dim (S, P)
    f_std = f_std, # dim (S, P)
    sigma2 = sigma2,
    x = x
  )
}

# Compute componentwise and total function posteriors
fp_gaussian.compute <- function(K, Ks, Kss, sigma2, delta, y) {

  # Helper function for linear algebra
  gp_posterior_helper <- function(Ly, Ks, Kss_diag, v) {
    P <- length(Kss_diag)
    A <- t(forwardsolve(Ly, t(Ks)))
    f_post <- matrix(0, P, 2)
    f_post[, 1] <- A %*% v # mean
    f_post[, 2] <- sqrt(Kss_diag - rowSums(A * A)) # sd
    return(f_post)
  }

  # Function that sums along a list
  listsum <- function(A) {
    out <- A[[1]]
    J <- length(A)
    if (J > 1) {
      for (j in 2:length(A)) out <- out + A[[j]]
    }
    return(out)
  }
  Kss_diag <- lapply(Kss, diag) # has the diagonals as elements now

  # Setup output arrays
  J <- length(Ks) # number of components
  P <- dim(Ks[[1]])[1] # number of output points
  N <- dim(Ks[[1]])[2] # number of data points
  F_MU <- matrix(0, P, J + 1)
  F_SD <- matrix(0, P, J + 1)

  # Compute Ky and Cholesky decompose it
  Ky <- listsum(K) + (delta + sigma2) * diag(N)
  Ly <- t(chol(Ky))
  v <- forwardsolve(Ly, y)

  # Component-wise means and sds
  for (j in seq_len(J)) {
    fp_j <- gp_posterior_helper(Ly, Ks[[j]], Kss_diag[[j]], v)
    F_MU[, j] <- fp_j[, 1]
    F_SD[, j] <- fp_j[, 2]
  }

  # Total mean and sd
  fp_sum <- gp_posterior_helper(Ly, listsum(Ks), listsum(Kss_diag), v)
  F_MU[, J + 1] <- fp_sum[, 1]
  F_SD[, J + 1] <- fp_sum[, 2]

  # Return
  list(mean = F_MU, sd = F_SD)
}

# Extrapolate function posterior draws using kernel regression
fp_extrapolate <- function(km, fit, x, reduce, draws, verbose, STREAM) {
  stop("NOT IMPLEMENTED")
  # Fetch delta and draws of f
  delta <- dollar(get_stan_input(fit), "delta")
  fp_at_data <- get_pred(fit)

  # Get kernel matrices
  K <- dollar(km, "K")
  Ks <- dollar(km, "Ks")
  Kss <- dollar(km, "Kss")

  # Create output arrays
  S <- length(Ks) # number of parameter sets
  P <- nrow(x) # number of output points
  J <- get_num_comps(fit) # number of components
  f_comp <- array(0, c(S, P, J))
  f_mean <- array(0, c(S, P))
  f_std <- array(0, c(S, P))

  # Setup
  progbar <- verbose && S > 1
  pb <- progbar_setup(L = S)
  hdr <- dollar(pb, "header")
  idx_print <- dollar(pb, "idx_print")
  log_progress(paste0(
    "Extrapolating function posterior draws using ",
    "kernel regression..."
  ), verbose)
  log_progress(hdr, progbar)

  # Loop through parameter sets
  for (idx in seq_len(S)) {

    # Perform computations for one parameter set
    K_i <- K[[idx]]
    Ks_i <- Ks[[idx]]
    Kss_i <- Kss[[idx]]
    fp_i <- fp_extrapolate.compute(K_i, Ks_i, Kss_i, delta, fp_draws)
    mean_i <- dollar(fp_i, "mean") # matrix with shape (P, J + 1)
    std_i <- dollar(fp_i, "sd") # matrix with shape (P, J + 1)

    # Store result and update progress
    f_comp_mean[idx, , ] <- mean_i[, 1:J]
    f_comp_std[idx, , ] <- std_i[, 1:J]
    f_mean[idx, ] <- mean_i[, J + 1]
    f_std[idx, ] <- std_i[, J + 1]

    if (progbar) progbar_print(idx, idx_print)
  }
  log_progress(" ", progbar)

  # Return
  comp_names <- dollar(km, "comp_names")
  f_comp_mean <- aperm(f_comp_mean, c(3, 1, 2)) # dim (S, P, J) -> (J, S, P)
  f_comp_std <- aperm(f_comp_std, c(3, 1, 2)) # dim (S, P, J) -> (J, S, P)
  list(
    f_comp_mean = arr3_to_list(f_comp_mean, comp_names), # list with len J
    f_comp_std = arr3_to_list(f_comp_std, comp_names), # list with len J
    f_mean = f_mean, # dim (S, P)
    f_std = f_std, # dim (S, P)
    sigma2 = sigma2,
    x = x
  )
}

# Extrapolate componentwise and total function posteriors
fp_extrapolate.compute <- function(K, Ks, Kss, delta, f_draws) {

  # Helper function for linear algebra
  gp_posterior_helper <- function(Ly, Ks, Kss_diag, v) {
    P <- length(Kss_diag)
    A <- t(forwardsolve(Ly, t(Ks)))
    f_post <- matrix(0, P, 2)
    f_post[, 1] <- A %*% v # mean
    f_post[, 2] <- sqrt(Kss_diag - rowSums(A * A)) # sd
    return(f_post)
  }

  # Function that sums along a list
  listsum <- function(A) {
    out <- A[[1]]
    J <- length(A)
    if (J > 1) {
      for (j in 2:length(A)) out <- out + A[[j]]
    }
    return(out)
  }
  Kss_diag <- lapply(Kss, diag) # has the diagonals as elements now

  # Setup output arrays
  J <- length(Ks) # number of components
  P <- dim(Ks[[1]])[1] # number of output points
  N <- dim(Ks[[1]])[2] # number of data points
  F_MU <- matrix(0, P, J + 1)
  F_SD <- matrix(0, P, J + 1)

  # Compute Ky and Cholesky decompose it
  Ky <- listsum(K) + (delta + sigma2) * diag(N)
  Ly <- t(chol(Ky))
  v <- forwardsolve(Ly, y)

  # Component-wise means and sds
  for (j in seq_len(J)) {
    fp_j <- gp_posterior_helper(Ly, Ks[[j]], Kss_diag[[j]], v)
    F_MU[, j] <- fp_j[, 1]
    F_SD[, j] <- fp_j[, 2]
  }

  # Total mean and sd
  fp_sum <- gp_posterior_helper(Ly, listsum(Ks), listsum(Kss_diag), v)
  F_MU[, J + 1] <- fp_sum[, 1]
  F_SD[, J + 1] <- fp_sum[, 2]

  # Return
  list(mean = F_MU, sd = F_SD)
}

# Safeguard
prevent_too_large_mats <- function(fit, x, reduce, draws, verbose, force) {
  S <- determine_num_paramsets(fit, draws, reduce)
  J <- length(component_names(fit@model))
  N <- get_num_obs(fit)
  P <- if (is.null(x)) N else nrow(x)
  msg <- paste0(
    "Computations will require creating and handling ",
    P, " x ", N, " matrices ", S, " x ", J + 1,
    " times."
  )
  P_LIMIT <- 6000
  log_info(msg, verbose)
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
