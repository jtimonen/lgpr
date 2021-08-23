# Computes the analytic conditional posterior distributions of each model
# component, and their sum.
posterior_f_gaussian <- function(fit, x, x_is_data,
                                 reduce, draws, verbose, STREAM) {
  kc <- create_kernel_computer(fit, x, x_is_data, reduce, draws, STREAM)
  y <- get_y(fit, original = FALSE) # normalized y
  d_sigma <- get_draws(fit, pars = "sigma[1]", reduce = reduce, draws = draws)
  sigma2 <- as.vector(d_sigma^2)
  fp_gaussian(kc, sigma2, y, verbose)
}

#  Extract function samples and compute their sum. If \code{x_is_data} is not
#  \code{TRUE}, the function draws are extrapolated to the points \code{x}
#  using kernel regression.
posterior_f_latent <- function(fit, x, x_is_data,
                               reduce, draws, verbose, STREAM) {
  fp_at_data <- get_f_draws(fit, reduce = reduce, draws = draws)
  if (!x_is_data) {
    kc <- create_kernel_computer(fit, x, x_is_data, reduce, draws, STREAM)
    out <- fp_extrapolate(kc, fp_at_data, verbose)
  } else {
    out <- fp_at_data
  }

  return(out)
}

# Analytic function posteriors
fp_gaussian <- function(kc, sigma2, y, verbose) {

  # Extract info
  input <- kc@input # shared input
  S <- num_paramsets(kc) # number of parameter sets
  P <- num_evalpoints(kc) # number of output points
  J <- num_components(kc) # number of components
  comp_names <- component_names(kc)
  delta <- dollar(input, "delta")

  # Create output arrays
  f_comp_mean <- array(0.0, c(S, P, J))
  f_comp_std <- array(0.0, c(S, P, J))
  f_mean <- array(0.0, c(S, P))
  f_std <- array(0.0, c(S, P))

  # Setup
  progbar <- verbose && S > 1
  pb <- progbar_setup(L = S)
  hdr <- dollar(pb, "header")
  idx_print <- dollar(pb, "idx_print")
  log_progress("Computing analytic function posteriors...", verbose)
  log_progress(hdr, progbar)

  # Loop through parameter sets
  for (idx in seq_len(S)) {

    # Compute full kernel matrices for one parameter set
    K_i <- kernel_all(kc@K_input, kc@input, idx, kc@STREAM)
    if (kc@no_separate_output_points) {
      Ks_i <- K_i
    } else {
      Ks_i <- kernel_all(kc@Ks_input, kc@input, idx, kc@STREAM)
    }
    Kss_i_diag <- kernel_all_diag(kc@Kss_input, kc@input, idx, kc@STREAM)

    # Perform computations for one parameter set
    fp_i <- fp_gaussian.compute(K_i, Ks_i, Kss_i_diag, sigma2[idx], delta, y)
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
  f_comp_mean <- aperm(f_comp_mean, c(3, 1, 2)) # dim (S, P, J) -> (J, S, P)
  f_comp_std <- aperm(f_comp_std, c(3, 1, 2)) # dim (S, P, J) -> (J, S, P)
  list(
    f_comp_mean = arr3_to_list(f_comp_mean, comp_names), # list with len J
    f_comp_std = arr3_to_list(f_comp_std, comp_names), # list with len J
    f_mean = f_mean, # dim (S, P)
    f_std = f_std, # dim (S, P)
    sigma2 = sigma2,
    x = kc@x
  )
}

# Compute componentwise and total function posteriors
fp_gaussian.compute <- function(K, Ks, Kss_diag, sigma2, delta, y) {

  # Helper function for linear algebra
  # See e.g. http://www.gaussianprocess.org/gpml/chapters/RW.pdf Algorithm 2.1
  gp_posterior_helper <- function(Ly, Ks, Kss_diag, v) {
    P <- length(Kss_diag)
    A <- t(forwardsolve(Ly, t(Ks))) # A is (P x N)
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
fp_extrapolate <- function(kc, fp_at_data, verbose) {

  # Extract info
  input <- kc@input # shared input
  S <- num_paramsets(kc) # number of parameter sets
  P <- num_evalpoints(kc) # number of output points
  J <- num_components(kc) # number of components
  comp_names <- component_names(kc)
  delta <- dollar(input, "delta")
  fp_comp_draws <- dollar(fp_at_data, "comp")
  take_row <- function(A, idx) A[idx, ]

  # Create output arrays
  f_ext_comp <- array(0, c(S, P, J))
  f_ext <- array(0, c(S, P))

  # Setup
  progbar <- verbose && S > 1
  pb <- progbar_setup(L = S)
  hdr <- dollar(pb, "header")
  idx_print <- dollar(pb, "idx_print")
  log_progress("Extrapolating function posterior draws...", verbose)
  log_progress(hdr, progbar)

  # Loop through parameter sets
  for (idx in seq_len(S)) {

    # Perform computations for one parameter set
    K_i <- kernel_all(kc@K_input, kc@input, idx, kc@STREAM)
    if (kc@no_separate_output_points) {
      Ks_i <- K_i
    } else {
      Ks_i <- kernel_all(kc@Ks_input, kc@input, idx, kc@STREAM)
    }
    fp_comp_i <- sapply(fp_comp_draws, take_row, idx = idx) # dim = (N, J)
    fp_i <- fp_extrapolate.compute(K_i, Ks_i, delta, fp_comp_i) # dim = (P, J)

    # Store result and update progress
    f_ext_comp[idx, , ] <- fp_i
    f_ext[idx, ] <- rowSums(fp_i)

    if (progbar) progbar_print(idx, idx_print)
  }
  log_progress(" ", progbar)

  # Return
  f_ext_comp <- aperm(f_ext_comp, c(3, 1, 2)) # dim (S, P, J) -> (J, S, P)
  out <- list(
    comp = arr3_to_list(f_ext_comp, comp_names), # list with len J
    sum = f_ext, # dim (S, P),
    extapolated = TRUE,
    x = kc@x
  )
  return(out)
}

# Extrapolated componentwise and total function posterior draws,
# (for one MCMC draw)
fp_extrapolate.compute <- function(K, Ks, delta, fp_comp) {
  P <- dim(Ks[[1]])[1]
  N <- dim(fp_comp)[1]
  J <- dim(fp_comp)[2]
  out <- matrix(0, P, J)
  DELTA <- delta * diag(N)

  # Loop through components
  for (j in seq_len(J)) {
    fj <- fp_comp[, j]
    k <- K[[j]] + DELTA
    ks <- Ks[[j]]
    fj_ext <- ks %*% solve(k, fj) # kernel regression
    out[, j] <- fj_ext
  }
  return(out)
}


# Gets draws of the total \code{f}
# Returns an array of shape \code{num_draws} x \code{num_obs}
get_f_draws <- function(fit, draws, reduce) {
  f_comp <- get_f_draws.comp(fit, draws, reduce)
  list(
    comp = f_comp,
    sum = STAN_matrix_array_sum(f_comp, get_stream()),
    extrapolated = FALSE,
    x = get_data(fit)
  )
}

# Get the draws of each component of \code{f}
# Returns a list of arrays of shape \code{num_draws} x \code{num_obs}
get_f_draws.comp <- function(fit, draws, reduce) {
  nams <- component_names(fit)
  D <- length(nams)
  fp <- get_draws(fit, pars = "f_latent", draws = draws, reduce = reduce)
  fp <- array_to_arraylist(fp, D)
  names(fp) <- nams
  return(fp)
}

# Safeguard
prevent_too_large_mats <- function(fit, x, reduce, draws, verbose, force) {
  S <- determine_num_paramsets(fit@stan_fit, draws, reduce)
  J <- length(component_names(fit@model))
  N <- get_num_obs(fit)
  P <- if (is.null(x)) N else nrow(x)

  P_LIMIT <- 6000 # __HARDCODED__
  if (P > P_LIMIT) {
    msg <- paste0(
      "Computations will require creating and handling ",
      P, " x ", N, " matrices ", S * (J + 1), " times."
    )
    log_info(msg, verbose)
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
