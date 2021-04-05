#' Kernel matrices for function posterior computation
#'
#' @inheritParams posterior_f
kernels_fpost <- function(fit,
                          x,
                          reduce,
                          draws,
                          verbose,
                          debug_dims,
                          STREAM = get_stream()) {
  vrb <- verbose
  dd <- debug_dims
  if (vrb) cat("Creating inputs...\n")
  input <- kernels_fpost.create_input(fit, x, reduce, draws)

  if (vrb) cat("Creating constant kernel matrices...\n")
  K_const <- kernel_const_all(input, FALSE, FALSE, STREAM)
  Ks_const <- kernel_const_all(input, TRUE, FALSE, STREAM)
  Kss_const <- kernel_const_all(input, TRUE, TRUE, STREAM)

  # Final parameter-dependent kernel matrices
  if (verbose) cat("Creating final kernel matrices (data vs. data)...\n")
  K <- kernel_all(input, K_const, FALSE, FALSE, vrb, dd, STREAM)
  if (is.null(x)) {
    Ks <- K
    Kss <- K
  } else {
    if (verbose) cat("Creating final kernel matrices (out vs. data)...\n")
    Ks <- kernel_all(input, Ks_const, TRUE, FALSE, vrb, dd, STREAM)
    if (verbose) cat("Creating final kernel matrices (out vs. out)...\n")
    Kss <- kernel_all(input, Kss_const, TRUE, TRUE, vrb, dd, STREAM)
  }

  # Return
  comp_names <- component_names(fit@model)
  list(
    K = K,
    Ks = Ks,
    Kss = Kss,
    comp_names = comp_names
  )
}

# Compute all constant kernel matrices of a model
# Input is a list returned by fp_input
kernel_const_all <- function(input, is_out1, is_out2, STREAM) {
  A1 <- if (is_out1) "x_cat_OUT" else "x_cat"
  A2 <- if (is_out2) "x_cat_OUT" else "x_cat"
  B1 <- if (is_out1) "x_cont_mask_OUT" else "x_cont_mask"
  B2 <- if (is_out2) "x_cont_mask_OUT" else "x_cont_mask"
  C1 <- if (is_out1) "num_OUT" else "num_obs"
  C2 <- if (is_out2) "num_OUT" else "num_obs"

  x1_cat <- matrix_to_list(dollar(input, A1))
  x2_cat <- matrix_to_list(dollar(input, A2))
  x1_cont_mask <- matrix_to_list(dollar(input, B1))
  x2_cont_mask <- matrix_to_list(dollar(input, B2))
  n1 <- dollar(input, C1)
  n2 <- dollar(input, C2)

  components <- matrix_to_list(dollar(input, "components"))
  x_cat_num_levels <- dollar(input, "x_cat_num_levels")
  STAN_kernel_const_all(
    n1, n2, x1_cat, x2_cat, x1_cont_mask, x2_cont_mask,
    x_cat_num_levels, components, STREAM
  )
}



# Compute all kernel matrices of a model
kernel_all <- function(input, K_const, is_out1, is_out2,
                       verbose, debug_dims, STREAM) {

  # Field names
  field_name <- function(is_out, base_name) {
    if (is_out) paste0(base_name, "_OUT") else base_name
  }
  A1 <- field_name(is_out1, "x_cont")
  A2 <- field_name(is_out2, "x_cont")
  B1 <- field_name(is_out1, "x_cont_unnorm")
  B2 <- field_name(is_out2, "x_cont_unnorm")
  C1 <- if (is_out1) "num_OUT" else "num_obs"
  C2 <- if (is_out2) "num_OUT" else "num_obs"
  D1 <- field_name(is_out1, "idx_expand")
  D2 <- field_name(is_out2, "idx_expand")

  # Get fields (possibly in list format)
  x1 <- matrix_to_list(dollar(input, A1))
  x2 <- matrix_to_list(dollar(input, A2))
  x1_unnorm <- matrix_to_list(dollar(input, B1))
  x2_unnorm <- matrix_to_list(dollar(input, B2))
  n1 <- dollar(input, C1)
  n2 <- dollar(input, C2)
  idx1_expand <- dollar(input, D1)
  idx2_expand <- dollar(input, D2)
  components <- matrix_to_list(dollar(input, "components"))
  vm_params <- dollar(input, "vm_params")
  teff_zero <- matrix_to_list(dollar(input, "teff_zero"))

  # Kernel hyper-parameters
  alpha <- dollar(input, "d_alpha")
  ell <- dollar(input, "d_ell")
  wrp <- dollar(input, "d_wrp")
  beta <- dollar(input, "d_beta") # has shape (S, num_heter > 1, num_bt)
  teff <- dollar(input, "d_teff") # has shape (S, num_uncrt > 1, num_bt)

  # Setup
  S <- dollar(input, "num_paramsets")
  K_out <- list()

  # Print dimensions
  if (debug_dims) {
    print_arr_dim(alpha)
    print_arr_dim(ell)
    print_arr_dim(wrp)
    print_arr_dim(beta)
    print_arr_dim(teff)
    print_arr_dim(K_out)
  }

  # Loop through parameter sets
  progbar <- verbose && S > 1
  pb <- progbar_setup(S)
  idx_print <- dollar(pb, "idx_print")
  if (progbar) cat(dollar(pb, "header"))

  for (idx in seq_len(S)) {

    # Get parameters in correct format
    alpha_idx <- alpha[idx, ]
    ell_idx <- ell[idx, ]
    wrp_idx <- wrp[idx, ]
    beta_idx <- list(beta[idx, , ])
    teff_idx <- list(teff[idx, , ])

    # Compute kernels for each component
    K_out[[idx]] <- STAN_kernel_all(
      n1, n2, K_const, components,
      x1, x2, x1_unnorm, x2_unnorm,
      alpha_idx, ell_idx, wrp_idx, beta_idx, teff_idx,
      vm_params, idx1_expand, idx2_expand, teff_zero, STREAM
    )

    # Store result and update progress
    if (progbar) progbar_print(idx, idx_print)
  }
  if (progbar) cat("\n")
  return(K_out)
}

# Create a list of things that will be used as input to the wrapped Stan
# kernel computation functions (after some formatting)
kernels_fpost.create_input <- function(fit, x, reduce, draws) {
  si <- get_stan_input(fit) # common
  si_x_OUT <- fp_input_x(fit, x) # output points (covariate values)
  si_draws <- fp_input_draws(fit, reduce, draws) # parameter values
  c(si, si_x_OUT, si_draws)
}

# covariate input
fp_input_x <- function(fit, x) {
  si <- get_stan_input(fit)
  if (is.null(x)) {
    out <- list(
      num_OUT = dollar(si, "num_obs"),
      x_cont_OUT = dollar(si, "x_cont"),
      x_cont_unnorm_OUT = dollar(si, "x_cont_unnorm"),
      x_cont_mask_OUT = dollar(si, "x_cont_mask"),
      x_cat_OUT = dollar(si, "x_cat"),
      idx_expand_OUT = dollar(si, "idx_expand")
    )
  } else {
    m <- get_model(fit)
    x_names <- unique(rhs_variables(m@model_formula@terms))
    check_df_with(x, x_names)
    x_cont_scl <- dollar(m@var_scalings, "x_cont")
    covariates <- stan_data_covariates(x, x_names, x_cont_scl)
    covs_stan <- dollar(covariates, "to_stan")
    comp_info <- get_component_encoding(m)
    expanding <- stan_data_expanding(covs_stan, comp_info)
    si <- c(covs_stan, dollar(expanding, "to_stan"))
    out <- list(
      num_OUT = nrow(x),
      x_cont_OUT = dollar(si, "x_cont"),
      x_cont_unnorm_OUT = dollar(si, "x_cont_unnorm"),
      x_cont_mask_OUT = dollar(si, "x_cont_mask"),
      x_cat_OUT = dollar(si, "x_cat"),
      idx_expand_OUT = dollar(si, "idx_expand")
    )
  }
  return(out)
}

# parameter draws input
fp_input_draws <- function(fit, reduce, draws) {

  # Get param sets
  d_common <- fp_input_draws.common(fit, reduce, draws)
  if (is_f_sampled(fit)) {
    d_add <- fp_input_draws.latent(fit, reduce, draws)
  } else {
    d_add <- fp_input_draws.marginal(fit, reduce, draws)
  }
  c(d_common, d_add)
}

# parameter draws input (marginal gp)
fp_input_draws.marginal <- function(fit, reduce, draws) {
  S <- determine_num_paramsets(fit, draws, reduce)
  list(
    d_sigma = get_draw_arr(fit, draws, reduce, "sigma", S, 1)
  )
}

# parameter draws input (latent gp)
fp_input_draws.latent <- function(fit, reduce, draws) {
  S <- determine_num_paramsets(fit, draws, reduce)
  si <- get_stan_input(fit)
  LH <- dollar(si, "obs_model")
  num_sigma <- as.numeric(LH == 1)
  num_phi <- as.numeric(LH == 3)
  num_gamma <- as.numeric(LH == 5)

  # Get draws
  list(
    num_sigma = num_sigma,
    num_phi = num_phi,
    num_gamma = num_gamma,
    d_f_latent = get_draws(fit, pars = "f_latent"), # S x (num_comps*num_obs)
    d_sigma = get_draw_arr(fit, draws, reduce, "sigma", S, num_sigma),
    d_phi = get_draw_arr(fit, draws, reduce, "phi", S, num_phi),
    d_gamma = get_draw_arr(fit, draws, reduce, "gamma", S, num_gamma)
  )
}

# parameter draws input (common fields)
fp_input_draws.common <- function(fit, reduce, draws) {

  # Get dimensions
  S <- determine_num_paramsets(fit, draws, reduce)
  si <- get_stan_input(fit)
  num_comps <- dollar(si, "num_comps")
  num_ell <- dollar(si, "num_ell")
  num_ns <- dollar(si, "num_ns")
  UNCRT <- dollar(si, "num_uncrt") > 0
  HETER <- dollar(si, "num_heter") > 0
  num_bt <- dollar(si, "num_bt")

  # Get draws
  list(
    num_paramsets = S,
    d_alpha = get_draw_arr(fit, draws, reduce, "alpha", S, num_comps),
    d_ell = get_draw_arr(fit, draws, reduce, "ell", S, num_ell),
    d_wrp = get_draw_arr(fit, draws, reduce, "wrp", S, num_ns),
    d_beta = get_draw_arr_vec(fit, draws, reduce, "beta", S, HETER, num_bt),
    d_teff = get_draw_arr_vec(fit, draws, reduce, "teff", S, UNCRT, num_bt)
  )
}

# Get an array of draws formatted suitably for Stan input
get_draw_arr <- function(fit, draws, reduce, par_name, S, D) {
  out <- array(0.0, dim = c(S, D))
  if (D > 0) {
    out <- get_draws(fit, draws = draws, reduce = reduce, pars = par_name)
  }
  return(out)
}

# Get an array of vector draws formatted suitably for Stan input
get_draw_arr_vec <- function(fit, draws, reduce, par_name, S, B, V) {
  D <- as.numeric(B)
  out <- array(0.0, dim = c(S, D, V))
  if (B) {
    tmp <- get_draws(fit, draws = draws, reduce = reduce, pars = par_name)
    # tmp has dim c(S, V)
    out[, 1, ] <- tmp
  }
  return(out)
}
