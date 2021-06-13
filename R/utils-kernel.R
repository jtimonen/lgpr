# Create a kernel computer
create_kernel_computer <- function(model,
                                   stan_fit,
                                   x,
                                   reduce,
                                   draws,
                                   STREAM) {

  # Settings
  if (!is.null(draws)) reduce <- NULL
  input <- kernelcomp.create_input(model, stan_fit, x, reduce, draws)

  # Constant kernel computations and covariate-dependent inputs
  K_init <- kernelcomp.init(input, FALSE, FALSE, STREAM)
  if (is.null(x)) {
    Ks_init <- NULL
    Kss_init <- NULL
    x <- get_data(model)
    P <- dollar(K_init, "n1") # number of output points
  } else {
    Ks_init <- kernelcomp.init(input, TRUE, FALSE, STREAM)
    Kss_init <- kernelcomp.init(input, TRUE, TRUE, STREAM)
    P <- dollar(Ks_init, "n1") # number of output points
  }
  N <- dollar(K_init, "n1")
  S <- dollar(input, "num_paramsets")
  J <- get_num_comps(model)
  comp_names <- component_names(model)
  init <- list(
    K_init = K_init, Ks_init = Ks_init, Kss_init = Kss_init,
    P = P, N = N, S = S, J = J, comp_names = comp_names
  )

  # Get kernel parameter draws
  param_draws <- list(
    alpha = dollar(input, "d_alpha"),
    ell = dollar(input, "d_ell"),
    wrp = dollar(input, "d_wrp"),
    beta = dollar(input, "d_beta"), # has shape (S, num_heter > 1, num_bt)
    teff = dollar(input, "d_teff") # has shape (S, num_uncrt > 1, num_bt)
  )

  # Return
  new("KernelComputer",
    init = init,
    input = input,
    param_draws = param_draws,
    STREAM = STREAM
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

# Compute all constant kernel matrices of a model given draw idx
kernel_all <- function(init, input, param_draws, idx, STREAM) {
  components <- matrix_to_list(dollar(input, "components"))
  vm_params <- dollar(input, "vm_params")
  teff_zero <- matrix_to_list(dollar(input, "teff_zero"))

  # Get parameters in correct format
  alpha_idx <- dollar(param_draws, "alpha")[idx, ]
  ell_idx <- dollar(param_draws, "ell")[idx, ]
  wrp_idx <- dollar(param_draws, "wrp")[idx, ]
  beta_idx <- list(dollar(param_draws, "beta")[idx, , ])
  teff_idx <- list(dollar(param_draws, "teff")[idx, , ])

  # Get covariate input in correct format
  K_const <- dollar(init, "K_const")
  n1 <- dollar(init, "n1")
  n2 <- dollar(init, "n2")
  x1 <- dollar(init, "x1")
  x2 <- dollar(init, "x2")
  x1_unnorm <- dollar(init, "x1_unnorm")
  x2_unnorm <- dollar(init, "x2_unnorm")
  idx1_expand <- dollar(init, "idx1_expand")
  idx2_expand <- dollar(init, "idx2_expand")

  # Compute kernels for each component (a list with length num_comps)
  K_all <- STAN_kernel_all(
    n1, n2, K_const, components,
    x1, x2, x1_unnorm, x2_unnorm,
    alpha_idx, ell_idx, wrp_idx, beta_idx, teff_idx,
    vm_params, idx1_expand, idx2_expand, teff_zero, STREAM
  )
  return(K_all)
}

# Initialize kernel matrix computations
kernelcomp.init <- function(input, is_out1, is_out2, STREAM) {

  # Compute constant kernel matrices
  K_const <- kernel_const_all(input, is_out1, is_out2, STREAM)

  # Covariate-input field names
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

  # Return
  list(
    K_const = K_const,
    x1 = x1,
    x2 = x2,
    x1_unnorm = x1_unnorm,
    x2_unnorm = x2_unnorm,
    n1 = n1,
    n2 = n2,
    idx1_expand = idx1_expand,
    idx2_expand = idx2_expand
  )
}

# Create a list of things that will be used as input to the wrapped Stan
# kernel computation functions (after some formatting)
kernelcomp.create_input <- function(model, stan_fit, x, reduce, draws) {
  si <- get_stan_input(model) # common
  si_x_OUT <- kernelcomp.input_x(model, x) # output points (covariate values)
  si_draws <- kernelcomp.input_draws(model, stan_fit, reduce, draws) # params
  c(si, si_x_OUT, si_draws)
}

# covariate input
kernelcomp.input_x <- function(model, x) {
  si <- get_stan_input(model)
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
    m <- model
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

# parameter draws input (common fields)
kernelcomp.input_draws <- function(model, stan_fit, reduce, draws) {

  # Get dimensions
  S <- determine_num_paramsets(stan_fit, draws, reduce)
  si <- get_stan_input(model)
  num_comps <- dollar(si, "num_comps")
  num_ell <- dollar(si, "num_ell")
  num_ns <- dollar(si, "num_ns")
  UNC <- dollar(si, "num_uncrt") > 0
  HET <- dollar(si, "num_heter") > 0
  num_bt <- dollar(si, "num_bt")

  # Get draws
  list(
    num_paramsets = S,
    d_alpha = get_draw_arr(stan_fit, draws, reduce, "alpha", S, num_comps),
    d_ell = get_draw_arr(stan_fit, draws, reduce, "ell", S, num_ell),
    d_wrp = get_draw_arr(stan_fit, draws, reduce, "wrp", S, num_ns),
    d_beta = get_draw_arr_vec(stan_fit, draws, reduce, "beta", S, HET, num_bt),
    d_teff = get_draw_arr_vec(stan_fit, draws, reduce, "teff", S, UNC, num_bt)
  )
}

# Get an array of draws formatted suitably for Stan input
get_draw_arr <- function(stan_fit, draws, reduce, par_name, S, D) {
  out <- array(0.0, dim = c(S, D))
  if (D > 0) {
    out <- get_draws(stan_fit, draws = draws, reduce = reduce, pars = par_name)
  }
  return(out)
}

# Get an array of vector draws formatted suitably for Stan input
get_draw_arr_vec <- function(stan_fit, draws, reduce, par_name, S, B, V) {
  D <- as.numeric(B)
  out <- array(0.0, dim = c(S, D, V))
  if (B) {
    tmp <- get_draws(stan_fit, draws = draws, reduce = reduce, pars = par_name)
    # tmp has dim c(S, V)
    out[, 1, ] <- tmp
  }
  return(out)
}
