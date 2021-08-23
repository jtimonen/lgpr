# Create a kernel computer
create_kernel_computer <- function(model,
                                   stan_fit,
                                   x,
                                   x_is_data,
                                   reduce,
                                   draws,
                                   STREAM) {

  # Settings
  check_not_null(x)
  full_covariance <- FALSE
  input <- kernelcomp.create_input(model, stan_fit, x, reduce, draws)

  # Constant kernel computations and covariate-dependent inputs
  K_input <- kernelcomp.init(input, FALSE, FALSE, STREAM, FALSE)
  no_separate_output_points <- x_is_data
  if (no_separate_output_points) {
    x <- get_data(model)
    Ks_input <- K_input
  } else {
    Ks_input <- kernelcomp.init(input, TRUE, FALSE, STREAM, FALSE)
  }
  Kss_input <- kernelcomp.init(input, TRUE, TRUE, STREAM, !full_covariance)

  # Return
  new("KernelComputer",
    x = x,
    input = input,
    K_input = K_input,
    Ks_input = Ks_input,
    Kss_input = Kss_input,
    comp_names = component_names(model),
    full_covariance = full_covariance,
    STREAM = STREAM,
    no_separate_output_points = no_separate_output_points
  )
}

# Compute all constant kernel matrices of a model
# Input is a list returned by fp_input
kernel_const_all <- function(input, is_out1, is_out2, STREAM) {
  A1 <- if (is_out1) "Z_OUT" else "Z"
  A2 <- if (is_out2) "Z_OUT" else "Z"
  B1 <- if (is_out1) "X_mask_OUT" else "X_mask"
  B2 <- if (is_out2) "X_mask_OUT" else "X_mask"
  C1 <- if (is_out1) "N_OUT" else "N"
  C2 <- if (is_out2) "N_OUT" else "N"

  Z1 <- matrix_to_list(dollar(input, A1))
  Z2 <- matrix_to_list(dollar(input, A2))
  X1_mask <- matrix_to_list(dollar(input, B1))
  X2_mask <- matrix_to_list(dollar(input, B2))
  N1 <- dollar(input, C1)
  N2 <- dollar(input, C2)

  components <- matrix_to_list(dollar(input, "components"))
  Z_M <- dollar(input, "Z_M")
  STAN_kernel_const_all(
    N1, N2, Z1, Z2, X1_mask, X2_mask, Z_M, components, STREAM
  )
}

# Compute diagonals of all constant kernel matrices of a model,
# between output points
# Input is a list returned by fp_input
kernel_const_all_diag <- function(input, STREAM) {

  # Get input
  Z <- matrix_to_list(dollar(input, "Z_OUT"))
  X_mask <- matrix_to_list(dollar(input, "X_mask_OUT"))
  n <- dollar(input, "N_OUT")
  components <- matrix_to_list(dollar(input, "components"))

  # Call Stan function
  STAN_kernel_const_all_diag(n, Z, X_mask, components, STREAM)
}

# Compute all kernel matrices given draw idx
kernel_all <- function(init, input, idx, STREAM) {
  components <- matrix_to_list(dollar(input, "components"))
  vm_params <- dollar(input, "vm_params")
  teff_zero <- matrix_to_list(dollar(input, "teff_zero"))

  # Get parameters in correct format
  alpha_idx <- dollar(input, "d_alpha")[idx, ]
  ell_idx <- dollar(input, "d_ell")[idx, ]
  wrp_idx <- dollar(input, "d_wrp")[idx, ]
  beta_idx <- list(dollar(input, "d_beta")[idx, , ])
  teff_idx <- list(dollar(input, "d_teff")[idx, , ])

  # Get covariate input in correct format
  K_const <- dollar(init, "K_const")
  N1 <- dollar(init, "N1")
  N2 <- dollar(init, "N2")
  X1 <- dollar(init, "X1")
  X2 <- dollar(init, "X2")
  X_scale <- dollar(init, "X_scale")
  beta_idx1 <- matrix_to_list(dollar(init, "BETA_IDX1"))
  beta_idx2 <- matrix_to_list(dollar(init, "BETA_IDX2"))
  teff_idx1 <- matrix_to_list(dollar(init, "TEFF_IDX1"))
  teff_idx2 <- matrix_to_list(dollar(init, "TEFF_IDX2"))

  # Compute kernels for each component (a list with length num_comps)
  K_all <- STAN_kernel_all(
    N1, N2, K_const, components,
    X1, X2, X_scale,
    alpha_idx, ell_idx, wrp_idx, beta_idx, teff_idx, vm_params,
    beta_idx1, beta_idx2, teff_idx1, teff_idx2, teff_zero, STREAM
  )
  return(K_all)
}


# Compute all kernel matrices' diagonals given draw idx
kernel_all_diag <- function(init, input, idx, STREAM) {
  components <- matrix_to_list(dollar(input, "components"))
  vm_params <- dollar(input, "vm_params")
  teff_zero <- matrix_to_list(dollar(input, "teff_zero"))

  # Get parameters in correct format
  alpha_idx <- dollar(input, "d_alpha")[idx, ]
  wrp_idx <- dollar(input, "d_wrp")[idx, ]
  beta_idx <- list(dollar(input, "d_beta")[idx, , ])
  teff_idx <- list(dollar(input, "d_teff")[idx, , ])

  # Get covariate input in correct format
  K_const_diag <- dollar(init, "K_const")
  N <- dollar(init, "N1")
  X <- dollar(init, "X1")
  X_scale <- dollar(init, "X_scale")
  BETA_IDX <- matrix_to_list(dollar(init, "BETA_IDX1"))
  TEFF_IDX <- matrix_to_list(dollar(init, "TEFF_IDX1"))

  # Compute kernel diagonals for each component (a list with length J)
  K_all_diag <- STAN_kernel_all_diag(
    N, K_const_diag, components,
    X, X_scale, alpha_idx, wrp_idx, beta_idx, teff_idx,
    vm_params, BETA_IDX, TEFF_IDX, teff_zero, STREAM
  )
  return(K_all_diag)
}


# Initialize kernel matrix computations
kernelcomp.init <- function(input, is_out1, is_out2, STREAM, diag) {

  # Compute constant kernel matrices
  if (diag) {
    K_const <- kernel_const_all_diag(input, STREAM)
  } else {
    K_const <- kernel_const_all(input, is_out1, is_out2, STREAM)
  }

  # Covariate-input field names
  field_name <- function(is_out, base_name) {
    if (is_out) paste0(base_name, "_OUT") else base_name
  }
  A1 <- field_name(is_out1, "X")
  A2 <- field_name(is_out2, "X")
  B1 <- field_name(is_out1, "N")
  B2 <- field_name(is_out2, "N")
  C1 <- field_name(is_out1, "BETA_IDX")
  C2 <- field_name(is_out2, "BETA_IDX")
  D1 <- field_name(is_out1, "TEFF_IDX")
  D2 <- field_name(is_out2, "TEFF_IDX")

  # Get fields (possibly in list format)
  list(
    K_const = K_const,
    X1 = matrix_to_list(dollar(input, A1)),
    X2 = matrix_to_list(dollar(input, A2)),
    X_scale = dollar(input, "X_scale"),
    N1 = dollar(input, B1),
    N2 = dollar(input, B2),
    BETA_IDX1 = dollar(input, C1),
    BETA_IDX2 = dollar(input, C2),
    TEFF_IDX1 = dollar(input, D1),
    TEFF_IDX2 = dollar(input, D2),
    is_diag = diag
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
      N_OUT = dollar(si, "N"),
      X_OUT = dollar(si, "X"),
      X_mask_OUT = dollar(si, "X_mask"),
      Z_OUT = dollar(si, "Z"),
      BETA_IDX_OUT = dollar(si, "BETA_IDX"),
      TEFF_IDX_OUT = dollar(si, "TEFF_IDX")
    )
  } else {
    form <- model@model_formula
    covs_stan <- standata_base.covariates(x, form)
    comps_stan <- standata_base.components(form, covs_stan)
    expanding <- standata_base.expanding(covs_stan, comps_stan)
    out <- list(
      N_OUT = nrow(x),
      X_OUT = dollar(covs_stan, "X"),
      X_mask_OUT = dollar(covs_stan, "X_mask"),
      Z_OUT = dollar(covs_stan, "Z"),
      BETA_IDX_OUT = dollar(expanding, "BETA_IDX"),
      TEFF_IDX_OUT = dollar(expanding, "TEFF_IDX")
    )
  }
  return(out)
}

# parameter draws input (common fields)
kernelcomp.input_draws <- function(model, stan_fit, reduce, draws) {

  # Get dimensions
  S <- determine_num_paramsets(stan_fit, draws, reduce)
  si <- get_stan_input(model)
  J <- dollar(si, "J")
  num_ell <- dollar(si, "num_ell")
  num_wrp <- dollar(si, "num_wrp")
  DH <- dollar(si, "num_het") > 0
  DU <- dollar(si, "num_unc") > 0
  num_beta <- dollar(si, "num_beta")
  num_teff <- dollar(si, "num_teff")

  # Get draws
  list(
    num_paramsets = S,
    d_alpha = get_draw_arr(stan_fit, draws, reduce, "alpha", S, J),
    d_ell = get_draw_arr(stan_fit, draws, reduce, "ell", S, num_ell),
    d_wrp = get_draw_arr(stan_fit, draws, reduce, "wrp", S, num_wrp),
    d_beta = get_draw_arr_vec(stan_fit, draws, reduce, "beta", S, DH, num_beta),
    d_teff = get_draw_arr_vec(stan_fit, draws, reduce, "teff", S, DU, num_teff)
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
