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
kernel_all <- function(input, K_const, is_out1, is_out2, STREAM) {

  # List format
  A1 <- if (is_out1) "x_cont_OUT" else "x_cont"
  A2 <- if (is_out2) "x_cont_OUT" else "x_cont"
  B1 <- if (is_out1) "x_cont_unnorm_OUT" else "x_cont_unnorm"
  B2 <- if (is_out2) "x_cont_unnorm_OUT" else "x_cont_unnorm"
  C1 <- if (is_out1) "num_OUT" else "num_obs"
  C2 <- if (is_out2) "num_OUT" else "num_obs"
  D1 <- if (is_out1) "idx_expand_OUT" else "idx_expand"
  D2 <- if (is_out2) "idx_expand_OUT" else "idx_expand"

  x1 <- matrix_to_list(dollar(input, A1))
  x2 <- matrix_to_list(dollar(input, A2))
  x1_unnorm <- matrix_to_list(dollar(input, B1))
  x2_unnorm <- matrix_to_list(dollar(input, B2))
  n1 <- dollar(input, C1)
  n2 <- dollar(input, C2)
  idx1_expand <- dollar(input, D1) # list?
  idx2_expand <- dollar(input, D2) # list?

  components <- matrix_to_list(dollar(input, "components"))
  vm_params <- dollar(input, "vm_params")
  teff_zero <- matrix_to_list(dollar(input, "teff_zero"))

  # TODO: params, list?

  # Return
  STAN_kernel_all(
    n1, n2, K_const, components,
    x1, x2, x1_unnorm, x2_unnorm,
    alpha, ell, wrp, beta, teff, vm_params,
    idx1_expand, idx2_expand, teff_zero, STREAM
  )
}


# Kernel matrices for function posterior computation
kernels_fpost <- function(fit, x, reduce, draws, verbose, STREAM) {
  if (verbose) cat("Creating inputs...\n")
  input <- fp_input(fit, x, reduce, draws)

  if (verbose) cat("Computing constant kernel matrices...\n")
  K_const <- kernel_const_all(input, FALSE, FALSE, STREAM)
  Ks_const <- kernel_const_all(input, TRUE, FALSE, STREAM)
  Kss_const <- kernel_const_all(input, TRUE, TRUE, STREAM)

  if (verbose) cat("Computing parameter-dependent kernel matrices...\n")
  K <- kernel_all(input, K_const, FALSE, FALSE, STREAM)
  Ks <- kernel_all(input, Ks_const, TRUE, FALSE, STREAM)
  Kss <- kernel_all(input, Kss_const, TRUE, TRUE, STREAM)

  # Return
  list(
    K = K,
    Ks = Ks,
    Kss = Kss
  )
}
