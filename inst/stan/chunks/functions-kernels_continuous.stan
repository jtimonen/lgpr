// Multiplier matrix to enable variance masking
matrix STAN_kernel_base_var_mask(
  data vector x1,
  data vector x2,
  real steepness,
  data real[] vm_params)
{
  int n1 = num_elements(x1);
  int n2 = num_elements(x2);
  real a = steepness * vm_params[2];
  real r = inv(a)*logit(vm_params[1]);
  matrix[1, n1] s1 = to_matrix(STAN_var_mask(x1 - r, a));
  matrix[1, n2] s2 = to_matrix(STAN_var_mask(x2 - r, a));
  matrix[n1, n2] K = s1 * transpose(s2);
  return(K);
}

// Compute a nonstationary kernel matrix
matrix STAN_kernel_nonstationary(
  vector x1,
  vector x2,
  real alpha,
  real ell,
  real steepness,
  data real[,] vm_params)
{

  // Input warping
  int n1 = num_elements(x1);
  int n2 = num_elements(x2);
  real w1[n1] = to_array_1d(STAN_warp_input(x1, steepness));
  real w2[n2] = to_array_1d(STAN_warp_input(x2, steepness));
  matrix[n1, n2] K = cov_exp_quad(w1, w2, alpha, ell);
  
  // Variance masking
  int is_var_masked = size(vm_params);
  if(is_var_masked){
    K = K .* STAN_kernel_base_var_mask(x1, x2, steepness, vm_params[1]);
  }
  
  return(K);
}
