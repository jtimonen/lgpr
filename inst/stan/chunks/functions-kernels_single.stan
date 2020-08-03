// Compute one fixed kernel matrix
// - Does not depend on parameters and therefore this function
//   never needs to be evaluated during sampling
matrix STAN_kernel_fixed(
  data int[] x1,
  data int[] x2,
  data int num_levels,
  data int ctype,
  data int ktype)
{
  int n1 = size(x1);
  int n2 = size(x2);
  matrix[n1, n2] K;
  if(ctype==0 || ctype==2){
    K = STAN_kernel_discrete(x1, x2, ktype, num_levels);
  }else if(ctype==3){
    K = STAN_kernel_base_disease_mask(x1, x2); 
  }else{
    reject("STAN_kernel_fixed: <ctype> should be 0, 2, or 3!")
  }
  return(K);
}

// Compute one stationary kernel matrix
// - This function needs to be evaluated repeatedly during sampling
matrix STAN_kernel_stationary(
  matrix K_fixed,
  data vector x1,
  data vector x2,
  data int ctype,
  real alpha,
  real ell)
{
  int n1 = num_elements(x1);
  int n2 = num_elements(x2);
  matrix[n1, n2] K;
  if(ctype == 0){
    K = square(alpha) * K_fixed;
  }else if(ctype==1 || ctype==2){
    K = cov_exp_quad(to_array_1d(x1), to_array_1d(x2), alpha, ell);
    if(ctype == 2){
      K = K .* K_fixed;
    }
  }else{
    reject("STAN_kernel_stationary: <ctype> must be 0, 1, or 2!")
  }
  return(K);
}


// Compute the diseaase kernel matrix
// - This function needs to be evaluated repeatedly during sampling
matrix STAN_kernel_disease(
  matrix K_fixed,
  data vector x1,
  data vector x2,
  real alpha,
  real ell,
  real wrp,
  vector[] beta,
  vector[] teff,
  data real[,] vm_params,
  data int[,] idx1_expand,
  data int[,] idx2_expand,
  data vector[] teff_obs)
{
  
  // Options and variable declarations
  int is_heter = size(beta);
  int is_uncrt = size(teff);
  int n1 = num_elements(x1);
  int n2 = num_elements(x2); 
  matrix[n1, n2] K = K_fixed;
  
  // Handle uncertainty in disease-related age
  vector[n1] t1;
  vector[n2] t2;
  if(is_uncrt){
    t1 = STAN_edit_dis_age(x1, idx1_expand[1], teff_obs[1], teff[1]);
    t2 = STAN_edit_dis_age(x2, idx2_expand[1], teff_obs[1], teff[1]);
  }else{
    t1 = x1;
    t2 = x2;
  }
  
  // Compute kernel
  K = K .* STAN_kernel_nonstationary(t1, t2, alpha, ell, wrp, vm_params);
  
  // Heterogeneity
  if(is_heter){
    vector[n1] b1 = beta[1][idx1_expand[1]];
    vector[n2] b2 = beta[1][idx2_expand[1]];
    matrix[n1, n2] K_beta = to_matrix(b1) * transpose(to_matrix(b2));
    K = K .* K_beta;
  }
  
  return(K);
}
