// Compute all fixed kernel matrices
//
// These do not depend on parameters and therefore this function
// never needs to be evaluated during sampling
matrix[] STAN_kernel_fixed_all(
    data int[,] x1, 
    data int[,] x2,
    data int[] num_levels,
    data int[,] components)
{
  int n1 = size(x1[1]);
  int n2 = size(x2[1]);
  int num_comps = size(components[1]);
  matrix[n1, n2] K_fixed[num_comps];
  for(j in 1:num_comps){
    int ctype = components[1][j];
    if(ctype != 1){
      int ktype = components[2][j];
      int idx = components[3][j];
      int x1_j[n1] = x1[idx];
      int x2_j[n1] = x2[idx];
      int n_levels = num_levels[idx];
      K_fixed[j] = STAN_kernel_fixed(x1_j, x2_j, n_levels, ctype, ktype); 
    }else{
      K_fixed[j] = rep_matrix(0, n1, n2);
    }
  }
  return(K_fixed);
}

// Compute all kernel matrices
//
// These depend on parameters and therefore this function
// needs to be evaluated over and over during sampling
matrix[] STAN_kernel_all(
    matrix[] K_fixed,
    data int[,] components,
    data vector[] x1,
    data vector[] x2,
    real[] alpha,
    real[] ell,
    real[] wrp,
    vector[] beta,
    vector[] teff,
    data real[,] vm_params,
    data int[,] idx1_expand,
    data int[,] idx2_expand,
    data vector[] teff_obs)
{
  int ell_idx = 0;
  real ell_j;
  real alpha_j;
  int n1 = num_elements(x1[1]);
  int n2 = num_elements(x2[1]);
  int num_comps = num_elements(components[1]);
  matrix[n1, n2] KX[num_comps];
  
  for(j in 1:num_comps){
    
    // Component type + covariate
    int ctype = components[1][j];
    int idx = components[4][j];
    
    // Decide correct lengthscale param or if this component needs it at all
    alpha_j = alpha[j];
    if(ctype > 0){
      ell_idx += 1;
      ell_j = ell[ell_idx];
    }else{
      ell_j = 0.0; // not needed
    }
    
    // Compute the kernel matrix
    if(ctype!=3){
      KX[j] = STAN_kernel_stationary(
        K_fixed[j], x1[idx], x2[idx], ctype, alpha_j, ell_j);
    }else{
      KX[j] = STAN_kernel_disease(
        K_fixed[j], x1[idx], x2[idx], alpha_j, ell_j, wrp[1],
        beta, teff, vm_params, idx1_expand, idx2_expand, teff_obs);
    }

  }
  return(KX);
}
