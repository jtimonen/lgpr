  // Categorical zero-sum kernel
  matrix STAN_kernel_zerosum(data int[] x1, data int[] x2, data int num_cat) {
    int n1 = size(x1); 
    int n2 = size(x2);
    matrix[n1, n2] K;
    for (i in 1:n1) {
      for (j in 1:n2) {
        if (x1[i] == x2[j]) {
          K[i,j] = 1;
        } else {
          K[i,j] = - inv(num_cat - 1); 
        }
      }
    }
    return(K);
  }
  
  // Categorical kernel
  matrix STAN_kernel_cat(data int[] x1, data int[] x2) {
    int n1 = size(x1);
    int n2 = size(x2);
    matrix[n1,n2] K;
    for (i in 1:n1) {
      for (j in 1:n2) {
        K[i,j] = (x1[i] == x2[j]);
      }
    }
    return(K);
  }
  
  // Binary mask kernel
  matrix STAN_kernel_bin(data int[] x1, data int[] x2) {
    int n1 = size(x1);
    int n2 = size(x2);
    matrix[n1,n2] K;
    for (i in 1:n1) {
      for (j in 1:n2) {
        K[i,j] = (x1[i] == 0) * (x2[j] == 0);
      }
    }
    return(K);
  }
  
  // Compute one constant kernel matrix. Does not depend on parameters and
  // therefore this function never needs to be evaluated during sampling.
  matrix STAN_kernel_const(data int[] z1, data int[] z2, 
    data int kernel_type,  data int num_cat) 
  {
    int n1 = num_elements(z1);
    int n2 = num_elements(z2);
    matrix[n1, n2] K;
    if (kernel_type == 1) {
      K = STAN_kernel_cat(z1, z2);
    } else if (kernel_type == 2) {
      K = STAN_kernel_zerosum(z1, z2, num_cat);
    } else if (kernel_type == 3){
      K = STAN_kernel_bin(z1, z2);
    } else {
      reject("kernel_type should be 1, 2 or 3, found =", kernel_type)
    }
    return(K);
  }
  
  // Compute all constant kernel matrices. These do not depend on parameters and
  // therefore this function never needs to be evaluated during sampling.
  matrix[] STAN_kernel_const_all(
    data int n1,           data int n2,
    data int[,] Z1,        data int[,] Z2,
    data int[,] X1_mask,   data int[,] X2_mask,
    data int[] num_levels, data int[,] components)
  {
    int J = size(components);
    matrix[n1, n2] K_const[J];
    for (j in 1:J) {
      matrix[n1, n2] Kj;
      int idx_cat = components[j, 1];
      int ker_cat = components[j, 2];
      int idx_cont = components[j, 3];
      
      // Compute mask kernel for continuous covariate
      if (idx_cont != 0) {
        // binary mask
        Kj = STAN_kernel_const(X1_mask[idx_cont], X2_mask[idx_cont], 3, 0);
      } else {
        Kj = rep_matrix(1, n1, n2);
      }
      
      // Compute kernel for categorical covariate
      if (ker_cat != 0) {
        int M = num_levels[idx_cat];
        Kj = Kj .* STAN_kernel_const(Z1[idx_cat], Z2[idx_cat], ker_cat, M);
      }
      K_const[j] = Kj;
    }
    return(K_const);
  }
  
  // Exponentiated quadratic kernel (with vector inputs)
  matrix STAN_kernel_eq(vector x1, vector x2, real alpha, real ell) {
    return(cov_exp_quad(to_array_1d(x1), to_array_1d(x2), alpha, ell));
  }
  
  // Multiplier matrix to enable variance masking
  matrix STAN_kernel_varmask(vector x1, vector x2, 
    real steepness, data real[] vm_params) 
  {
    real a = steepness * vm_params[2];
    real r = inv(a)*logit(vm_params[1]);
    return(
      to_matrix(to_matrix(STAN_var_mask(x1 - r, a))) *
      transpose(to_matrix(STAN_var_mask(x2 - r, a)))
    );
  }
  
  // Multiplier matrix to enable heterogeneous effects
  matrix STAN_kernel_beta(vector beta, int[] idx1_expand, int[] idx2_expand) {
    return(
      to_matrix(STAN_expand(sqrt(beta), idx1_expand)) *
      transpose(to_matrix(STAN_expand(sqrt(beta), idx2_expand)))
    );
  }
  
  /* 
    Compute all kernel matrices. These depend on parameters and
    therefore this function needs to be evaluated repeatedly during sampling.
  */
  matrix[] STAN_kernel_all(
    data int n1,
    data int n2,
    data matrix[] K_const,
    data int[,] components,
    data vector[] X1,
    data vector[] X2,
    data real[] X_scale,
    real[] alpha,
    real[] ell,
    real[] wrp,
    vector[] beta,
    vector[] teff,
    data real[] vm_params,
    data int[,] beta_idx1,
    data int[,] beta_idx2,
    data int idx_unc,
    data int[,] teff_idx1,
    data int[,] teff_idx2,
    data vector[] teff_zero)
  {
    int idx_ell = 0;
    int idx_wrp = 0;
    int idx_alpha = 0;
    int J = size(components);
    matrix[n1, n2] KX[J];
  
    // Loop through components
    for(j in 1:J){
      
      // 1. Initialize with constant part of the kernel matrix
      matrix[n1, n2] Kj = K_const[j];
      vector[n1] x1;
      vector[n2] x2;
  
      // 2. Get component properties
      int idx_cont = components[j, 3];
      int ker_cont = components[j, 4];
      int is_warped = components[j, 5];
      int is_heter = components[j, 6];
      
      // 3. Pick the possible continuous covariate of this component
      if(idx_cont != 0){
        x1 = X1[idx_cont];
        x2 = X2[idx_cont];
        
        // 3.1 Handle possible uncertainty in covariate
        if(idx_unc>0 && idx_unc==idx_cont) {
          x1 = STAN_edit_x_cont(x1, teff_idx1[1], teff_zero[1], teff[1]);
          x2 = STAN_edit_x_cont(x2, teff_idx2[1], teff_zero[1], teff[1]);
        }
        
        // 3.2 Normalize
        x1 = x1 / X_scale[idx_cont];
        x2 = x2 / X_scale[idx_cont];
      }
      
      // 4. Handle possible nonstationarity
      if(is_warped > 0){
        idx_wrp += 1;
        
        // 4.1 Variance masking
        if(is_warped==2){
          Kj = Kj .* STAN_kernel_varmask(x1, x2, wrp[idx_wrp], vm_params);
        }
        
        // 4.2 Input warping
        x1 = STAN_warp_input(x1, wrp[idx_wrp]);
        x2 = STAN_warp_input(x2, wrp[idx_wrp]);
      }
      
      // Compute the kernel matrix
      idx_alpha += 1;
      if(ker_cont != 0){
        idx_ell += 1;
        Kj = Kj .* STAN_kernel_eq(x1, x2, alpha[idx_alpha], ell[idx_ell]);
      } else {
        Kj = square(alpha[idx_alpha]) * Kj;
      }
      
      // Possible heterogeneity
      if(is_heter){
        Kj = Kj .* STAN_kernel_beta(beta[1], beta_idx1[1], beta_idx2[1]);
      }
      
      KX[j] = Kj; // store kernel matrix
    }
    
    return(KX);
  }
