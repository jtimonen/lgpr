// These are NOT used during inference of Stan model, only exported to R for
// prediction computations.

  // Compute one constant kernel matrix diagonal.
  vector STAN_kernel_const_diag(data int[] z, data int kernel_type) 
  {
    int P = num_elements(z);
    vector[P] K_diag = rep_vector(1.0, P);
    if (kernel_type == 3) {
      vector[P] is_zero;
      for(j in 1:P) {
        is_zero[j] = (z[j]==0);
      }
      K_diag = is_zero; // binary mask: one if both inputs are 0
    }
    return(K_diag);
  }
  
  // Compute all constant kernel matrices' diagonals.
  vector[] STAN_kernel_const_all_diag(
    data int P,
    data int[,] Z,
    data int[,] X_mask,  
    data int[,] components)
  {
    int J = size(components);
    vector[P] K_const_diag[J];
    for (j in 1:J) {
      vector[P] Kj_diag;
      int ker_cat = components[j, 1];
      int idx_cat = components[j, 3];
      int idx_cont = components[j, 4];
      
      // Compute mask kernel for continuous covariate
      if (idx_cont != 0) {
        Kj_diag = STAN_kernel_const_diag(X_mask[idx_cont], 3); // 3=binary mask
      } else {
        Kj_diag = rep_vector(1.0, P);
      }
      
      // Compute kernel for categorical covariate
      if (ker_cat != 0) {
        Kj_diag = Kj_diag .* STAN_kernel_const_diag(Z[idx_cat], ker_cat);
      }
      K_const_diag[j] = Kj_diag;
    }
    return(K_const_diag);
  }
  
  // Exponentiated quadratic kernel (with vector inputs)
  vector STAN_kernel_eq_diag(int n, real alpha) {
    return(rep_vector(square(alpha), n));
  }
  
  // Multiplier matrix to enable variance masking
  vector STAN_kernel_varmask_diag(vector x, real steepness,
    data real[] vm_params)
  {
    real a = steepness * vm_params[2];
    real r = inv(a)*logit(vm_params[1]);
    int n = num_elements(x);
    vector[n] f_vm = STAN_var_mask(x - r, a);
    return(f_vm .* f_vm);
  }
  
  // Multiplier matrix to enable heterogeneous effects
  vector STAN_kernel_beta_diag(vector beta, int[] idx_expand) {
    return(STAN_expand(beta, idx_expand));
  }
  
  /* 
    Compute diagonals of all kernel matrices.
  */
  vector[] STAN_kernel_all_diag(
    data int n,
    data vector[] K_const_diag,
    data int[,] components,
    data vector[] X,
    data real[] X_scale,
    real[] alpha,
    real[] wrp,
    vector[] beta,
    vector[] teff,
    data real[] vm_params,
    data int[,] beta_idx,
    data int idx_unc,
    data int[,] teff_idx,
    data vector[] teff_zero)
  {
    int idx_wrp = 0;
    int idx_alpha = 0;
    int J = size(components);
    vector[n] KX_diag[J];
  
    // Loop through components
    for(j in 1:J){
      
      // 1. Initialize with constant part of the kernel diagonal
      vector[n] Dj = K_const_diag[j];
      vector[n] x;
  
      // 2. Get component properties
      int ker_cont = components[j, 2];
      int idx_cont = components[j, 4];
      int is_warped = components[j, 5];
      int is_heter = components[j, 6];
      
      // 3. Continuous covariate
      if(idx_cont != 0){
        x = X[idx_cont];
        if(idx_unc>0 && idx_unc==idx_cont) {
          x = STAN_edit_x_cont(x, teff_idx[1], teff_zero[1], teff[1]);
        }
        x = x / X_scale[idx_cont];
      }
      
      // 4. Possible nonstationarity
      if(is_warped > 0){
        idx_wrp += 1;
        if(is_warped==2){
          Dj = Dj .* STAN_kernel_varmask_diag(x, wrp[idx_wrp], vm_params);
        }
        x = STAN_warp_input(x, wrp[idx_wrp]);
      }
      
      // Compute the kernel matrix  diagonal
      idx_alpha += 1;
      if(ker_cont != 0){
        Dj = Dj .* STAN_kernel_eq_diag(n, alpha[idx_alpha]);
      } else {
        Dj = square(alpha[idx_alpha]) * Dj;
      }
      
      // Possible heterogeneity
      if(is_heter){
        Dj = Dj .* STAN_kernel_beta_diag(beta[1], beta_idx[1]);
      }
      
      KX_diag[j] = Dj; // store kernel matrix diagonal
    }
    
    return(KX_diag);
  }
