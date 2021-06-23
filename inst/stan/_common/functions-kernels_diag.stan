// These are NOT used during inference of Stan model, only exported to R for
// prediction computations.

  // Compute one constant kernel matrix diagonal.
  vector STAN_kernel_const_diag(data int[] x, data int kernel_type) 
  {
    int n = num_elements(x);
    vector[n] K_diag = rep_vector(1.0, n);
    if (kernel_type == 2) {
      vector[n] is_zero;
      for(j in 1:n) {
        is_zero[j] = (x[j]==0);
      }
      K_diag = is_zero; // binary mask: one if both inputs are 0
    }
    return(K_diag);
  }
  
  // Compute all constant kernel matrices' diagonals.
  vector[] STAN_kernel_const_all_diag(
    data int n,
    data int[,] x,
    data int[,] x_mask,  
    data int[,] components)
  {
    int num_comps = size(components);
    vector[n] K_const_diag[num_comps];
    for (j in 1:num_comps) {
      vector[n] K_diag;
      int opts[9] = components[j];
      int ctype = opts[1];
      int ktype = opts[2];
      int idx_cat = opts[8];
      int idx_cont = opts[9];
      
      // Compute mask kernel for continuous covariate
      if (idx_cont != 0) {
        K_diag = STAN_kernel_const_diag(x_mask[idx_cont], 2); // 2=binary mask
      } else {
        K_diag = rep_vector(1.0, n);
      }
      
      // Compute kernel for categorical covariate
      if (ctype == 0 || ctype == 2) {
        K_diag = K_diag .* STAN_kernel_const_diag(x[idx_cat], ktype);
      }
      K_const_diag[j] = K_diag;
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
    data vector[] x,
    data vector[] x_unnorm,
    real[] alpha,
    real[] wrp,
    vector[] beta,
    vector[] teff,
    data real[] vm_params,
    data int[] idx_expand,
    data vector[] teff_zero)
  {
    int idx_wrp = 0;
    int idx_alpha = 0;
    int num_comps = size(components);
    vector[n] KX_diag[num_comps];
  
    // Loop through components
    for(j in 1:num_comps){
      
      // 1. Initialize with constant part of the kernel diagonal
      vector[n] K_diag = K_const_diag[j];
      vector[n] X;
  
      // 2. Get component properties
      int opts[9] = components[j];
      int ctype = opts[1];
      int idx_cont = opts[9];
      int is_heter = opts[4];
      int is_warped = opts[5];
      int is_var_masked = opts[6];
      int is_uncrt = opts[7];
      
      // 3. Pick the possible continuous covariate of this component
      if(ctype != 0){
        if(is_warped){
          X = x_unnorm[idx_cont];
        }else{
          X = x[idx_cont];
        }
      }
      
      // 4. Handle possible nonstationarity
      if(is_warped){
        real s;
        idx_wrp += 1;
        
        // 4.1 Handle possible uncertainty
        if(is_uncrt){
          X = STAN_edit_x_cont(X, idx_expand, teff_zero[1], teff[1]);
        }
        
        // 4.2 Variance masking
        s = wrp[idx_wrp];
        if(is_var_masked){
          K_diag = K_diag .* STAN_kernel_varmask_diag(X, s, vm_params);
        }
        
        // 4.3 Input warping
        X = STAN_warp_input(X, s);
      }
      
      // Compute the kernel matrix  diagonal
      idx_alpha += 1;
      if(ctype != 0){
        K_diag = K_diag .* STAN_kernel_eq_diag(n, alpha[idx_alpha]);
      } else {
        K_diag = square(alpha[idx_alpha]) * K_diag;
      }
      
      // Possible heterogeneity
      if(is_heter){
        K_diag = K_diag .* STAN_kernel_beta_diag(beta[1], idx_expand);
      }
      
      KX_diag[j] = K_diag; // store kernel matrix diagonal
    }
    
    return(KX_diag);
  }
