  // Input warping function
  vector STAN_warp_input(vector x, real a){
    return( -1 + 2*inv(1+exp(-a*x)) );
  }
  
  // Variance masking function
  vector STAN_var_mask(vector x, real a){
    return( inv(1+exp(-a*x)) );
  }
  
  // Expand a vector
  vector STAN_expand(vector v, data int[] idx_expand){
    int L = num_elements(v);
    vector[L+1] v_add0 = rep_vector(0.0, L+1);
    v_add0[2:(L+1)] = v;
    return(v_add0[idx_expand]);
  }
  
  // Edit a continuous covariate according to sampled uncertainty
  vector STAN_edit_x_cont(
    vector x_cont,
    data int[] idx_expand,
    data vector teff_obs,
    vector teff)
  {
    int n = num_elements(x_cont);
    vector[n] x_teff_obs = STAN_expand(teff_obs, idx_expand);
    vector[n] x_teff = STAN_expand(teff, idx_expand);
    return(x_cont + x_teff_obs - x_teff);
  }
  
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
  
  // BINARY MASK KERNEL
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
  matrix STAN_kernel_const(data int[] x1, data int[] x2, 
    data int kernel_type,  data int ncat) 
  {
    int n1 = num_elements(x1);
    int n2 = num_elements(x2);
    matrix[n1, n2] K;
    if (kernel_type == 1) {
      K = STAN_kernel_cat(x1, x2);
    } else if (kernel_type == 2) {
      K = STAN_kernel_bin(x1, x2);
    } else {
      // kernel_type should be 0
      K = STAN_kernel_zerosum(x1, x2, ncat);
    }
    return(K);
  }
  
  // Compute all constant kernel matrices. These do not depend on parameters and
  // therefore this function never needs to be evaluated during sampling.
  matrix[] STAN_kernel_const_all(
    data int n1,           data int n2,
    data int[,] x1,        data int[,] x2,
    data int[,] x1_mask,   data int[,] x2_mask,
    data int[] num_levels, data int[,] components)
  {
    int num_comps = size(components);
    matrix[n1, n2] K_const[num_comps];
    for (j in 1:num_comps) {
      matrix[n1, n2] K;
      int opts[9] = components[j];
      int ctype = opts[1];
      int ktype = opts[2];
      int idx_cat = opts[8];
      int idx_cont = opts[9];
      
      // Compute mask kernel for continuous covariate
      if (idx_cont != 0) {
        K = STAN_kernel_const(x1_mask[idx_cont], x2_mask[idx_cont], 2, 0);
      } else {
        K = rep_matrix(1, n1, n2);
      }
      
      // Compute kernel for categorical covariate
      if (ctype == 0 || ctype == 2) {
        int M = num_levels[idx_cat];
        K = K .* STAN_kernel_const(x1[idx_cat], x2[idx_cat], ktype, M);
      }
      K_const[j] = K;
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
    data vector[] x1,
    data vector[] x2,
    data vector[] x1_unnorm,
    data vector[] x2_unnorm,
    real[] alpha,
    real[] ell,
    real[] wrp,
    vector[] beta,
    vector[] teff,
    data real[] vm_params,
    data int[] idx1_expand,
    data int[] idx2_expand,
    data vector[] teff_zero)
  {
    int idx_ell = 0;
    int idx_wrp = 0;
    int idx_alpha = 0;
    int num_comps = size(components);
    matrix[n1, n2] KX[num_comps];
  
    // Loop through components
    for(j in 1:num_comps){
      
      // 1. Initialize with constant part of the kernel matrix
      matrix[n1, n2] K = K_const[j];
      vector[n1] X1;
      vector[n2] X2;
  
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
          X1 = x1_unnorm[idx_cont];
          X2 = x2_unnorm[idx_cont];
        }else{
          X1 = x1[idx_cont];
          X2 = x2[idx_cont];
        }
      }
      
      // 4. Handle possible nonstationarity
      if(is_warped){
        real s;
        idx_wrp += 1;
        
        // 4.1 Handle possible uncertainty
        if(is_uncrt){
          X1 = STAN_edit_x_cont(X1, idx1_expand, teff_zero[1], teff[1]);
          X2 = STAN_edit_x_cont(X2, idx2_expand, teff_zero[1], teff[1]);
        }
        
        // 4.2 Variance masking
        s = wrp[idx_wrp];
        if(is_var_masked){
          K = K .* STAN_kernel_varmask(X1, X2, s, vm_params);
        }
        
        // 4.3 Input warping
        X1 = STAN_warp_input(X1, s);
        X2 = STAN_warp_input(X2, s);
      }
      
      // Compute the kernel matrix
      idx_alpha += 1;
      if(ctype != 0){
        idx_ell += 1;
        K = K .* STAN_kernel_eq(X1, X2, alpha[idx_alpha], ell[idx_ell]);
      } else {
        K = square(alpha[idx_alpha]) * K;
      }
      
      // Possible heterogeneity
      if(is_heter){
        K = K .* STAN_kernel_beta(beta[1], idx1_expand, idx2_expand);
      }
      
      KX[j] = K; // store kernel matrix
    }
    
    return(KX);
  }
