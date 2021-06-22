  // Compute total signal by summing components
  vector STAN_vectorsum(vector[] vecs, data int L){
    int num_vecs = size(vecs);
    vector[L] s = rep_vector(0, L);
    for (j in 1:num_vecs){
      s += vecs[j];
    }
    return(s);
  }
  
  // Sum an array of matrices
  matrix STAN_matrix_array_sum(matrix[] K){
    int n1 = rows(K[1]);
    int n2 = cols(K[1]);
    matrix[n1, n2] K_sum = K[1];
    for(j in 2:size(K)){
      K_sum += K[j];
    }
    return(K_sum);
  }
  
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
