// INPUT WARP
vector STAN_warp_input(vector x, real a){
  return( -1 + 2*inv(1+exp(-a*x)) );
}

// VARIANCE MASK
vector STAN_var_mask(vector x, real a){
  return( inv(1+exp(-a*x)) );
}

// COMPUTE X_TILDE (corrected diseaseAge)
vector STAN_get_x_tilde(vector x_disAge, vector T_effect, vector T_observed, int[,] mapping, int[] lengths){
  int n_tot = num_elements(x_disAge);
  int N_cases = num_elements(lengths);
  vector[n_tot] x_tilde = rep_vector(0.0, n_tot);
  for(k in 1:N_cases){
    int inds[lengths[k]] = mapping[k, 1:lengths[k]];
    x_tilde[inds] = x_disAge[inds] + T_observed[k] - T_effect[k];
  }
  return(x_tilde);
}
