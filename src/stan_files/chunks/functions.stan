
// The input warping function 
vector warp_input(vector x, real a, real b, real c){
  int n_tot = num_elements(x);
  vector[n_tot] w = 2*c*(-0.5 + rep_vector(1, n_tot)./(1+exp(-a*(x-b))));
  return(w);
}

// The variance mask function
vector var_mask(vector x, real a){
  int n_tot = num_elements(x);
  vector[n_tot] s = rep_vector(1, n_tot)./(1+exp(-a*x));
  return(s);
}

// Compute X_tilde
vector get_x_tilde(vector x_disAge, vector T_onset, vector T_observed, int[,] mapping, int[] lengths){
  int n_tot = num_elements(x_disAge);
  int N_cases = num_elements(lengths);
  vector[n_tot] x_tilde = rep_vector(0.0, n_tot);
  for(k in 1:N_cases){
    int inds[lengths[k]] = mapping[k, 1:lengths[k]];
    x_tilde[inds] = x_disAge[inds] + T_observed[k] - T_onset[k];
  }
  return(x_tilde);
}
