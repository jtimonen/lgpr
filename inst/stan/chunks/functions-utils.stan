// Check that number is a positive real number
void STAN_check_real_positive(real a){
  if(a<=0){ reject("argument must be positive!"); }
  if(is_nan(a)){ reject("argument must not be NaN!"); }
  if(is_inf(a)){ reject("argument must be finite!"); }
}

// Check that number is a non-zero probability
void STAN_check_prob_positive(real a){
  if(a<=0 || a > 1){ reject("argument must be on the interval (0, 1]!"); }
  if(is_nan(a)){ reject("argument must not be NaN!"); }
}

// Input warping function
vector STAN_warp_input(vector x, real a){
  STAN_check_real_positive(a);
  return( -1 + 2*inv(1+exp(-a*x)) );
}

// Variance masking function
vector STAN_var_mask(vector x, real a){
  STAN_check_real_positive(a);
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
