// INPUT WARP
vector STAN_warp_input(vector x, real a){
  return( -1 + 2*inv(1+exp(-a*x)) );
}

// VARIANCE MASK FUNCION
vector STAN_var_mask(vector x, real a){
  return( inv(1+exp(-a*x)) );
}

// INPUT WARP
vector STAN_expand(vector v, data int[] idx_expand){
  int L = num_elements(v);
  vector[L+1] v_add0 = rep_vector(0.0, L+1);
  v_add0[1:(L+1)] = v;
  return(v_add0[idx_expand]);
}

// EDIT DISEASE AGE ACCORDING TO UNCERTAINTY
vector STAN_edit_dis_age(
  data vector x_dis_age,
  data int[] idx_expand,
  data vector teff_obs,
  vector teff)
{
  int n = num_elements(x_dis_age);
  vector[n] x_teff_obs = STAN_expand(teff_obs, idx_expand);
  vector[n] x_teff = STAN_expand(teff, idx_expand);
  vector[n] x_tilde = x_dis_age + x_teff_obs - x_teff;
  return(x_dis_age);
}
