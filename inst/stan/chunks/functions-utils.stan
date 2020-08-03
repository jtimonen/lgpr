// INPUT WARP
vector STAN_warp_input(vector x, real a){
  return( -1 + 2*inv(1+exp(-a*x)) );
}

// VARIANCE MASK FUNCION
vector STAN_var_mask(vector x, real a){
  return( inv(1+exp(-a*x)) );
}

// EDIT DISEASE AGE ACCORDING TO UNCERTAINTY
vector STAN_edit_dis_age(
  data vector x_dis_age,
  data int[] idx_expand,
  data vector teff_obs,
  data vector teff)
{
  int n_tot = num_elements(x_dis_age);
  //int N_cases = num_elements(lengths);
  //vector[n_tot] x_tilde = rep_vector(0.0, n_tot);
  //for(k in 1:N_cases){
  //  int inds[lengths[k]] = mapping[k, 1:lengths[k]];
  //  x_tilde[inds] = x_disAge[inds] + T_observed[k] - T_effect[k];
  //}
  return(x_dis_age);
}
