// Count how many times <val> occurs in <arr>
int STAN_count(int[] arr, int val){
  int L = size(arr);
  int n = 0;
  for(i in 1:L){
    if(arr[i]==val){
      n += 1;
    }
  }
  return(n);
}

// Count required number of different parameters
int[] STAN_count_num_params(int[ , ] components, int obs_model){
  int num_comps = size(components[1]);
  int ct[num_comps] = components[1];
  
  int num_alpha = size(ct);
  int num_ell = STAN_count(ct, 1) + STAN_count(ct, 2) + STAN_count(ct, 3);
  int num_steepness = STAN_count(ct, 3);
  int num_sigma = (obs_model == 1);
  int num_phi = (obs_model == 3);
  
  int num_params[5] = {
    num_alpha, 
    num_ell,
    num_steepness,
    num_sigma,
    num_phi
  };
  return(num_params);
}

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
