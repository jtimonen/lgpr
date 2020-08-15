// Precompute fixed kernel matrices
matrix[num_obs, num_obs] K_const[num_comps] = 
  STAN_kernel_const_all(num_obs, num_obs,
    x_cat, x_cat, x_cont_mask, x_cont_mask, x_cat_num_levels, components);

// Delta vector for diagonal jitter
vector[num_obs] delta_vec = rep_vector(delta, num_obs);

// Other options
int is_generated_done = (1 - is_f_sampled) * (1 - is_generated_skipped);

// Print info
if(is_verbose){
  print("");
  print("Stan info:");
  print(" - number of observations = ", num_obs);
  print(" - number of components = ", num_comps);
  print(" - number of continuous covariates = ", num_cov_cont);
  print(" - number of categorical covariates = ", num_cov_cat);
  print(" - observation model = ", obs_model);
  print("");
}

