// Precompute fixed kernel matrices
matrix[num_obs, num_obs] K_fixed[num_comps] = 
  STAN_kernel_fixed_all(x_disc, x_disc, num_levels, components);

// Delta vector for diagonal jitter
vector[num_obs] delta_vec = rep_vector(delta, num_obs);

// Other options
int is_generated_done = (1 - is_f_sampled) * (1 - is_generated_skipped);

