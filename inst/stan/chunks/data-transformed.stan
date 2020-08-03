// Precompute fixed kernel matrices
matrix[num_obs, num_obs] K_fixed[num_comps] = 
  STAN_kernel_fixed_all(x_disc, x_disc, num_levels, components);

// Other options
int is_gen_quant_done = (1 - is_f_sampled) * (1 - is_gen_quant_skipped);
