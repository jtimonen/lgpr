// Precompute fixed kernel matrices
matrix[num_obs, num_obs] K_const[num_comps] = 
  STAN_kernel_const_all(num_obs, num_obs,
    x_cat, x_cat, x_cont_mask, x_cont_mask, x_cat_num_levels, components);

// Delta vector for diagonal jitter
vector[num_obs] delta_vec = rep_vector(delta, num_obs);

// Constant kernel matrix ranks and eigenvalues
int ranks[num_comps] = STAN_ranks(components, x_cat_num_levels);
int R = sum(ranks);
int RM = R * num_basisfun;
vector[R] bfa_delta = STAN_delta_matrix(K_const, ranks, components);
matrix[num_obs, R] bfa_theta = STAN_theta_matrix(K_const, ranks, components);

