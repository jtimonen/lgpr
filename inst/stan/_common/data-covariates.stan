  // Covariates
  vector[num_obs] x_cont[num_cov_cont];
  vector[num_obs] x_cont_unnorm[num_cov_cont];
  int x_cont_mask[num_cov_cont, num_obs];
  int x_cat[num_cov_cat, num_obs];
  int<lower=1, upper=num_bt+1> idx_expand[num_obs]; // expands beta and t_eff.
