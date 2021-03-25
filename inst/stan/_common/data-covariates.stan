  // Covariates
  vector[num_obs] x_cont[num_cov_cont];
  vector[num_obs] x_cont_unnorm[num_cov_cont];
  int x_cont_mask[num_cov_cont, num_obs];
  int x_cat[num_cov_cat, num_obs];
  
  // Expanding beta or teff to a vector of length equal to num_obs.
  // The value idx_expand[j]-1 tells the index of the beta or teff parameter
  // that should be used for observation j. If observation j doesn't
  //  correspond to any beta parameter, then idx_expand[j] should be 1.
  int<lower=1, upper=num_bt+1> idx_expand[num_obs];
