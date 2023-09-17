  // Covariates
  array[num_cov_cont] vector[num_obs] x_cont;
  array[num_cov_cont] vector[num_obs] x_cont_unnorm;
  array[num_cov_cont, num_obs] int x_cont_mask;
  array[num_cov_cat, num_obs] int x_cat;
  
  // Expanding beta or teff to a vector of length equal to num_obs.
  // The value idx_expand[j]-1 tells the index of the beta or teff parameter
  // that should be used for observation j. If observation j doesn't
  //  correspond to any beta parameter, then idx_expand[j] should be 1.
  array[num_obs] int<lower=1, upper=num_bt+1> idx_expand;
