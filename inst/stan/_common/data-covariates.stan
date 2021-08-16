  // Continuous covariates
  vector[N] X[num_X];
  int<lower=0, upper=1> X_mask[num_X, N]; // missing? (0 = no, 1 = yes)
  real X_scale[num_X];
  
  // Categorical covariates
  int<lower=1> Z[num_Z, N];
  int<lower=0> Z_M[num_Z]; // numbers of categories
  
  // Expanding beta or teff to a vector of length equal to N.
  // The value *_IDX[n]-1 tells the index of the beta or teff parameter
  // that should be used for observation n. If observation n doesn't
  //  correspond to any beta parameter, then *_IDX[n] should be 1.
  int<lower=1, upper=num_beta+1> BETA_IDX[num_het>0, N];
  int<lower=1, upper=num_teff+1> TEFF_IDX[idx_unc>0, N];
