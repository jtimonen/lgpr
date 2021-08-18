  // Priors types and hyperparameters
  int<lower=0> prior_alpha[J, 2];  // {prior_type, transform}
  int<lower=0> prior_ell[num_ell, 2];      // {prior_type, transform}
  int<lower=0> prior_wrp[num_wrp, 2];       // {prior_type, transform}
  int<lower=0> prior_teff[num_unc>0, 2]; // {prior type, is_backwards}
  real hyper_alpha[J, 3];
  real hyper_ell[num_ell, 3];
  real hyper_wrp[num_wrp, 3];
  real hyper_teff[num_unc>0, 3];
  real hyper_beta[num_het>0, 2];
