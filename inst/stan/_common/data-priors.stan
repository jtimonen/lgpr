  // Priors types and hyperparameters
  int<lower=0> prior_alpha[num_comps, 2];  // {prior_type, transform}
  int<lower=0> prior_ell[num_ell, 2];      // {prior_type, transform}
  int<lower=0> prior_wrp[num_ns, 2];       // {prior_type, transform}
  int<lower=0> prior_teff[num_uncrt>0, 2]; // {prior type, is_backwards}
  real hyper_alpha[num_comps, 3];
  real hyper_ell[num_ell, 3];
  real hyper_wrp[num_ns, 3];
  real hyper_teff[num_uncrt>0, 3];
  real hyper_beta[num_heter>0, 2];
