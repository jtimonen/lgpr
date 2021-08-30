  // Priors types and hyperparameters
  int<lower=0> prior_alpha[J, 2];  // {prior_type, transform}
  int<lower=0> prior_ell[num_ell, 2];      // {prior_type, transform}
  int<lower=0> prior_wrp[num_wrp, 2];       // {prior_type, transform}
  int<lower=0> prior_xpar[num_unc>0, 2]; // {prior type, is_backwards}
  real hyper_alpha[J, 3];
  real hyper_ell[num_ell, 3];
  real hyper_wrp[num_wrp, 3];
  real hyper_xpar[num_unc>0, 3];
  real hyper_beta[num_het>0, 2];
  
  // Observed effect times and uncertainty bounds
  vector[num_xpar] xpar_zero[num_unc>0];
  vector[num_xpar] xpar_lb[num_unc>0];
  vector[num_xpar] xpar_ub[num_unc>0];
