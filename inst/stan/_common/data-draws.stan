  // Parameter draws
  int<lower=1> S;
  real<lower=0> d_alpha[S, num_comps];
  real<lower=0> d_ell[S, num_ell];
  real<lower=0> d_wrp[S, num_ns]; 
  vector<lower=0, upper=1>[num_bt] d_beta[S, num_heter>0];
  vector[num_bt] d_teff[S, num_uncrt>0]; // transformed effect time draws
