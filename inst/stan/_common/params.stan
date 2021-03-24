  real<lower=0> alpha[num_comps];
  real<lower=0> ell[num_ell];
  real<lower=0> wrp[num_ns];
  vector<lower=0, upper=1>[num_bt] beta[num_heter>0];
  vector<lower=0, upper=1>[num_bt] teff_raw[num_uncrt>0];
