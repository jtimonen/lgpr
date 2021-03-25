  real<lower=1e-12> alpha[num_comps];
  real<lower=1e-12> ell[num_ell];
  real<lower=1e-12> wrp[num_ns];
  vector<lower=1e-12, upper=1-1e-12>[num_bt] beta[num_heter>0];
  vector<lower=1e-12, upper=1-1e-12>[num_bt] teff_raw[num_uncrt>0];
