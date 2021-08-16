  real<lower=1e-12> alpha[J];
  real<lower=1e-12> ell[num_ell];
  real<lower=1e-12> wrp[num_wrp];
  vector<lower=1e-12, upper=1-1e-12>[num_beta] beta[num_het>0];
  vector<lower=1e-12, upper=1-1e-12>[num_teff] teff_raw[idx_unc>0];
