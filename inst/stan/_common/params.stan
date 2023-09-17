  array[num_comps] real<lower=1e-12> alpha;
  array[num_ell] real<lower=1e-12> ell;
  array[num_ns] real<lower=1e-12> wrp;
  array[num_heter>0] vector<lower=1e-12, upper=1-1e-12>[num_bt] beta;
  array[num_uncrt>0] vector<lower=1e-12, upper=1-1e-12>[num_bt] teff_raw;
