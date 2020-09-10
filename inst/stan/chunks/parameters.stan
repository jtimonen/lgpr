// Positive parameters
real<lower=0> alpha[num_comps];
real<lower=0> ell[num_ell];
real<lower=0> wrp[num_ns];
real<lower=0> sigma[obs_model==1];
real<lower=0> phi[obs_model==3];

// Parameters on [0, 1]
real<lower=0, upper=1> gamma[obs_model==5];
vector<lower=0, upper=1>[num_bt] beta[num_heter>0];
vector<lower=0, upper=1>[num_bt] teff_raw[num_uncrt>0];

// Isotropic versions of function components
vector[num_obs] eta[is_f_sampled, num_comps]; 
