// Kernel hyperparameters and observation model
real<lower=0> alpha[num_comps];
real<lower=0> ell[num_ell];
real<lower=0> wrp[num_ns];
real<lower=0> sigma[obs_model==1];
real<lower=0> phi[obs_model==3];

// Effect time uncertainty
vector<lower=0, upper=1>[num_bt] teff_raw[num_uncrt>0];

// Heterogeneous effect params
vector<lower=0, upper=1>[num_bt] beta[num_heter>0];

// Isotropic versions of function components
vector[num_obs] eta[is_f_sampled, num_comps]; 
