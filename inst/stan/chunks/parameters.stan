// Kernel hyperparameters and observation model
real<lower=0> alpha[num_comps];
real<lower=0> ell[num_ell];
real<lower=0> wrp[num_dis];
real<lower=0> sigma[obs_model==1];
real<lower=0> phi[obs_model==3];

// Effect time uncertainty
vector<lower=0>[num_cases] teff_raw[is_uncrt];

// Heterogeneous effect params
vector<lower=0>[num_cases] beta[is_heter];

// Isotropic versions of function components
vector[num_obs] eta[is_f_sampled, num_comps]; 
