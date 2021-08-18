  int<lower=1,upper=5> obs_model; // 1-5: Gaussian, Poisson, NB, Binom, BB
  int<lower=0> prior_sigma[obs_model==1, 2];
  int<lower=0> prior_phi[obs_model==3, 2];
  real hyper_sigma[obs_model==1, 3];
  real hyper_phi[obs_model==3, 3];
  real hyper_gamma[obs_model==5, 2];
