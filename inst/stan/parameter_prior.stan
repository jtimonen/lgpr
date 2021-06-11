#include _common/licence.stan

functions{
#include _common/functions-prior.stan
}

data {
#include _common/data-general.stan
#include _common/data-covariates.stan
#include _common/data-priors.stan
  int<lower=1,upper=5> obs_model; // 1-5: Gaussian, Poisson, NB, Bin, BB
  int<lower=0> prior_sigma[obs_model==1, 2];
  int<lower=0> prior_phi[obs_model==3, 2];
  real hyper_sigma[obs_model==1, 3];
  real hyper_phi[obs_model==3, 3];
  real hyper_gamma[obs_model==5, 2];
}

parameters {
#include _common/params.stan
  real<lower=1e-12> sigma[obs_model==1];
  real<lower=1e-12> phi[obs_model==3];
  real<lower=1e-12, upper=1-1e-12> gamma[obs_model==5];
}

transformed parameters {
#include _common/tparams.stan
}

model {
  // Priors
#include _common/priors.stan
  if(obs_model==1){
    target += STAN_log_prior(sigma[1], prior_sigma[1], hyper_sigma[1]);
  }else if(obs_model==3){
    target += STAN_log_prior(phi[1], prior_phi[1], hyper_phi[1]);
  }else if(obs_model==5){
    target += beta_lpdf(gamma[1] | hyper_gamma[1][2], hyper_gamma[1][2]);
  }
}
