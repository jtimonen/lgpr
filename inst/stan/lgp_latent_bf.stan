#include _common/licence.stan

functions{
#include _common/functions-utils.stan
#include _common/functions-prior.stan
#include _bf/functions-basisfun.stan
#include _bf/functions-kernels_decomp.stan
}

data {
#include _common/data-general.stan
#include _common/data-covariates.stan
#include _common/data-priors.stan
#include _latent/data-general.stan
#include _latent/data-response.stan

  // Categorical decomposition stuff
  int<lower=1> len_eigvals;
  int<lower=1> len_eigvecs;
  int<lower=0> C_sizes[J];
  vector<lower=0>[len_eigvals] C_eigvals;
  vector[len_eigvecs] C_eigvecs;
  
  // Basis function stuff
  real<lower=0> scale_bf;
  int<lower=1> num_bf;
  int<lower=0> num_xi;
}

transformed data{
}

parameters {
#include _common/params.stan
#include _latent/params.stan
  vector[num_xi] xi; // bf multipliers
}

transformed parameters {
  vector[N] f_latent[J];
#include _common/tparams.stan
}

model {
  vector[N] f_sum = rep_vector(0.0, N); //STAN_vectorsum(f_latent, N) + c_hat;
#include _common/priors.stan
#include _latent/priors.stan
#include _latent/likelihood.stan
}
