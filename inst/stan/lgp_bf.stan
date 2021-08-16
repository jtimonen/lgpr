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
  int<lower=1,upper=5> obs_model; // 1-5: Gaussian, Poisson, NB, Bin, BB
  int<lower=0> y_int[obs_model>1, N]; // response variable (int)
  real y_real[obs_model==1, N]; // response variable (real)
  int<lower=1> y_num_trials[obs_model>3, N]; // for Bin or BB model
  int<lower=0> prior_sigma[obs_model==1, 2];
  int<lower=0> prior_phi[obs_model==3, 2];
  real hyper_sigma[obs_model==1, 3];
  real hyper_phi[obs_model==3, 3];
  real hyper_gamma[obs_model==5, 2];
  vector[N] c_hat;
  
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
  real<lower=1e-12> sigma[obs_model==1];
  real<lower=1e-12> phi[obs_model==3];
  real<lower=1e-12, upper=1-1e-12> gamma[obs_model==5];
  vector[num_xi] xi; // bf multipliers
}

transformed parameters {
  vector[N] f_latent[J];
#include _common/tparams.stan
  {
    //matrix[N, N] KX[J] = STAN_kernel_all(
    //  N, N, K_const, components, 
    //  x_cont, x_cont, x_cont_unnorm, x_cont_unnorm,
    //  alpha, ell, wrp, beta, teff, 
    //  vm_params, idx_expand, idx_expand, teff_zero
    //);
    //for(j in 1:J){
    //  f_latent[j] = 0.0;
    //}
  }
}

model {
  vector[N] f_sum = rep_vector(0.0, N); //STAN_vectorsum(f_latent, N) + c_hat;

  // Priors
#include _common/priors.stan
  xi ~ normal_lpdf(0, 1);
  if(obs_model==1){
    target += STAN_log_prior(sigma[1], prior_sigma[1], hyper_sigma[1]);
  }else if(obs_model==3){
    target += STAN_log_prior(phi[1], prior_phi[1], hyper_phi[1]);
  }else if(obs_model==5){
    target += beta_lpdf(gamma[1] | hyper_gamma[1][2], hyper_gamma[1][2]);
  }

  // Likelihood
  if(obs_model==1 && is_likelihood_skipped==0) {
    // 1. Gaussian
    real MU[N] = to_array_1d(f_sum); // means
    target += normal_lpdf(y_real[1] | MU, sigma[1]);
  }else if(obs_model==2 && is_likelihood_skipped==0){
    // 2. Poisson
    real LOG_MU[N] = to_array_1d(f_sum); // means (log-scale)
    target += poisson_log_lpmf(y_int[1] | LOG_MU);
  }else if(obs_model==3 && is_likelihood_skipped==0){
    // 3. Negative binomial
    real LOG_MU[N] = to_array_1d(f_sum); // means (log-scale)
    real PHI[N] = to_array_1d(rep_vector(phi[1], N)); // dispersion
    target += neg_binomial_2_log_lpmf(y_int[1] | LOG_MU, PHI);
  }else if(obs_model==4 && is_likelihood_skipped==0){
    // 4. Binomial
    real LOGIT_P[N] = to_array_1d(f_sum); // p success (logit-scale)
    target += binomial_logit_lpmf(y_int[1] | y_num_trials[1], LOGIT_P);
  }else if(obs_model==5 && is_likelihood_skipped==0){
    // 5. Beta-binomial
    real tgam = inv(gamma[1]) - 1.0;
    vector[N] P = inv_logit(f_sum); // p success
    real aa[N] = to_array_1d(P * tgam);
    real bb[N] = to_array_1d((1.0 - P) * tgam);
    target += beta_binomial_lpmf(y_int[1] | y_num_trials[1], aa, bb);
  }
}
