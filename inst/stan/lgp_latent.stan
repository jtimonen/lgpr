#include _common/licence.stan

functions{
#include _common/functions-utils.stan
#include _common/functions-prior.stan
#include _common/functions-kernels.stan
}

data {
#include _common/data-general.stan
#include _common/data-covariates.stan
#include _common/data-priors.stan
  int<lower=1,upper=5> obs_model; // 1-5: Gaussian, Poisson, NB, Bin, BB
  array[obs_model>1, num_obs] int<lower=0> y_int; // response variable (int)
  array[obs_model==1, num_obs] real y_real; // response variable (real)
  array[obs_model>3, num_obs] int<lower=1> y_num_trials; // for Bin or BB model
  array[obs_model==1, 2] int<lower=0> prior_sigma;
  array[obs_model==3, 2] int<lower=0> prior_phi;
  array[obs_model==1, 3] real hyper_sigma;
  array[obs_model==3, 3] real hyper_phi;
  array[obs_model==5, 2] real hyper_gamma;
  vector[num_obs] c_hat; // GP mean vector
}

transformed data{
#include _common/tdata.stan
}

parameters {
#include _common/params.stan
  array[obs_model==1] real<lower=1e-12> sigma;
  array[obs_model==3] real<lower=1e-12> phi;
  array[obs_model==5] real<lower=1e-12, upper=1-1e-12> gamma;
  array[num_comps] vector[num_obs] eta; // isotropic versions of func components
}

transformed parameters {
  array[num_comps] vector[num_obs] f_latent;
#include _common/tparams.stan
  {
    matrix[num_obs, num_obs] Delta = diag_matrix(delta_vec);
    array[num_comps] matrix[num_obs, num_obs] KX = STAN_kernel_all(
      num_obs, num_obs, K_const, components, 
      x_cont, x_cont, x_cont_unnorm, x_cont_unnorm,
      alpha, ell, wrp, beta, teff, 
      vm_params, idx_expand, idx_expand, teff_zero
    );
    for(j in 1:num_comps){
      f_latent[j] = cholesky_decompose(KX[j] + Delta) * eta[j];
    }
  }
}

model {
  vector[num_obs] f_sum = STAN_vectorsum(f_latent, num_obs) + c_hat;

  // Priors
#include _common/priors.stan
  for(j in 1:num_comps){ target += normal_lpdf(eta[j] | 0, 1); }
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
    array[num_obs] real MU = to_array_1d(f_sum); // means
    target += normal_lpdf(y_real[1] | MU, sigma[1]);
  }else if(obs_model==2 && is_likelihood_skipped==0){
    // 2. Poisson
    array[num_obs] real LOG_MU = to_array_1d(f_sum); // means (log-scale)
    target += poisson_log_lpmf(y_int[1] | LOG_MU);
  }else if(obs_model==3 && is_likelihood_skipped==0){
    // 3. Negative binomial
    array[num_obs] real LOG_MU = to_array_1d(f_sum); // means (log-scale)
    array[num_obs] real PHI = to_array_1d(rep_vector(phi[1], num_obs)); // dispersion
    target += neg_binomial_2_log_lpmf(y_int[1] | LOG_MU, PHI);
  }else if(obs_model==4 && is_likelihood_skipped==0){
    // 4. Binomial
    array[num_obs] real LOGIT_P = to_array_1d(f_sum); // p success (logit-scale)
    target += binomial_logit_lpmf(y_int[1] | y_num_trials[1], LOGIT_P);
  }else if(obs_model==5 && is_likelihood_skipped==0){
    // 5. Beta-binomial
    real tgam = inv(gamma[1]) - 1.0;
    vector[num_obs] P = inv_logit(f_sum); // p success
    array[num_obs] real aa = to_array_1d(P * tgam);
    array[num_obs] real bb = to_array_1d((1.0 - P) * tgam);
    target += beta_binomial_lpmf(y_int[1] | y_num_trials[1], aa, bb);
  }
}
