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
  vector[num_obs] y_norm;
  array[1, 2] int<lower=0> prior_sigma;
  array[1, 3] real hyper_sigma;
}

transformed data{
  vector[num_obs] m0 = rep_vector(0.0, num_obs);
#include _common/tdata.stan
}

parameters {
#include _common/params.stan
  array[1] real<lower=1e-12> sigma;
}

transformed parameters {
#include _common/tparams.stan
}

model {
  // Priors
#include _common/priors.stan
  target += STAN_log_prior(sigma[1], prior_sigma[1], hyper_sigma[1]);

  // Likelihood
  if (is_likelihood_skipped == 0) {
    matrix[num_obs, num_obs] Ky = diag_matrix(delta_vec);
    array[num_comps] matrix[num_obs, num_obs] KX = STAN_kernel_all(num_obs, num_obs,
      K_const, components, x_cont, x_cont, x_cont_unnorm, x_cont_unnorm,
      alpha, ell, wrp, beta, teff,
      vm_params, idx_expand, idx_expand, teff_zero
    );
    for(j in 1:num_comps){
      Ky += KX[j];
    }
    Ky = add_diag(Ky, square(sigma[1]));
    y_norm ~ multi_normal_cholesky(m0, cholesky_decompose(Ky));
  }
}
