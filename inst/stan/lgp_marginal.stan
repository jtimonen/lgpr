#include _common/licence.stan

functions{
#include _common/functions-utils.stan
#include _common/functions-prior.stan
#include _common/functions-kernels.stan
}

data {
#include _common/data.stan
  vector[num_obs] y;
  int<lower=0> prior_sigma[2];
  real hyper_sigma[3];
}

transformed data{
#include _common/tdata.stan
}

parameters {
#include _common/params.stan
  real<lower=0> sigma;
}

transformed parameters {
#include _common/tparams.stan
}

model {
  // Priors
#include _common/priors.stan
  target += STAN_log_prior(sigma, prior_sigma, hyper_sigma);

  // Likelihood
  {
    vector[num_obs] c_hat = rep_vector(0.0, num_obs);
    vector[num_obs] sigma2_vec = rep_vector(square(sigma), num_obs);
    matrix[num_obs, num_obs] Ky = diag_matrix(delta_vec);
    matrix[num_obs, num_obs] KX[num_comps] = STAN_kernel_all(num_obs, num_obs,
      K_const, components, x_cont, x_cont, x_cont_unnorm, x_cont_unnorm,
      alpha, ell, wrp, beta, teff,
      vm_params, idx_expand, idx_expand, teff_zero
    );
    for(j in 1:num_comps){
      Ky += KX[j];
    }
    Ky = Ky + diag_matrix(sigma2_vec);
    y ~ multi_normal_cholesky(c_hat, cholesky_decompose(Ky));
  }
}
