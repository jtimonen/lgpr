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
  vector[N] y[1];
  int<lower=0> prior_sigma[1, 2];
  real hyper_sigma[1, 3];
}

transformed data{
  vector[N] m0 = rep_vector(0.0, N);
#include _common/tdata.stan
}

parameters {
#include _common/params.stan
  real<lower=1e-12> sigma[1];
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
    matrix[N, N] Ky = diag_matrix(delta_vec);
    matrix[N, N] KX[J] = STAN_kernel_all(N, N,
      K_const, components, X, X, X_scale,
      alpha, ell, wrp, beta, xpar, vm_params, 
      BETA_IDX, BETA_IDX, XPAR_IDX, XPAR_IDX, xpar_zero
    );
    for(j in 1:J){
      Ky += KX[j];
    }
    Ky = add_diag(Ky, square(sigma[1]));
    y[1] ~ multi_normal_cholesky(m0, cholesky_decompose(Ky));
  }
}
