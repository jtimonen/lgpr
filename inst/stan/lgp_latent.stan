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
#include _latent/data-general.stan
#include _latent/data-response.stan
}

transformed data{
#include _common/tdata.stan
}

parameters {
#include _common/params.stan
#include _latent/params.stan
  vector[N] eta[J]; // isotropic versions of func components
}

transformed parameters {
  vector[N] f_latent[J];
#include _common/tparams.stan
  {
    matrix[N, N] Delta = diag_matrix(delta_vec);
    matrix[N, N] KX[J] = STAN_kernel_all(
      N, N, K_const, components, X, X, X_scale,
      alpha, ell, wrp, beta, teff, vm_params, 
      BETA_IDX, BETA_IDX, TEFF_IDX, TEFF_IDX, teff_zero
    );
    for(j in 1:J){
      f_latent[j] = cholesky_decompose(KX[j] + Delta) * eta[j];
    }
  }
}

model {
  vector[N] f_sum = STAN_vectorsum(f_latent, N);
  if(obs_model>1) {f_sum = f_sum + c_hat[1];}
  for(j in 1:J){ target += normal_lpdf(eta[j] | 0, 1); }
#include _common/priors.stan
#include _latent/priors.stan
#include _latent/likelihood.stan
}
