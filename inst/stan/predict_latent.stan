#include _common/licence.stan

functions{
#include _common/functions-kernels.stan
#include _common/functions-kr.stan
}

data {
#include _common/data-general.stan
#include _common/data-covariates.stan
#include _common/data-pred.stan
#include _common/data-draws.stan
  int<lower=1,upper=5> obs_model; // 1-5: Gaussian, Poisson, NB, Bin, BB
  vector<lower=0>[num_pred] d_f_latent[S, num_comps];
  real<lower=0> d_sigma[S, obs_model==1];
  real<lower=0> d_phi[S, obs_model==1];
  real<lower=0> d_gamma[S, obs_model==1];
}

transformed data{
#include _common/tdata.stan
#include _common/tdata-pred.stan
}

generated quantities {
  vector[num_pred] F_KR[S, num_comps];
  for (is in 1:S) {
    
    // Compute kernel matrices
    matrix[num_obs, num_obs] KX[num_comps] = STAN_kernel_all(
      num_obs, num_obs, K_const, components, 
      x_cont, x_cont, x_cont_unnorm, x_cont_unnorm,
      d_alpha[is], d_ell[is], d_wrp[is], d_beta[is], d_teff[is],
      vm_params, idx_expand, idx_expand, teff_zero
    );
    matrix[num_pred, num_obs] KX_s[num_comps] = STAN_kernel_all(
      num_pred, num_obs, K_const_s, components, 
      x_cont_PRED, x_cont, x_cont_unnorm_PRED, x_cont_unnorm,
      d_alpha[is], d_ell[is], d_wrp[is], d_beta[is], d_teff[is],
      vm_params, idx_expand_PRED, idx_expand, teff_zero
    );
    matrix[num_pred, num_pred] KX_ss[num_comps] = STAN_kernel_all(
      num_pred, num_pred, K_const_ss, components,
      x_cont_PRED, x_cont_PRED, x_cont_unnorm_PRED, x_cont_unnorm_PRED,
      d_alpha[is], d_ell[is], d_wrp[is], d_beta[is], d_teff[is],
      vm_params, idx_expand_PRED, idx_expand_PRED, teff_zero
    );
    
    // Compute the posterior means and variances
    F_KR[is] = STAN_kr(d_f_latent[is], KX, KX_s, KX_ss, delta);
  }
}
