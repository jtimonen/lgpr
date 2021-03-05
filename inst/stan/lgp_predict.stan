#include _common/licence.stan

functions{
#include _common/functions-utils.stan
#include _common/functions-kernels.stan
#include _common/functions-posterior.stan
}

data {
#include _common/data-general.stan
#include _common/data-covariates.stan
#include _common/data-pred.stan
#include _common/data-draws.stan
  real<lower=0> d_sigma[S];
  vector[num_obs] y_norm;
}

transformed data{
#include _common/tdata.stan
#include _common/tdata-pred.stan
}

generated quantities {
  vector[num_pred] F_POST[S, 2*(num_comps+1)];
  for (is in 1:S) {
    vector[num_pred] f_post[2*(num_comps+1)];
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
    F_POST[is] = STAN_gp_posterior(KX, KX_s, KX_ss, y_norm, delta, d_sigma[is]);
  }
}
