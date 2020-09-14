/* Component-wise and total posterior means and stds. Let D be the
  number of components (<num_comps>). Then the elemens of <f_post> are
    - [1:D] = component means
    - [D+1] = total mean
    - [(D+2):(2*D+1)] = component stds
    - [2*D + 2] = total std
*/
vector[num_obs] f_post[is_f_sampled==0, 2*(num_comps+1)];
int y_rng_disc[is_f_sampled && (obs_model != 1), num_obs];
real y_rng_cont[is_f_sampled && (obs_model == 1), num_obs];

if(is_f_sampled==0){
  
  // Compute all kernel matrices
  matrix[num_obs, num_obs] KX[num_comps] = STAN_kernel_all(num_obs, num_obs,
      K_const, components, x_cont, x_cont, x_cont_unnorm, x_cont_unnorm,
      alpha, ell, wrp, beta, teff,
      vm_params, idx_expand, idx_expand, teff_zero);
      
  // Compute component-wise and total function posteriors
  f_post[1] = STAN_gp_posterior(KX, KX, KX, y_cont[1], delta, sigma[1]);
  
} else {
  
  vector[num_obs] f_sum = STAN_vectorsum(f_latent[1], num_obs) + c_hat;
  if(obs_model==2){
    // 2. Poisson
    real LOG_MU[num_obs] = to_array_1d(f_sum); // means (log-scale)
    y_rng_disc[1] = poisson_log_rng(LOG_MU);
  }else if(obs_model==3){
    // 3. Negative binomial
    real LOG_MU[num_obs] = to_array_1d(f_sum); // means (log-scale)
    real PHI[num_obs] = to_array_1d(rep_vector(phi[1], num_obs)); // dispersion
    y_rng_disc[1] = neg_binomial_2_log_rng(LOG_MU, PHI);
  }else if(obs_model==4){
    // 4. Binomial
    real P[num_obs] = to_array_1d(inv_logit(f_sum)); // p success
    y_rng_disc[1] = binomial_rng(y_num_trials[1], P);
  }else if(obs_model==5){
    // 5. Beta-binomial
    real tgam = inv(gamma[1]) - 1.0;
    vector[num_obs] P = inv_logit(f_sum); // p success
    real aa[num_obs] = to_array_1d(P * tgam);
    real bb[num_obs] = to_array_1d((1.0 - P) * tgam);
    y_rng_disc[1] = beta_binomial_rng(y_num_trials[1], aa, bb);
  }else{
    // 1. Gaussian (obs_model should be 1)
    real MU[num_obs] = to_array_1d(f_sum); // means
    real SIGMA[num_obs] = to_array_1d(rep_vector(sigma[1], num_obs)); // stds
    y_rng_cont[1] = normal_rng(MU, SIGMA);
  }
}
