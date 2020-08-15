if(is_f_sampled){
  
  // Compute total signal + offset c_hat
  vector[num_obs] f_sum = c_hat;
  for (j in 1:num_comps){
    f_sum += f_latent[1,j];
  }

  // Compute likelihood
  if(obs_model==1){
    // 1. Gaussian observation model
    real MU[num_obs] = to_array_1d(f_sum); // means
    real SIGMA[num_obs] = to_array_1d(rep_vector(sigma[1], num_obs)); // stds
    target += normal_lpdf(y_cont[1] | MU, SIGMA);
  }else if(obs_model==2){
    // 2. Poisson observation model
    real LOG_MU[num_obs] = to_array_1d(f_sum); // means (log-scale)
    target += poisson_log_lpmf(y_disc[1] | LOG_MU);
  }else if(obs_model==3){
    // 3. Negative binomial observation model
    real LOG_MU[num_obs] = to_array_1d(f_sum); // means (log-scale)
    real PHI[num_obs] = to_array_1d(rep_vector(phi[1], num_obs)); // dispersion
    target += neg_binomial_2_log_lpmf(y_disc[1] | LOG_MU, PHI);
  }else if(obs_model==4){
    // 4. Bernoulli or binomial observation model
    real LOGIT_P[num_obs] = to_array_1d(f_sum); // p success (logit-scale)
    target += binomial_logit_lpmf(y_disc[1] | y_num_trials[1], LOGIT_P);
  }else{
    reject("<obs_model> must be 1, 2, 3 or 4!");
  }

}else{
  // F NOT SAMPLED
  vector[num_obs] sigma2_vec = rep_vector(square(sigma[1]), num_obs);
  matrix[num_obs, num_obs] Ky = diag_matrix(num_comps * delta_vec);
  matrix[num_obs, num_obs] KX[num_comps] = STAN_kernel_all(num_obs, num_obs,
      K_const, components, x_cont, x_cont, alpha, ell, wrp, beta, teff,
      vm_params, idx_expand, idx_expand, teff_obs);

  if(obs_model!=1){
    reject("<obs_model> must be 1 if the latent functions are not sampled!")
  }
  for(j in 1:num_comps){
    Ky += KX[j];
  }
  Ky = Ky + diag_matrix(sigma2_vec);
  target += multi_normal_lpdf(y_cont[1] | c_hat, Ky);
}
