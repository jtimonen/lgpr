  if(obs_model==1 && is_likelihood_skipped==0) {
    // 1. Gaussian
    target += normal_lpdf(y[1] | f_sum, sigma[1]);
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
