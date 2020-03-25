// Compute likelihood
if(LH==1){
  // 1. Gaussian observation model
  real SIGMA[n] = to_array_1d(rep_vector(sigma_n[1], n)); // means
  real MU[n] = to_array_1d(F_ss);                         // stds
  target += normal_lpdf(y | MU, SIGMA);
}else if(LH==2){
  // 2. Poisson observation model
  real LOG_MU[n] = to_array_1d(F_ss); // means (log-scale)
  target += poisson_log_lpmf(y_int | LOG_MU);
}else if(LH==3){
  // 3. Negative binomial observation model
  real LOG_MU[n] = to_array_1d(F_ss);               // means (log-scale)
  real PHI[n] = to_array_1d(rep_vector(phi[1], n)); // dispersion param
  target += neg_binomial_2_log_lpmf(y_int | LOG_MU, PHI);
}else if(LH==4){
  // 4. Bernoulli or binomial observation model
  real LOGIT_P[n] = to_array_1d(F_ss); // p success (logit-scale)
  target += binomial_logit_lpmf(y_int | N_trials, LOGIT_P);
}else{
  reject("Unknown observation model!")
}
