
if(F_IS_SAMPLED){
  
  // Compute total signal + offset c_hat
  vector[n] F_ss = C_hat;
  for (i in 1:n){
    F_ss[i] += sum(F[1,,i]);
  }

  // Compute likelihood
  if(LH==1){
    // 1. Gaussian observation model
    real SIGMA[n] = to_array_1d(rep_vector(sigma_n[1], n)); // means
    real MU[n] = to_array_1d(F_ss);                   // stds
    target += normal_lpdf(y | MU, SIGMA);
  }else if(LH==2){
    // 2. Poisson observation model
    real LOG_MU[n] = to_array_1d(F_ss); // means (log-scale)
    target += poisson_log_lpmf(y_int | LOG_MU);
  }else if(LH==3){
    // 3. Negative binomial observation model
    real LOG_MU[n] = to_array_1d(F_ss); // means (log-scale)
    real PHI[n] = to_array_1d(rep_vector(phi[1], n)); // dispersion param
    target += neg_binomial_2_log_lpmf(y_int | LOG_MU, PHI);
  }else if(LH==4){
    // 4. Bernoulli or binomial observation model
    real LOGIT_P[n] = to_array_1d(F_ss); // p success (logit-scale)
    target += binomial_logit_lpmf(y_int | N_trials, LOGIT_P);
  }else{
    reject("Unknown observation model!")
  }

}else{
  // F NOT SAMPLED
  matrix[n,n] Ky;
  matrix[n,n] Ly;
  matrix[n,n] Kx = diag_matrix(rep_vector(DELTA, n));
  matrix[n,n] KX[sum_D] = STAN_compute_kernel_matrices(X, caseID_to_rows, row_to_caseID, caseID_nrows, KF, T_effect, T_observed, D, UNCRT, HMGNS, USE_VAR_MASK, vm_params, alpha_idAge, alpha_sharedAge,  alpha_diseaseAge, alpha_continuous, alpha_categAge, alpha_categOffset, ell_idAge, ell_sharedAge, ell_diseaseAge, ell_continuous, ell_categAge, warp_steepness, beta);
  if(LH!=1){
    reject("Likelihood must be Gaussian if F is not sampled!")
  }
  for(j in 1:sum_D){
    Kx += KX[j];
  }
  Ky = Kx + diag_matrix(rep_vector(square(sigma_n[1]), n));
  Ly = cholesky_decompose(Ky);
  y ~ multi_normal_cholesky(C_hat, Ly);
}
