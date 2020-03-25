{
  int ix;
  int ikf;
  int r = 0;
  if(D[1]==1){
    r += 1;
    KX[r] = cov_exp_quad(x_age, alpha_idAge[1], ell_idAge[1]) .* KF[1];
  }
  if(D[2]==1){
    r += 1;
    KX[r] = cov_exp_quad(x_age, alpha_sharedAge[1], ell_sharedAge[1]);
  }
  if(D[3]==1){
    r += 1;
    KX[r] = STAN_K_disease(X[3], KF[2], alpha_diseaseAge, ell_diseaseAge, warp_steepness, vm_params, beta, USE_VAR_MASK, UNCRT, HMGNS, T_effect, T_observed, caseID_to_rows, row_to_caseID_plus1, caseID_nrows);
  }
  for(j in 1:D[4]){
    r += 1;
    ix = 2 + D[3] + j;
    KX[r] = cov_exp_quad(to_array_1d(X[ix]), alpha_continuous[j], ell_continuous[j]);
  }
  for(j in 1:D[5]){
    r += 1;
    ikf = 1 + D[3] + j;
    KX[r] = cov_exp_quad(x_age, alpha_categAge[j], ell_categAge[j]) .* KF[ikf];
  }
  for(j in 1:D[6]){
    r += 1;
    ikf = 1 + D[3] + D[5] + j;
    KX[r] = square(alpha_categOffset[j]) * KF[ikf];
  }
}
