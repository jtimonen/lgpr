// Build the additive function components
// To be used only when the function components are being sampled
matrix[n_tot,n_tot] Kxr;
matrix[n_tot,n_tot] Lxr;
matrix[n_tot,n_tot] EYE = diag_matrix(rep_vector(DELTA, n_tot));
int r = 0;
if(D[1]==1){
  r += 1;
  Kxr = cov_exp_quad(x_age, alpha_idAge[1], lengthscale_idAge[1]) .* KF[1];
  Lxr = cholesky_decompose(Kxr + EYE);
  F[1,r,] = Lxr*ETA[1,r,];
}
if(D[2]==1){
  r += 1;
  Kxr = cov_exp_quad(x_age, alpha_sharedAge[1], lengthscale_sharedAge[1]);
  Lxr = cholesky_decompose(Kxr + EYE);
  F[1,r,] = Lxr*ETA[1,r,];
}
if(D[3]==1){
  real alp = alpha_diseaseAge[1];
  real ell = lengthscale_diseaseAge[1];
  real stp = warp_steepness[1];
  vector[n_tot] x_tilde;
  real w[n_tot];
  r += 1;
    
  // Handle diseaseAge uncertainty
  if(UNCRT==0){
    x_tilde = X[3];
  }else{
    x_tilde = get_x_tilde(X[3], T_onset[1], T_observed, caseID_to_rows, caseID_nrows); 
  }
  w = to_array_1d(warp_input(x_tilde, stp, 0.0, 1.0));
    
  // Create disease effect kernel
  Kxr = KF[2] .* cov_exp_quad(w, alp, ell);
  if(HMGNS==0){
    Kxr = K_beta(beta[1], row_to_caseID) .* Kxr;
  }
  if(USE_VAR_MASK==1){
    Kxr = K_var_mask(x_tilde, stp) .* Kxr;
  }
  Lxr = cholesky_decompose(Kxr + EYE);
  F[1,r,] = Lxr*ETA[1,r,];
}
for(j in 1:D[4]){
  r += 1;
  Kxr = cov_exp_quad(to_array_1d(X[2+D[3]+j]), alpha_continuous[j], lengthscale_continuous[j]);
  Lxr = cholesky_decompose(Kxr + EYE);
  F[1,r,] = Lxr*ETA[1,r,];
}
for(j in 1:D[5]){
  r += 1;
  Kxr = cov_exp_quad(x_age, alpha_categAge[j], lengthscale_categAge[j]) .* KF[1+D[3]+j];
  Lxr = cholesky_decompose(Kxr + EYE);
  F[1,r,] = Lxr*ETA[1,r,];
}
for(j in 1:D[6]){
  r += 1;
  Kxr = square(alpha_categOffset[j]) * KF[1+D[3]+D[5]+j];
  Lxr = cholesky_decompose(Kxr + EYE);
  F[1,r,] = Lxr*ETA[1,r,];
}

