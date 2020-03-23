
// COMPUTE FIXED KERNEL MATRICES (do not depend on parameters)
matrix[] STAN_kernels_fixed(vector[] X, int[] X_nn, int[] D, int N_tot, int[] N_cat)
{
  int n = num_elements(X[1]);
  int nf = 1 + D[3] + D[5] + D[6];
  matrix[n,n] KF[nf];
  KF[1] = STAN_K_zerosum(X[1], X[1], N_cat[1]);
  for(j in 1:D[3]){
    KF[1+j] = STAN_K_bin(X_nn, X_nn, 1);
  }
  for(j in 1:D[5]){
    int ix = 2 + D[3] + D[4] + j;
    KF[1+D[3]+j] = STAN_K_zerosum(X[ix], X[ix], N_cat[1+j]);
  }
  for(j in 1:D[6]){
    int ix = 2 + D[3] + D[4] + D[5] + j;
    KF[1+D[3]+D[5]+j] = STAN_K_zerosum(X[ix], X[ix], N_cat[1+D[5]+j]);
  }
  return(KF);
}

}

// COMPUTE ALL OTHER KERNEL MATRICES
matrix[] STAN_kernels(vector[] X, int[,]  matrix[] KF, int[] D, real[] alpha_idAge, real[] alpha_sharedAge, real[] alpha_continuous, real[] alpha_categAge, real[] alpha_categOffset, real[] ell_idAge, real[] ell_sharedAge, real[] ell_continuous, real[] ell_categAge){
    
  int n = num_elements(X[1]);
  real x_age[n] = to_array_1d(X[2]);   // age covariate as an array
  int sum_D = sum(D);
  matrix[n,n] KX[sum_D];
  int r = 0;
  if(D[1]==1){
    r += 1;
    KX[r] = cov_exp_quad(x_age, alpha_idAge[1], ell_idAge[1]) .* KF[1];
  }
  if(D[2]==1){
    r += 1;
    KX[r] = cov_exp_quad(x_age, alpha_sharedAge[1], ell_sharedAge[1]);
  }
  for(j in 1:D[4]){
    real alp = alpha_continuous[j];
    real ell = ell_continuous[j];
    int ix  = 2 + D[3] + j;
    r += 1;
    KX[r] = cov_exp_quad(to_array_1d(X[ix]), alp, ell);
  }
  for(j in 1:D[5]){
    real alp = alpha_categAge[j];
    real ell = ell_categAge[j];
    int ikf = 1 + D[3] + j;
    r += 1;
    KX[r] = cov_exp_quad(x_age, alp, ell) .* KF[ikf];
  }
  for(j in 1:D[6]){
    int ikf = 1 + D[3] + D[5] + j;
    real alp = alpha_categOffset[j];
    r += 1;
    KX[r] = square(alp) * KF[ikf];
  }
  return(KX);
}


// COMPUTE DISEASE KERNEL MATRIX
matrix STAN_disease_kernel(vector x, matrix Kf, real[] alpha_diseaseAge, real[] ell_diseaseAge, real[] warp_steepness, real[] vm_params, vector[] beta, int USE_VAR_MASK, int UNCRT, int HMGNS, vector[] T_effect, vector T_observed, caseID_to_rows, int[] row_to_caseID_plus1, int[] caseID_nrows){
  real alp = alpha_diseaseAge[1];
  real ell = ell_diseaseAge[1];
  real stp = warp_steepness[1];
  vector[n] x_tilde;
  real w[n];

  // Handle diseaseAge uncertainty
  if(UNCRT==0){
    x_tilde = x;
  }else{
    x_tilde = STAN_get_x_tilde(x, T_effect[1], T_observed, caseID_to_rows, caseID_nrows);
  }
  
  // Input warping
  w = to_array_1d(STAN_warp_input(x_tilde, stp));

  // Create disease effect kernel
  Kx = Kf .* cov_exp_quad(w, alp, ell);
  if(HMGNS==0){ Kx = STAN_K_beta(beta[1], row_to_caseID_plus1) .* Kx; }
  if(USE_VAR_MASK==1){ Kx = STAN_K_var_mask(x_tilde, stp, vm_params) .* Kx;}
  return(Kx)
}
  