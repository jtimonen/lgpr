// COMPUTE FIXED KERNEL MATRICES (do not depend on parameters)
matrix[] STANFUNC_compute_fixed_kernel_matrices(vector[] X1, vector[] X2, int[] X1_nn, int[] X2_nn, int[] D, int cat_interact_kernel){
  
  int n1 = num_elements(X1[1]);
  int n2 = num_elements(X2[1]);
  int nf = 1 + D[3] + D[5] + D[6];
  matrix[n1,n2] KF[nf];
  KF[1] = STANFUNC_K_cat(X1[1], X2[1]);
  for(j in 1:D[3]){
    KF[1+j] = STANFUNC_K_bin_int(X1_nn, X2_nn, 1);
  }
  for(j in 1:D[5]){
    int ix = 2 + D[3] + D[4] + j;
    if(cat_interact_kernel == 1){
      KF[1+D[3]+j] = STANFUNC_K_cat(X1[ix], X2[ix]);
    }else{
      KF[1+D[3]+j] = STANFUNC_K_bin_real(X1[ix], X2[ix], 1);
    }
  }
  for(j in 1:D[6]){
    int ix = 2+D[3]+D[4]+D[5]+j;
    KF[1+D[3]+D[5]+j] = STANFUNC_K_cat(X1[ix], X2[ix]);
  }
  return(KF);
}


// COMPUTE ALL KERNEL MATRICES
matrix[] STANFUNC_compute_kernel_matrices(vector[] X1, vector[] X2, int[,] caseID_to_rows_1, int[,] caseID_to_rows_2, int[] row_to_caseID_1, int[] row_to_caseID_2, int[] caseID_nrows_1, int[] caseID_nrows_2, matrix[] KF, vector[] T_onset, vector T_observed, int[] D, int UNCRT, int HMGNS, int USE_VAR_MASK, real[] vm_params, real[] alpha_idAge, real[] alpha_sharedAge, real[] alpha_diseaseAge, real[] alpha_continuous, real[] alpha_categAge, real[] alpha_categOffset, real[] lengthscale_idAge, real[] lengthscale_sharedAge, real[] lengthscale_diseaseAge, real[] lengthscale_continuous, real[] lengthscale_categAge, real[] warp_steepness, vector[] beta){
  
  int n1 = num_elements(X1[1]);
  int n2 = num_elements(X2[1]);
  real x1_age[n1] = to_array_1d(X1[2]);   // age covariate as an array
  real x2_age[n2] = to_array_1d(X2[2]);   // age covariate as an array
  int sum_D = sum(D);
  matrix[n1,n2] KX[sum_D];
  int r = 0;
  if(D[1]==1){
    real alp = alpha_idAge[1];
    real ell = lengthscale_idAge[1];
    r += 1;
    KX[r] = cov_exp_quad(x1_age, x2_age, alp, ell) .* KF[1];
  }
  if(D[2]==1){
    real alp = alpha_sharedAge[1];
    real ell = lengthscale_sharedAge[1];
    r += 1;
    KX[r] = cov_exp_quad(x1_age, x2_age, alp, ell);
  }
  if(D[3]==1){
    real alp = alpha_diseaseAge[1];
    real ell = lengthscale_diseaseAge[1];
    real stp = warp_steepness[1];
    vector[n1] x1_tilde;
    vector[n2] x2_tilde;
    real w1[n1];
    real w2[n2];
    r += 1;

    // Handle diseaseAge uncertainty
    if(UNCRT==0){
      x1_tilde = X1[3];
      x2_tilde = X2[3];
    }else{
      x1_tilde = STANFUNC_get_x_tilde(X1[3], T_onset[1], T_observed, caseID_to_rows_1, caseID_nrows_1);
      x2_tilde = STANFUNC_get_x_tilde(X2[3], T_onset[1], T_observed, caseID_to_rows_2, caseID_nrows_2);
    }
    w1 = to_array_1d(STANFUNC_warp_input(x1_tilde, stp, 0.0, 1.0));
    w2 = to_array_1d(STANFUNC_warp_input(x2_tilde, stp, 0.0, 1.0));

    // Create disease effect kernel
    KX[r] = KF[2] .* cov_exp_quad(w1, w2, alp, ell);
    if(HMGNS==0){
      KX[r] = STANFUNC_K_beta(beta[1], row_to_caseID_1, row_to_caseID_2) .* KX[r];
    }
    if(USE_VAR_MASK==1){
      KX[r] = STANFUNC_K_var_mask(x1_tilde, x2_tilde, stp, vm_params) .* KX[r];
    }
  }
  for(j in 1:D[4]){
    real alp = alpha_continuous[j];
    real ell = lengthscale_continuous[j];
    int ix  = 2 + D[3] + j;
    r += 1;
    KX[r] = cov_exp_quad(to_array_1d(X1[ix]), to_array_1d(X2[ix]), alp, ell);
  }
  for(j in 1:D[5]){
    real alp = alpha_categAge[j];
    real ell = lengthscale_categAge[j];
    int ikf = 1 + D[3] + j;
    r += 1;
    KX[r] = cov_exp_quad(x1_age, x2_age, alp, ell) .* KF[ikf];
  }
  for(j in 1:D[6]){
    int ikf = 1 + D[3] + D[5] + j;
    real alp = alpha_categOffset[j];
    r += 1;
     KX[r] = square(alp) * KF[ikf];
  }
  return(KX);
}
  
  
// COMPUTE ALL KERNEL MATRICES (SYMMETRIC)
matrix[] STANFUNC_compute_kernel_matrices_symmetric(vector[] X, int[,] caseID_to_rows, int[] row_to_caseID, int[] caseID_nrows, matrix[] KF, vector[] T_onset, vector T_observed, int[] D, int UNCRT, int HMGNS, int USE_VAR_MASK, real[] vm_params, real[] alpha_idAge, real[] alpha_sharedAge, real[] alpha_diseaseAge, real[] alpha_continuous, real[] alpha_categAge, real[] alpha_categOffset, real[] lengthscale_idAge, real[] lengthscale_sharedAge, real[] lengthscale_diseaseAge, real[] lengthscale_continuous, real[] lengthscale_categAge, real[] warp_steepness, vector[] beta){
    
  int n = num_elements(X[1]);
  real x_age[n] = to_array_1d(X[2]);   // age covariate as an array
  int sum_D = sum(D);
  matrix[n,n] KX[sum_D];
  int r = 0;
  if(D[1]==1){
    real alp = alpha_idAge[1];
    real ell = lengthscale_idAge[1];
    r += 1;
    KX[r] = cov_exp_quad(x_age, alp, ell) .* KF[1];
  }
  if(D[2]==1){
    real alp = alpha_sharedAge[1];
    real ell = lengthscale_sharedAge[1];
    r += 1;
    KX[r] = cov_exp_quad(x_age, alp, ell);
  }
  if(D[3]==1){
    real alp = alpha_diseaseAge[1];
    real ell = lengthscale_diseaseAge[1];
    real stp = warp_steepness[1];
    vector[n] x_tilde;
    real w[n];
    r += 1;

    // Handle diseaseAge uncertainty
    if(UNCRT==0){
      x_tilde = X[3];
    }else{
      x_tilde = STANFUNC_get_x_tilde(X[3], T_onset[1], T_observed, caseID_to_rows, caseID_nrows);
    }
    w = to_array_1d(STANFUNC_warp_input(x_tilde, stp, 0.0, 1.0));

    // Create disease effect kernel
    KX[r] = KF[2] .* cov_exp_quad(w, alp, ell);
    if(HMGNS==0){
      KX[r] = STANFUNC_K_beta_symmetric(beta[1], row_to_caseID) .* KX[r];
    }
    if(USE_VAR_MASK==1){
      KX[r] = STANFUNC_K_var_mask(x_tilde, x_tilde, stp, vm_params) .* KX[r];
    }
  }
  for(j in 1:D[4]){
    real alp = alpha_continuous[j];
    real ell = lengthscale_continuous[j];
    int ix  = 2 + D[3] + j;
    r += 1;
    KX[r] = cov_exp_quad(to_array_1d(X[ix]), alp, ell);
  }
  for(j in 1:D[5]){
    real alp = alpha_categAge[j];
    real ell = lengthscale_categAge[j];
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
