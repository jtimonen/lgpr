
vector[n] F_mean_cmp[1 - F_IS_SAMPLED, sum_D];
vector[n] F_var_cmp[ 1 - F_IS_SAMPLED, sum_D];
vector[n] F_mean_tot[1 - F_IS_SAMPLED];
vector[n] F_var_tot[ 1 - F_IS_SAMPLED];

if(F_IS_SAMPLED==0){
  matrix[n,n] A;
  vector[n] v;
  matrix[n,n] Ky;
  matrix[n,n] Ly;
  matrix[n,n] Kx = diag_matrix(rep_vector(DELTA, n));
  matrix[n,n] KX[sum_D] = STAN_compute_kernel_matrices(X, caseID_to_rows, row_to_caseID, caseID_nrows, KF, T_effect, T_observed, D, UNCRT, HMGNS, USE_VAR_MASK, vm_params, alpha_idAge, alpha_sharedAge,  alpha_diseaseAge, alpha_continuous, alpha_categAge, alpha_categOffset, ell_idAge, ell_sharedAge, ell_diseaseAge, ell_continuous, ell_categAge, warp_steepness, beta);
     
  for(j in 1:sum_D){
    Kx += KX[j];
  }
  Ky = Kx + diag_matrix(rep_vector(square(sigma_n[1]), n));
  Ly = cholesky_decompose(Ky);
  v = mdivide_left_tri_low(Ly, y);
  for(j in 1:sum_D){
    A  = mdivide_left_tri_low(Ly, transpose(KX[j]));
    F_mean_cmp[1,j] = transpose(A)*v;
    F_var_cmp[1,j] = diagonal(KX[j] - crossprod(A));
  }
  A = mdivide_left_tri_low(Ly, transpose(Kx));
  F_mean_tot[1] = transpose(A)*v;
  F_var_tot[1] = diagonal(Kx - crossprod(A));
}
