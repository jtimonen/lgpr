
vector[N_cases] T_onset[UNCRT];
vector[n] F[F_is_sampled, sum_D];
if(UNCRT){
  T_onset[1] = L_ons[1] + (U_ons[1] - L_ons[1]) .* T_raw[1];
}
if(F_is_sampled){
  matrix[n,n] KX[sum_D] = STANFUNC_compute_kernel_matrices_symmetric(X, caseID_to_rows, row_to_caseID, caseID_nrows, KF, T_onset, T_observed, D, UNCRT, HMGNS, USE_VAR_MASK, vm_params, alpha_idAge, alpha_sharedAge,  alpha_diseaseAge, alpha_continuous, alpha_categAge, alpha_categOffset, lengthscale_idAge, lengthscale_sharedAge, lengthscale_diseaseAge, lengthscale_continuous, lengthscale_categAge, warp_steepness, beta);
  for(r in 1:sum_D){
    matrix[n,n] EYE = diag_matrix(rep_vector(DELTA, n));
    matrix[n,n] Lxr = cholesky_decompose(KX[r] + EYE);
    F[1,r,] = Lxr*ETA[1,r,];
  }
}
