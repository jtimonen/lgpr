
vector[N_cases] T_effect[UNCRT];
vector[n] F[F_IS_SAMPLED, sum_D];
if(UNCRT){
  T_effect[1] = L_ons[1] + (U_ons[1] - L_ons[1]) .* T_raw[1];
}
if(F_IS_SAMPLED){
  
  matrix[n,n] KX[sum_D] = STAN_compute_kernel_matrices(X, caseID_to_rows,
      row_to_caseID, caseID_nrows, KF, T_effect, T_observed, D, UNCRT, HMGNS,
      USE_VAR_MASK, vm_params, alpha_idAge, alpha_sharedAge,  alpha_diseaseAge,
      alpha_continuous, alpha_categAge, alpha_categOffset,
      ell_idAge, ell_sharedAge, ell_diseaseAge, ell_continuous, ell_categAge,
      warp_steepness, beta);
      
  for(r in 1:sum_D){
    matrix[n,n] EYE = diag_matrix(rep_vector(DELTA, n));
    matrix[n,n] Lxr = cholesky_decompose(KX[r] + EYE);
    F[1,r,] = Lxr*ETA[1,r,];
  }
  
}
