  // Precompute fixed kernel matrices
  matrix[N, N] K_const[J] = STAN_kernel_const_all(
    N, N, Z, Z, X_mask, X_mask, Z_M, components
  );

  // Delta vector for diagonal jitter
  vector[N] delta_vec = rep_vector(delta, N);
