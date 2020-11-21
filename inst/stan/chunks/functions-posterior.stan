// Sum an array of vectors
vector STAN_vector_array_sum(vector[] k){
  int n = num_elements(k[1]);
  vector[n] k_sum = k[1];
  for(j in 2:size(k)){
    k_sum += k[j];
  }
  return(k_sum);
}

// Sum an array of matrices
matrix STAN_matrix_array_sum(matrix[] K){
  int n1 = rows(K[1]);
  int n2 = cols(K[1]);
  matrix[n1, n2] K_sum = K[1];
  for(j in 2:size(K)){
    K_sum += K[j];
  }
  return(K_sum);
}

// Take diagonal from each matrix in an array
vector[] STAN_matrix_array_diagonals(matrix[] K){
 int L = size(K);
 int n = rows(K[1]);
 vector[n] k[L];
 for(j in 1:L) {
   k[j] = diagonal(K[j]);
 }
 return(k);
}

// Helper for the below function
vector[] STAN_gp_posterior_helper(matrix Ly, matrix K_s, vector k_ss, vector v){
  int n = num_elements(v);
  int p = num_elements(k_ss);
  matrix[p, n] A = transpose(mdivide_left_tri_low(Ly, transpose(K_s)));
  vector[p] f_post[2];
  f_post[1] = A*v; // mean
  f_post[2] = sqrt(k_ss - rows_dot_self(A)); // std
  return(f_post);
}
// GP posterior with Gaussian observation model
vector[] STAN_gp_posterior(
  matrix[] KX,
  matrix[] KX_s,
  matrix[] KX_ss,
  data vector y,
  real delta,
  real sigma)
{
  
  // Declare variables
  int num_comps = size(KX);
  int n = rows(KX[1]);
  int p = rows(KX_s[1]);
  int J = num_comps + 1;
  int inds[2];
  vector[p] F_POST[2*J];
  vector[n] v;
  matrix[n, n] Ly;
  
  // Take only diagonals of each KX_ss
  vector[p] kX_ss[num_comps] = STAN_matrix_array_diagonals(KX_ss);
  
  // Compute Ky
  vector[n] delta_vec = rep_vector(delta, n);
  vector[n] sigma2_vec = rep_vector(square(sigma), n);
  matrix[n, n] Ky = STAN_matrix_array_sum(KX) + diag_matrix(sigma2_vec) +
    diag_matrix(num_comps*delta_vec);
  
  // Cholesky decomposition for Ky
  Ly = cholesky_decompose(Ky);
  v = mdivide_left_tri_low(Ly, y);
  
  // Component-wise means and stds
  for(j in 1:num_comps){
    inds = {j, J + j};
    F_POST[inds] = STAN_gp_posterior_helper(Ly, KX_s[j], kX_ss[j], v);
  }
  
  // Total mean and std
  inds = {J, 2*J};
  F_POST[inds] = STAN_gp_posterior_helper(
    Ly, STAN_matrix_array_sum(KX_s), STAN_vector_array_sum(kX_ss), v
  );
  return(F_POST);

}
