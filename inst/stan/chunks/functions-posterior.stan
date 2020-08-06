// Helper for the below function
vector[] STAN_gp_posterior_helper(matrix Ly, matrix K, vector v){
  int n = num_elements(v);
  matrix[n, n] A = mdivide_left_tri_low(Ly, transpose(K));
  vector[n] f_post[2];
  f_post[1] = transpose(A)*v;
  f_post[2] = diagonal(K - crossprod(A));
  return(f_post);
}

// GP posterior with Gaussian observation model
vector[] STAN_gp_posterior(
  matrix[] KX,
  data vector y,
  real delta,
  real sigma)
{
  
  // Declare variables
  int num_comps = size(KX);
  int n = num_elements(y);
  int J = num_comps + 1;
  int inds[2];
  vector[n] F_POST[2*J];
  vector[n] v;
  matrix[n, n] Ky;
  matrix[n, n] Ly;
  
  // Compute total sum kernel
  matrix[n, n] Kx = rep_matrix(0.0, n, n);
  vector[n] delta_vec = rep_vector(delta, n);
  vector[n] sigma2_vec = rep_vector(square(sigma), n);
  for(j in 1:num_comps){
    Kx += KX[j];
  }
  Ky = Kx + diag_matrix(sigma2_vec) + diag_matrix(num_comps*delta_vec);
  
  // Cholesky decomposition for total sum kernel
  Ly = cholesky_decompose(Ky);
  v = mdivide_left_tri_low(Ly, y);
  
  // Component-wise means and variances
  for(j in 1:num_comps){
    inds = {j, J + j};
    F_POST[inds] = STAN_gp_posterior_helper(Ly, KX[j], v);
  }
  
  // Total mean and variance
  inds = {J, 2*J};
  F_POST[inds] = STAN_gp_posterior_helper(Ly, Kx, v);
  return(F_POST);

}
