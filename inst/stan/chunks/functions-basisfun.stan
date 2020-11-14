// EIGENFUNCTIONS OF DIRICHLET PROBLEM
vector STAN_phi(vector x, int m, data real L){
  real A = inv(sqrt(L));
  real B = pi()*m/(2.0*L);
  return(A*sin(B*(x+L)));
}

// EIGENVALUES OF DIRICHLET PROBLEM
real STAN_lambda(int m, data real L){
  real A = pi()*m/(2.0*L);
  return(square(A));
}

// SPECTRAL DENSITY OF EQ KERNEL
real STAN_spd_eq(real w, real alpha, real ell){
  real A = square(alpha)*ell*sqrt(2.0*pi());
  real B = 2.0*square(pi()*ell);
  return(A*exp(-B*square(w)));
}

// COMPUTE QUADRATIC FORM OF INVERSE
real STAN_quad_form_inv(vector x, matrix A){
  int n = num_elements(x);
  matrix[n, n] L = cholesky_decompose(A);
  vector[n] a = mdivide_left_tri_low(L, x);
  return dot_self(a);
}

// LOG PROB OF MULTIVARIATE NORMAL WITH LOW RANK COVARIANCE
real STAN_multi_normal_bfa_lpdf(vector y, matrix V, vector D_diag, real sigma){
  int n = num_elements(y);
  int RM = num_elements(D_diag);
  real t1 = n*log(2.0*pi());
  real t2; // log det
  real t3; // quadratic form
  real inv_s2 = inv_square(sigma);
  vector[RM] z = transpose(V)*y;
  matrix[RM,RM] Z = diag_matrix(inv(D_diag)) + inv_s2*crossprod(V);
  t2 = inv_s2*dot_self(y) + square(inv_s2)*STAN_quad_form_inv(z, Z);
  t3 = log_determinant(Z) + sum(log(D_diag)) + 2*n*log(sigma);
  return(-0.5*(t1 + t2 + t3));
}

// Compute basis function matrices Phi
matrix[] STAN_phi_matrix(data vector[] x, data int M, data real L){
  int n = num_elements(x[1]);
  int J = size(x);
  matrix[n, M] PHI[J];
  for(j in 1:J) {
    matrix[n, M] PHI_j = rep_matrix(0.0, n, M);
    for(m in 1:M) {
      // TODO: use warp-transformed x[j] for some components
      PHI_j[:,m] = STAN_phi(x[j], m, L);
    }
    PHI[j] = PHI_j;
  }
  return(PHI);
}

// Compute diagonals of diagonal matrices Lambda
vector[] STAN_lambda_matrix(real[] alpha, real[] ell, data int M, data real L){
  int J = size(alpha);
  vector[M] Lambda[J];
  for(j in 1:J) {
    vector[M] Lambda_j = rep_vector(0.0, M);
    for(m in 1:M) {
      // TODO: all components don't have ell
      real w = (m*pi())/(2.0*L);
      Lambda_j[m] = STAN_spd_eq(w, alpha[j], ell[j]);
    }
    Lambda[j] = Lambda_j;
  }
  return(Lambda);
}

// Get ranks of the constant kernel matrices
int[] STAN_ranks(data int[,] components, data int[] x_cat_num_levels){
  int J = size(components);
  int idx_cat[J] = components[,8];
  int ranks[J] = rep_array(0, J);
  for (j in 1:J) {
    int idx = idx_cat[j];
    if (idx > 0) {
      ranks[j] = x_cat_num_levels[idx] - 1; // TODO: prove this correct
    }
  }
  return(ranks);
}

// Compute diagonals of Delta
vector STAN_delta_matrix(data matrix[] K_const, data int[] ranks){
  int J = size(ranks);
  int n = cols(K_const[1]);
  int R = sum(ranks);
  vector[R] Delta = rep_vector(0.0, R);
  int idx = 1;
  for (j in 1:J) {
    int r = ranks[j];
    if (r > 0) {
      vector[n] lam = eigenvalues_sym(K_const[j]); // eigenvals in asc. order
      Delta[idx:(idx+r-1)] = lam[(n-r+1):n]; // only last r should be nonzero
    }
    idx = idx + r;
  }
  return(Delta);
}

// Compute eigenvector matrix Theta
matrix STAN_theta_matrix(data matrix[] K_const, data int[] ranks){
  int J = size(ranks);
  int n = cols(K_const[1]);
  int R = sum(ranks);
  matrix[n, R] Theta = rep_matrix(0.0, n, R);
  int idx = 1;
  for (j in 1:J) {
    int r = ranks[j];
    if (r > 0) {
      matrix[n,n] v = eigenvectors_sym(K_const[j]); // eigenvecs in asc. order
      Theta[:, idx:(idx+r-1)] = v[:, (n-r+1):n]; // take only last r vecs
    }
    idx = idx + r;
  }
  return(Theta);
}

