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
// - TODO: add explanation
matrix[] STAN_phi_matrix(data vector[] x, data int M, data real L, 
    data int[,] components){
  int n = num_elements(x[1]);
  int J = size(components);
  matrix[n, M] PHI[J];
  for(j in 1:J) {
    matrix[n, M] PHI_j = rep_matrix(1.0, n, M); // is this ok?
    if(components[j,1] > 0) {
      int idx_cont = components[j,9];
      for(m in 1:M) {
        // TODO: use warp-transformed x_cont_unnorm for some components
        PHI_j[:,m] = STAN_phi(x[idx_cont], m, L);
      }
      PHI[j] = PHI_j;
    }
  }
  return(PHI);
}

// Compute diagonals of diagonal matrices Lambda
vector[] STAN_lambda_matrix(real[] alpha, real[] ell, data int M, data real L,
    data int[,] components){
  int J = size(alpha);
  vector[M] Lambda[J];
  int j_ell = 0;
  for(j in 1:J) {
    vector[M] Lambda_j = rep_vector(1.0, M);
    if (components[j,1] > 0) {
      j_ell = j_ell + 1;
      for(m in 1:M) {
        real w = (m*pi())/(2.0*L);
        Lambda_j[m] = STAN_spd_eq(w, alpha[j], ell[j_ell]);
      }
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

// Create D
// - bfa_lambda = array of vectors with shape [num_basisfun], length num_comps
// - bfa_delta = a vector of shape [R]
vector STAN_D_matrix(real[] alpha, vector[] bfa_lambda, 
    vector bfa_delta, int[] ranks) {
  int M = num_elements(bfa_lambda[1]);
  int RM = sum(ranks)*M;
  int J = size(ranks);
  vector[RM] D_diag;
  int q0 = 0;
  int k0 = 0;
  for(j in 1:J) {
    // Create D_j
    int r = ranks[j];
    if (r > 0) {
      for (k in 1:r) {
        for (m in 1:M) {
          int q = m + (k-1)*M;
          D_diag[q0 + q] = alpha[j]*bfa_delta[k0+k]*bfa_lambda[j][m];
        }
      }
    }
    k0 = k0 + r;
    q0 = q0 + r*M;
  }
  return D_diag;
}

// Create V
// - bfa_phi = array of matrices with shape [num_obs, num_basisfun], length 
//     num_comps
// - bfa_theta = matrix of shape [num_obs, R]
matrix STAN_V_matrix(matrix[] bfa_phi, matrix bfa_theta, int[] ranks) {
  int J = size(bfa_phi);
  int n = rows(bfa_phi[1]);
  int M = cols(bfa_phi[1]);
  int RM = sum(ranks)*M;
  matrix[n, RM] V = rep_matrix(1.0, n, RM);
  int q0 = 0;
  int k0 = 0;
  for(j in 1:J) {
    // Create V_j
    int r = ranks[j];
    if (r > 0) {
      for (k in 1:r) {
        for (m in 1:M) {
          int q = m + (k-1)*M;
          V[:,q0 + q] = bfa_theta[:,k0+k] .* bfa_phi[j][:,m];
        }
      }
    }
    k0 = k0 + r;
    q0 = q0 + r*M;
  }
  return V;
}

