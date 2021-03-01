// Evaluate m'th eigenfunction of the Dirichlet boundary value problem
// - x = points where to evaluate the function
// - m = index of the function
// - L = domain width
vector STAN_phi(vector x, int m, data real L){
  real A = inv(sqrt(L));
  real B = pi()*m/(2.0*L);
  return(A*sin(B*(x+L)));
}

// Evaluate m'th eigenvalue of the Dirichlet boundary value problem
// - m = index of the function
// - L = domain width
real STAN_lambda(int m, data real L){
  real A = pi()*m/(2.0*L);
  return(square(A));
}

// Evaluate spectral density of the exponentiated quadratic kernel
// - w = frequency
// - ell = lengthscale of the kernel
real STAN_spd_eq(real w, real ell){
  real A = ell*sqrt(2.0*pi());
  real B = - 0.5*square(w*ell);
  return(A*exp(B) + 1e-9);
}

// Compute num_comps basis function matrices Phi, which each have shape
// [num_obs, num_basisfun]
matrix[] STAN_phi_matrix(data vector[] x, data int M, data real L, 
    data int[,] components){
  int n = num_elements(x[1]);
  int J = size(components);
  matrix[n, M] PHI[J];
  for(j in 1:J) {
    matrix[n, M] PHI_j = rep_matrix(1.0, n, M); //
    if(components[j,1] > 0) {
      int idx_cont = components[j,9];
      for(m in 1:M) {
        // TODO: use warp-transformed x_cont_unnorm for some components
        PHI_j[:,m] = STAN_phi(x[idx_cont], m, L);
      }
    }
    PHI[j] = PHI_j;
  }
  return(PHI);
}

// Compute diagonals of diagonal matrices Lambda
vector[] STAN_lambda_matrix(real[] ell, data int M, data real L,
    data int[,] components){
  int J = size(components);
  vector[M] Lambda[J];
  int j_ell = 0;
  for(j in 1:J) {
    if (components[j,1] > 0) {
      j_ell = j_ell + 1;
      for(m in 1:M) {
        Lambda[j][m] = STAN_spd_eq((m*pi())/(2.0*L), ell[j_ell]);
      }
    } else {
      Lambda[j] = rep_vector(1.0/M, M);
    }
  }
  return(Lambda);
}

// Get ranks of the constant kernel matrices
int[] STAN_ranks(data int[,] components, data int[] x_cat_num_levels){
  int J = size(components);
  int ranks[J] = rep_array(1, J);
  for (j in 1:J) {
    int idx_cat = components[j,8];
    if (idx_cat > 0) {
      ranks[j] = x_cat_num_levels[idx_cat] - 1; // TODO: prove this correct
    }
  }
  return(ranks);
}

// Compute diagonals of Delta
vector STAN_delta_matrix(data matrix[] K_const, data int[] ranks, 
    data int[,] components){
  int J = size(ranks);
  int n = cols(K_const[1]);
  int R = sum(ranks);
  vector[R] Delta = rep_vector(1.0, R);
  int idx = 1;
  for (j in 1:J) {
    int r = ranks[j];
    int idx_cat = components[j,8];
    if (idx_cat > 0) {
      vector[n] lam = eigenvalues_sym(K_const[j]); // eigenvals in asc. order
      Delta[idx:(idx+r-1)] = lam[(n-r+1):n]; // only last r should be nonzero
      //print("DELTA, j = ", j, ", wrote to Delta[",idx,":",idx+r-1,"]");
    }
    idx = idx + r;
  }
  return(Delta);
}

// Compute eigenvector matrix Theta
matrix STAN_theta_matrix(data matrix[] K_const, data int[] ranks,
    data int[,] components){
  int J = size(ranks);
  int n = cols(K_const[1]);
  int R = sum(ranks);
  matrix[n, R] Theta = rep_matrix(1.0, n, R);
  int idx = 1;
  for (j in 1:J) {
    int r = ranks[j];
    int idx_cat = components[j,8];
    if (idx_cat > 0) {
      matrix[n,n] v = eigenvectors_sym(K_const[j]); // eigenvecs in asc. order
      Theta[:, idx:(idx+r-1)] = v[:, (n-r+1):n]; // take only last r vecs
      //print("THETA, j = ", j, ", wrote to Delta[",idx,":",idx+r-1,"]");
    }
    idx = idx + r;
  }
  return(Theta);
}

// Create D, a vector of length num_basisfun*sum(ranks)
// - bfa_lambda = a an array of num_comps vectors of shape [num_basisfun]
// - bfa_delta = a vector of shape [R]
vector STAN_D_matrix(real[] alpha, vector[] bfa_lambda, 
    vector bfa_delta, int[] ranks) {
  int M = num_elements(bfa_lambda[1]);
  int RM = sum(ranks)*M;
  int J = size(ranks);
  vector[RM] alpha_diag;
  vector[RM] lambda_diag;
  vector[RM] delta_diag = STAN_rep_vector_each(bfa_delta, M);
  int i1 = 1;
  int i2 = 1;
  for(j in 1:J) {
    int r = ranks[j]*M;
    i2 = i1 + r - 1;
    alpha_diag[i1:i2] = rep_vector(square(alpha[j]), r);
    lambda_diag[i1:i2] = STAN_rep_vector_times(bfa_lambda[j], ranks[j]);
    i1 = i1 + r;
  }
  return alpha_diag .* delta_diag .* lambda_diag;
}

// Create V, a matrix of shape [n, sum(ranks)*num_basisfun]
// - bfa_phi = array of matrices with shape [num_obs, num_basisfun], length 
//     num_comps
// - bfa_theta = matrix of shape [num_obs, R]
matrix STAN_V_matrix(matrix[] bfa_phi, matrix bfa_theta, int[] ranks) {
  int J = size(bfa_phi);
  int n = rows(bfa_phi[1]);
  int M = cols(bfa_phi[1]);
  int RM = sum(ranks)*M;
  matrix[n, RM] THETA = STAN_rep_cols_each(bfa_theta, M);
  matrix[n, RM] PHI;
  int i1 = 1;
  int i2 = 1;
  for(j in 1:J) {
    int r = ranks[j]*M;
    i2 = i1 + r - 1;
    PHI[:,i1:i2] = STAN_rep_cols_times(bfa_phi[j], ranks[j]);
    i1 = i1 + r;
  }
  return THETA .* PHI;
}

// Compute log density of multivariate normal N(y | 0, K), where
// - 0 = vector of zeros
// - K = V*diag(D_diag)*V^T + sigma^2*I
real STAN_multi_normal_bfa_logpdf(vector y, matrix V, vector D_diag, 
    real sigma){
  int n = num_elements(y);
  int RM = num_elements(D_diag);
  vector[RM] z = transpose(V)*y;
  matrix[RM,RM] Z = diag_matrix(square(sigma)*inv(D_diag)) + crossprod(V);
  real t1 = n*log(2.0*pi());
  real t2 = inv_square(sigma)*(dot_self(y) - STAN_quad_form_inv(z, Z));
  real t3 = 2*(n-RM)*log(sigma) + log_determinant(Z) + sum(log(D_diag));
  return(-0.5*(t1 + t2 + t3));
}
