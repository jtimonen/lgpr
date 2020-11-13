// BASIS FUNCTION
vector STAN_bfa_phi(vector x, int m, data real L){
  real A = inv(sqrt(L));
  real B = pi()*m/(2.0*L);
  return(A*sin(B*(x+L)));
}

// EIGENVALUE
real STAN_bfa_lambda(int m, data real L){
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
real STAN_bfa_multi_normal_lpdf(vector y, matrix V, vector D_diag, real sigma){
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
matrix[] STAN_bfa_Phi(data vector[] x, data int M, data real L){
  int n = num_elements(x[1]);
  int J = size(x);
  matrix[n, M] PHI[J];
  for(j in 1:J) {
    matrix[n, M] PHI_j = rep_matrix(0.0, n, M);
    for(m in 1:M) {
      // TODO: use warp-transformed x[j] for some components
      PHI_j[:,m] = STAN_bfa_phi(x[j], m, L);
    }
    PHI[j] = PHI_j;
  }
  return(PHI);
}

// Compute diagonals of diagonal matrices Lambda
vector[] STAN_bfa_Lambda(real[] alpha, real[] ell, data int M, data real L){
  int J = size(alpha);
  vector[M] LAMBDA[J];
  for(j in 1:J) {
    vector[M] LAMBDA_j = rep_vector(0.0, M);
    for(m in 1:M) {
      // TODO: all components don't have ell
      real w = (m*pi())/(2.0*L);
      LAMBDA_j[m] = STAN_spd_eq(w, alpha[j], ell[j]);
    }
    LAMBDA[j] = LAMBDA_j;
  }
  return(LAMBDA);
}

