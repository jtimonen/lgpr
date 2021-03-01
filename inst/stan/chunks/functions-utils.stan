// Compute total signal by summing components
vector STAN_vectorsum(vector[] vecs, data int L){
  int num_vecs = size(vecs);
  vector[L] s = rep_vector(0, L);
  for (j in 1:num_vecs){
    s += vecs[j];
  }
  return(s);
}

// Input warping function
vector STAN_warp_input(vector x, real a){
  return( -1 + 2*inv(1+exp(-a*x)) );
}

// Variance masking function
vector STAN_var_mask(vector x, real a){
  return( inv(1+exp(-a*x)) );
}

// Expand a vector
vector STAN_expand(vector v, data int[] idx_expand){
  int L = num_elements(v);
  vector[L+1] v_add0 = rep_vector(0.0, L+1);
  v_add0[2:(L+1)] = v;
  return(v_add0[idx_expand]);
}

// Edit a continuous covariate according to sampled uncertainty
vector STAN_edit_x_cont(
  vector x_cont,
  data int[] idx_expand,
  data vector teff_obs,
  vector teff)
{
  int n = num_elements(x_cont);
  vector[n] x_teff_obs = STAN_expand(teff_obs, idx_expand);
  vector[n] x_teff = STAN_expand(teff, idx_expand);
  return(x_cont + x_teff_obs - x_teff);
}

// Repeat vector x, J times
vector STAN_rep_vector_times(vector x, int J) {
  int N = num_elements(x);
  vector[J*N] y = rep_vector(0.0, J*N);
  int idx = 1;
  for(j in 1:J) {
    y[idx:(idx+N-1)] = x;
    idx = idx + N;
  }
  return(y);
}

// Repeat each element of vector x, J times
vector STAN_rep_vector_each(vector x, int J) {
  int N = num_elements(x);
  vector[J*N] y = rep_vector(0.0, J*N);
  int idx = 1;
  for(n in 1:N) {
    y[idx:(idx+J-1)] = rep_vector(x[n], J);
    idx = idx + J;
  }
  return(y);
}

// Repeat each column of x, J times
matrix STAN_rep_cols_times(matrix X, int J) {
  int R = rows(X);
  int N = cols(X);
  matrix[R, J*N] Y = rep_matrix(0.0, R, J*N);
  int idx = 1;
  for(j in 1:J) {
    Y[:, idx:(idx+N-1)] = X;
    idx = idx + N;
  }
  return(Y);
}

// Repeat each element of vector x, J times
matrix STAN_rep_cols_each(matrix X, int J) {
  int R = rows(X);
  int N = cols(X);
  matrix[R, J*N] Y = rep_matrix(0.0, R, J*N);
  int idx = 1;
  for(n in 1:N) {
    for(j in 1:J) {
        Y[:, idx+j-1] = X[:,n];
    }
    idx = idx + J;
  }
  return(Y);
}

// Compute the quadractic form x * inv(A) * x^T
real STAN_quad_form_inv(vector x, matrix A){
  int n = num_elements(x);
  matrix[n, n] L = cholesky_decompose(A);
  vector[n] a = mdivide_left_tri_low(L, x);
  return dot_self(a);
}
