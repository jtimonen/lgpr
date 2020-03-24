// *model1.stan*
#include comments/header.stan
#include comments/license.stan

functions{
#include functions/utils.stan
#include functions/prior.stan
#include functions/kernels.stan
}

data {
#include chunks/data.stan
  int t_SIG[2];
  real p_SIG[3];
}

transformed data{
  vector[n] mu = rep_vector(0.0, n);
  real x_age[n] = to_array_1d(X[2]);
  int sum_D = sum(D);
#include covariance/fixed.stan
#include chunks/info.stan
}

parameters{
  real<lower=0> sigma_n;
#include chunks/parameters.stan
}

transformed parameters{
#include chunks/effect_time.stan
}

model{
  matrix[n,n] Ky = diag_matrix(rep_vector(DELTA + square(sigma_n), n));
  target += STAN_log_prior(sigma_n, t_SIG, p_SIG);
#include chunks/priors.stan
  if(SKIP_LH==0){
#include covariance/additive.stan
    for(j in 1:sum_D){ Ky += KX[j];}
    target += multi_normal_lpdf(y | mu, Ky);
  }
}

generated quantities{
  vector[n] F_mean_cmp[SKIP_GQ==0, sum_D];
  vector[n] F_var_cmp[SKIP_GQ==0, sum_D];
  vector[n] F_mean_tot[SKIP_GQ==0];
  vector[n] F_var_tot[SKIP_GQ==0];
  if(SKIP_GQ==0){
    matrix[n,n] A;
    vector[n] v;
    matrix[n,n] Ky;
    matrix[n,n] Ly;
    matrix[n,n] Kx = diag_matrix(rep_vector(DELTA, n));
#include covariance/additive.stan
    for(j in 1:sum_D){Kx += KX[j];}
    Ky = Kx + diag_matrix(rep_vector(square(sigma_n), n));
    Ly = cholesky_decompose(Ky);
    v = mdivide_left_tri_low(Ly, y);
    for(j in 1:sum_D){
      A  = mdivide_left_tri_low(Ly, transpose(KX[j]));
      F_mean_cmp[1,j] = transpose(A)*v;
      F_var_cmp[1,j] = diagonal(KX[j] - crossprod(A));
    }
    A = mdivide_left_tri_low(Ly, transpose(Kx));
    F_mean_tot[1] = transpose(A)*v;
    F_var_tot[1] = diagonal(Kx - crossprod(A));
  }
}
