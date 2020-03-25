// *model2.stan*
#include comments/header.stan
#include comments/license.stan

functions{
#include functions/utils.stan
#include functions/prior.stan
#include functions/kernels.stan
}

data {
#include chunks/data.stan
  int       y_int[n];       // the response variable (as an array of integers)
  vector[n] C_hat;          // GP mean vector
  int<lower=1,upper=4> LH;  // observation model
  int t_SIG[2];             // for Gaussian noise std
  real p_SIG[3];            // for Gaussian noise std
  int t_PHI[2];             // for precision parameter phi
  real p_PHI[3];            // for precision parameter phi
  int<lower=1> N_trials[n]; // #trials (ones for bernoulli model)
}

transformed data{
  real x_age[n] = to_array_1d(X[2]);
  int sum_D = sum(D);
  int nf = 1 + D[3] + D[5] + D[6];
  matrix[n,n] KF[nf];
#include covariance/fixed.stan
#include chunks/info.stan
}

parameters {
  real<lower=0> sigma_n[LH==1]; // noise std for Gaussian likelihood
  real<lower=0> phi[LH==3];     // parameter for NB likelihood
  vector[n] ETA[sum_D];         // isotropic versions of F
#include chunks/parameters.stan
}

transformed parameters {
  matrix[n,n] KX[sum_D];
  vector[n] F[1,sum_D];
  vector[N_cases] T_effect[UNCRT];
  if(UNCRT){ T_effect[1] = L_ons[1] + (U_ons[1] - L_ons[1]) .* T_raw[1];}
#include covariance/additive.stan
  for(j in 1:sum_D){
    matrix[n,n] EYE = diag_matrix(rep_vector(DELTA, n));
    F[1,j] = cholesky_decompose(KX[j] + EYE) * ETA[j];
  }
}

model {
  // Prior
  for(j in 1:sum_D){ target += normal_lpdf(ETA[j] | 0, 1); }
  if(LH==1){
    target += STAN_log_prior(sigma_n[1], t_SIG[1:2], p_SIG[1:3]);
  }else if(LH==3){
    target += STAN_log_prior(phi[1], t_PHI[1:2], p_PHI[1:3]);
  }
#include chunks/priors.stan

  // Likelihood
  if(SKIP_LH==0){
    vector[n] F_ss = C_hat;
    for (i in 1:n){
      F_ss[i] += sum(F[1,,i]);
    }
#include likelihood/sampled.stan
  }
}
