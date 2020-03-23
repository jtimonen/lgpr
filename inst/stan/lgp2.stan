// *lgp2.stan*
#include info/header.stan
#include info/license.stan

functions{
#include common/prior.stan
#include common/kernels.stan
#include common/matrices.stan
}

data {
#include common/data.stan
  int<lower=0,upper=4> LH;  // observation model
  int t_SIG[2];             // for Gaussian noise std
  real p_SIG[3];            // for Gaussian noise std
  int t_PHI[2];             // for precision parameter phi
  real p_PHI[3];            // for precision parameter phi
  int<lower=1> N_trials[n]; // #trials (ones for bernoulli model)
}

transformed data{
#include chunks/transformed_data.stan
}

parameters {
#include chunks/parameters.stan

// Miscellaneous
real<lower=0> sigma_n[LH==1 || LH==0];  // noise std for Gaussian likelihood
vector[n] ETA[F_IS_SAMPLED, sum_D];     // isotropic versions of F
real<lower=0> phi[LH==3 || LH==0];      // parameter for NB likelihood


}

transformed parameters {
#include chunks/transformed_parameters.stan

vector[n] F[sum_D];
matrix[n,n] KX[sum_D] = STAN_kernels(X, caseID_to_rows, row_to_caseID, caseID_nrows, KF, T_effect, T_observed, D, UNCRT, HMGNS, USE_VAR_MASK, vm_params, alpha_idAge, alpha_sharedAge,  alpha_diseaseAge, alpha_continuous, alpha_categAge, alpha_categOffset, ell_idAge, ell_sharedAge, ell_diseaseAge, ell_continuous, ell_categAge, warp_steepness, beta);
for(r in 1:sum_D){
  matrix[n,n] EYE = diag_matrix(rep_vector(DELTA, n));
  matrix[n,n] Lxr = cholesky_decompose(KX[r] + EYE);
  F[1,r,] = Lxr*ETA[1,r,];
}

}

model {
  // Isotropic normals for auxiliary variables when F is sampled
  if(F_IS_SAMPLED){
    for(j in 1:sum_D){
      target += normal_lpdf(ETA[1,j] | 0, 1);
    }
  }

  // Noise level parameters
  if(LH==1){
    target += STAN_log_prior(sigma_n[1], t_SIG[1:2], p_SIG[1:3]);
  }else if(LH==3){
    target += STAN_log_prior(phi[1], t_PHI[1:2], p_PHI[1:3]);
  }else if(LH==0){
    target += STAN_log_prior(sigma_n[1], t_SIG[1:2], p_SIG[1:3]);
    target += STAN_log_prior(phi[1], t_PHI[1:2], p_PHI[1:3]);
  }

#include chunks/model_priors.stan
#include chunks/model_likelihood.stan
}

generated quantities {
#include chunks/generated_quantities.stan
}

