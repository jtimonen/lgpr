// *lgp1.stan*
#include comments/header.stan
#include comments/license.stan

functions{
#include common/prior.stan
#include common/kernels.stan
#include common/matrices.stan
}

data {
#include common/data.stan
}

transformed data{
  int nf = 1 + D[3] + D[5] + D[6];
  int sum_D = sum(D);
  matrix[n,n] KF[nf] = STAN_kernels_fixed(X, X_notnan, D, N_tot, N_cat);
  int DO_GEN_QUANT = (1 - SKIP_GQ);
#include common/info.stan
}

parameters {
#include common/parameters.stan
}

transformed parameters {
#include common/effect_time.stan
}

model {
#include chunks/model_priors.stan
  if(LH!=0){
#include chunks/model_likelihood.stan
  }
}


generated quantities {
#include chunks/generated_quantities.stan
}

