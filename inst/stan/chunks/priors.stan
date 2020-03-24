// Kernel hyperparameters for idAge component
if(D[1]==1){
  target += STAN_log_prior(alpha_idAge[1], t_ID[1,1:2], p_ID[1,1:3]);
  target += STAN_log_prior(ell_idAge[1], t_ID[1,3:4], p_ID[1,4:6]);
}

// Kernel hyperparameters for sharedAge component
if(D[2]==1){
  target += STAN_log_prior(alpha_sharedAge[1], t_A[1,1:2], p_A[1,1:3]);
  target += STAN_log_prior(ell_sharedAge[1], t_A[1,3:4], p_A[1,4:6]);
}

// Kernel hyperparameters for diseaseAge component
if(D[3]==1){
  target += STAN_log_prior(alpha_diseaseAge[1], t_D[1,1:2], p_D[1,1:3]);
  target += STAN_log_prior(ell_diseaseAge[1], t_D[1,3:4], p_D[1,4:6]);
  target += STAN_log_prior(warp_steepness[1], t_D[1,5:6], p_D[1,7:9]);
}

// Kernel hyperparameters for other continuous components
for(j in 1:D[4]){
  target += STAN_log_prior(alpha_continuous[j], t_CNT[j,1:2], p_CNT[j,1:3]);
  target += STAN_log_prior(ell_continuous[j], t_CNT[j,3:4], p_CNT[j,4:6]);
}

// Kernel hyperparameters for categorical * age components
for(j in 1:D[5]){
  target += STAN_log_prior(alpha_categAge[j], t_CAT[j,1:2], p_CAT[j,1:3]);
  target += STAN_log_prior(ell_categAge[j], t_CAT[j,3:4], p_CAT[j,4:6]);
}

// Kernel hyperparameters for categorical offset components
for(j in 1:D[6]){
  target += STAN_log_prior(alpha_categOffset[j], t_OFS[j,1:2], p_OFS[j,1:3]);
}

// DiseaseAge uncertainty prior
if(UNCRT){
  real tx;
  for(k in 1:N_cases){
    if(RELATIVE==1){
      tx = - T_observed[k] + T_effect[1,k];
    }else if(BACKWARDS==1){
      tx = T_observed[k] - T_effect[1,k];
    }else{
      tx = T_effect[1,k];
    }
    target += STAN_log_prior(tx, t_ONS[k,1:2], p_ONS[k,1:3]);
  }
}

// Heterogeneous disease effect parameters
if(HMGNS==0){ target += beta_lpdf(beta[1] | p_BET[1], p_BET[2]); }
