
  // Kernel hyperparameters for idAge component
  if(D[1]==1){
    target += log_prior(alpha_idAge[1],       t_ID[1,1:2], p_ID[1,1:3]);
    target += log_prior(lengthscale_idAge[1], t_ID[1,3:4], p_ID[1,4:6]);
  }

  // Kernel hyperparameters for sharedAge component
  if(D[2]==1){
    target += log_prior(alpha_sharedAge[1],       t_A[1,1:2], p_A[1,1:3]);
    target += log_prior(lengthscale_sharedAge[1], t_A[1,3:4], p_A[1,4:6]);
  }

  // Kernel hyperparameters for diseaseAge component
  if(D[3]==1){
    target += log_prior(alpha_diseaseAge[1],       t_D[1,1:2], p_D[1,1:3]);
    target += log_prior(lengthscale_diseaseAge[1], t_D[1,3:4], p_D[1,4:6]);
    target += log_prior(warp_steepness[1],         t_D[1,5:6], p_D[1,7:9]);
  }

  // Kernel hyperparameters for other continuous components
  for(j in 1:D[4]){
    target += log_prior(alpha_continuous[j],       t_CNT[j,1:2], p_CNT[j,1:3]);
    target += log_prior(lengthscale_continuous[j], t_CNT[j,3:4], p_CNT[j,4:6]);
  }

  // Kernel hyperparameters for categorical * age components
  for(j in 1:D[5]){
    target += log_prior(alpha_categAge[j],       t_CAT[j,1:2], p_CAT[j,1:3]);
    target += log_prior(lengthscale_categAge[j], t_CAT[j,3:4], p_CAT[j,4:6]);
  }

  // Kernel hyperparameters for categorical offset components
  for(j in 1:D[6]){
    target += log_prior(alpha_categOffset[j], t_OFS[j,1:2], p_OFS[j,1:3]);
  }

  // Noise level parameters
  if(LH==1){
    target += log_prior(sigma_n[1], t_SIG[1:2], p_SIG[1:3]);
  }else if(LH==3){
    target += log_prior(phi[1], t_PHI[1:2], p_PHI[1:3]);
  }else if(LH==0){
    target += log_prior(sigma_n[1], t_SIG[1:2], p_SIG[1:3]);
    target += log_prior(phi[1],     t_PHI[1:2], p_PHI[1:3]);
  }

  // DiseaseAge uncertainty prior
  if(UNCRT){
    real tx;
    for(k in 1:N_cases){
      if(relative==1){
        tx = - T_observed[k] + T_onset[1,k];
      }else if(backwards==1){
        tx = T_observed[k] - T_onset[1,k];
      }else{
        tx = T_onset[1,k];
      }
      target += log_prior(tx, t_ONS[k,1:2], p_ONS[k,1:3]);
    }
  }
  
  // Heterogeneous disease effect parameters
  if(HMGNS==0){
    beta[1] ~ beta(p_BET[1], p_BET[2]);
  }
  
  // Isotropic normals
  if(F_is_sampled){
    for(j in 1:sum_D){
      ETA[1,j] ~ normal(0,1);
    }
  }
  
