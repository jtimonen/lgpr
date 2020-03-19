
// INPUT WARP
vector STAN_warp_input(vector x, real a){
  return( -1 + 2*inv(1+exp(-a*x)) );
}

// VARIANCE MASK
vector STAN_var_mask(vector x, real a){
  return( inv(1+exp(-a*x)) );
}

// COMPUTE X_TILDE
vector STAN_get_x_tilde(vector x_disAge, vector T_effect, vector T_observed, int[,] mapping, int[] lengths){
  int n_tot = num_elements(x_disAge);
  int N_cases = num_elements(lengths);
  vector[n_tot] x_tilde = rep_vector(0.0, n_tot);
  for(k in 1:N_cases){
    int inds[lengths[k]] = mapping[k, 1:lengths[k]];
    x_tilde[inds] = x_disAge[inds] + T_observed[k] - T_effect[k];
  }
  return(x_tilde);
}

// CATEGORICAL ZERO-SUM KERNEL
matrix STAN_K_zerosum(vector x1, vector x2, int M){
  int n1 = num_elements(x1);
  int n2 = num_elements(x2);
  matrix[n1,n2] K;
  for(i in 1:n1){
    for(j in 1:n2){
      if(x1[i]==x2[j]){
        K[i,j] = 1;
      }else{
       K[i,j] = -inv(M-1); 
      }
    }
  }
  return(K);
}


// BINARY (MASK) KERNEL
matrix STAN_K_bin(int[] x1, int[] x2, int c){
  int n1 = num_elements(x1);
  int n2 = num_elements(x2);
  matrix[n1,n2] K;
  for(i in 1:n1){
    for(j in 1:n2){
      K[i,j] = (x1[i]==c)*(x2[j]==c);
    }
  }
  return(K);
}


// Multiplier matrix to enable variance masking
matrix STAN_K_var_mask(vector x_tilde, real stp, real[] vm_params){
  int n = num_elements(x_tilde);
  real a = stp * vm_params[2];
  real r = inv(a)*logit(vm_params[1]);
  vector[n] s = STAN_var_mask(x_tilde - r, a);
  matrix[n,n] MASK = tcrossprod(to_matrix(s));
  return(MASK);
}


// Multiplier matrix to enable heterogeneous diseaseAge effect
matrix STAN_K_beta(vector beta, int[] row_to_caseID){
  int n = num_elements(row_to_caseID);
  int i_caseID = 0;
  int j_caseID = 0;
  real tmp;
  matrix[n,n] BETA;
  for(i in 1:(n-1)){
    i_caseID = row_to_caseID[i];
    for(j in (i+1):n){
      j_caseID = row_to_caseID[j];
      if(i_caseID*j_caseID > 0){
        tmp = sqrt(beta[i_caseID]*beta[j_caseID]);
      }else{
        tmp = 0;
      }
      BETA[i,j] = tmp;
      BETA[j,i] = tmp;
    }
    if(i_caseID > 0){
      BETA[i,i] = beta[i_caseID];
    }else{
      BETA[i,i] = 0;
    }
  }
  i_caseID = row_to_caseID[n];
  if(i_caseID > 0){
    BETA[n,n] = beta[i_caseID];
  }else{
    BETA[n,n] = 0;
  }
  return(BETA);
}
