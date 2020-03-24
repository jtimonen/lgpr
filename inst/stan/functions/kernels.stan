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

// DISEASE KERNEL
matrix STAN_K_disease(vector x, matrix Kf, real[] alpha, real[] ell, real[] steepness, real[] vm_params, vector[] beta, int USE_VAR_MASK, int UNCRT, int HMGNS, vector[] T_effect, vector T_observed, int[,] caseID_to_rows, int[] row_to_caseID_plus1, int[] caseID_nrows){
  
  // Setup
  int n = num_elements(x);
  vector[n] x_tilde;
  real w[n];
  matrix[n,n] Kx;

  // Handle diseaseAge uncertainty
  if(UNCRT==0){
    x_tilde = x;
  }else{
    x_tilde = STAN_get_x_tilde(x, T_effect[1], T_observed, caseID_to_rows, caseID_nrows);
  }
  
  // Input warping
  w = to_array_1d(STAN_warp_input(x_tilde, steepness[1]));

  // Create disease effect kernel
  Kx = Kf .* cov_exp_quad(w, alpha[1], ell[1]);
  if(HMGNS==0){ Kx = STAN_K_beta(beta[1], row_to_caseID_plus1) .* Kx; }
  if(USE_VAR_MASK==1){ Kx = STAN_K_var_mask(x_tilde, steepness[1], vm_params) .* Kx;}
  return(Kx);
}
