// Categorical kernel (not evaluated during sampling)
matrix K_cat(vector x1, vector x2){
  int n1 = num_elements(x1);
  int n2 = num_elements(x2);
  matrix[n1,n2] K;
  for(i in 1:n1){
    for(j in 1:n2){
      K[i,j] = x1[i]==x2[j];
    }
  }
  return(K);
}

// Binary (mask) kernel (not evaluated during sampling)
// - this is currently needed only for X_notnan
matrix K_bin(int[] x1, int[] x2, int c){
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

// Multiplier matrix to enable heterogeneous diseaseAge effect
matrix K_var_mask(vector x_tilde, real stp, real[] vm_params){
  int n = num_elements(x_tilde);
  real a = stp * vm_params[2];
  real h = vm_params[1];
  real r = 1/a*log(h/(1-h));
  vector[n] s = var_mask(x_tilde - r, a);
  matrix[n,n] MASK = tcrossprod(to_matrix(s));
  return(MASK);
}

// Multiplier matrix to enable heterogeneous diseaseAge effect
matrix K_beta(vector beta, int[] row_to_caseID){
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
