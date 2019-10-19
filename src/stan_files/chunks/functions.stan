functions{

  // LOG PRIOR TO BE ADDED TO TARGET
  real STANFUNC_log_prior(real x, int[ ] types, real[ ] hp){
    real lp;
    real a = hp[1]; // prior hyperparameter 1
    real b = hp[2]; // prior hyperparameter 2
    real c = hp[3]; // prior hyperparameter 3 (currently not used)
    real theta;

    // Possible transform and log of its absolute derivative
    if (types[2]==0){
      lp = 0;
      theta = x;
    }else if (types[2]==1){
      lp = log(fabs(2*x));
      theta = square(x);
    }else{
      reject("invalid value of types[2]!")
    }
  
    // Value of pdf
    if (types[1]==1){
      // do nothing
    }else if (types[1]==2){
      lp += normal_lpdf(theta|a,b);
    }else if (types[1]==3){
      lp += student_t_lpdf(theta|a,0,b);
    }else if (types[1]==4){
      lp += gamma_lpdf(theta|a,b);
    }else if (types[1]==5){
      lp += inv_gamma_lpdf(theta|a,b);
    }else if (types[1]==6){
      lp += lognormal_lpdf(theta|a,b);
    }else{
      reject("types[1] must be an integer between 1 and 6; found =", types[1]);
    }
    
    return(lp);
  }
  
  // INPUT WARP
  vector STANFUNC_warp_input(vector x, real a, real b, real c){
    int n_tot = num_elements(x);
    vector[n_tot] w = 2*c*(-0.5 + rep_vector(1, n_tot)./(1+exp(-a*(x-b))));
    return(w);
  }

  // VARIANCE MASK
  vector STANFUNC_var_mask(vector x, real a){
    int n_tot = num_elements(x);
    vector[n_tot] s = rep_vector(1, n_tot)./(1+exp(-a*x));
    return(s);
  }

  // COMPUTE X_TILDE
  vector STANFUNC_get_x_tilde(vector x_disAge, vector T_onset, vector T_observed, int[,] mapping, int[] lengths){
    int n_tot = num_elements(x_disAge);
    int N_cases = num_elements(lengths);
    vector[n_tot] x_tilde = rep_vector(0.0, n_tot);
    for(k in 1:N_cases){
      int inds[lengths[k]] = mapping[k, 1:lengths[k]];
      x_tilde[inds] = x_disAge[inds] + T_observed[k] - T_onset[k];
    }
    return(x_tilde);
  }

  // CATEGORICAL KERNEL 
  matrix STANFUNC_K_cat(vector x1, vector x2){
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

  // BINARY (MASK) KERNEL 
  matrix STANFUNC_K_bin_int(int[] x1, int[] x2, int c){
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

  // Binary (mask) kernel but with real inputs
  matrix STANFUNC_K_bin_real(vector x1, vector x2, int c){
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
  matrix STANFUNC_K_var_mask(vector x1_tilde, vector x2_tilde, real stp, real[] vm_params){
    int n1 = num_elements(x1_tilde);
    int n2 = num_elements(x2_tilde);
    real a = stp * vm_params[2];
    real h = vm_params[1];
    real r = 1/a*log(h/(1-h));
    vector[n1] s1 = STANFUNC_var_mask(x1_tilde - r, a);
    vector[n2] s2 = STANFUNC_var_mask(x2_tilde - r, a);
    matrix[n1,n2] S1 = rep_matrix(s1, n2);
    matrix[n1,n2] S2 = transpose(rep_matrix(s2, n1));
    matrix[n1,n2] MASK = S1 .* S2;
    return(MASK);
  }

  // Multiplier matrix to enable heterogeneous diseaseAge effect
  matrix STANFUNC_K_beta_symmetric(vector beta, int[] row_to_caseID){
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
  
    // Multiplier matrix to enable heterogeneous diseaseAge effect
  matrix STANFUNC_K_beta(vector beta, int[] row_to_caseID_1, int[] row_to_caseID_2){
    int n1 = num_elements(row_to_caseID_1);
    int n2 = num_elements(row_to_caseID_2);
    int i_caseID = 0;
    int j_caseID = 0;
    real tmp;
    matrix[n1,n2] BETA;
    for(i in 1:n1){
      i_caseID = row_to_caseID_1[i];
      for(j in 1:n2){
        j_caseID = row_to_caseID_2[j];
        if(i_caseID*j_caseID > 0){
          BETA[i,j] = sqrt(beta[i_caseID]*beta[j_caseID]);
        }else{
          BETA[i,j] = 0;
        }
      }
    }
    return(BETA);
  }
  
    
  // COMPUTE FIXED KERNEL MATRICES (do not depend on parameters)
  matrix[] STANFUNC_compute_fixed_kernel_matrices(vector[] X1, vector[] X2, int[] X1_nn, int[] X2_nn, int[] D, int cat_interact_kernel){
    int n1 = num_elements(X1[1]);
    int n2 = num_elements(X2[1]);
    int nf = 1 + D[3] + D[5] + D[6]; 
    matrix[n1,n2] KF[nf];
    KF[1] = STANFUNC_K_cat(X1[1], X2[1]);
    for(j in 1:D[3]){
      KF[1+j] = STANFUNC_K_bin_int(X1_nn, X2_nn, 1);
    }
    for(j in 1:D[5]){
      int ix = 2 + D[3] + D[4] + j;
      if(cat_interact_kernel == 1){
        KF[1+D[3]+j] = STANFUNC_K_cat(X1[ix], X2[ix]);
      }else{
        KF[1+D[3]+j] = STANFUNC_K_bin_real(X1[ix], X2[ix], 1);
      }
    }
    for(j in 1:D[6]){
      int ix = 2+D[3]+D[4]+D[5]+j;
      KF[1+D[3]+D[5]+j] = STANFUNC_K_cat(X1[ix], X2[ix]);
    }
    return(KF);
  }
  
  
  // COMPUTE ALL KERNEL MATRICES
  matrix[] STANFUNC_compute_kernel_matrices(vector[] X1, vector[] X2, int[,] caseID_to_rows_1, int[,] caseID_to_rows_2, int[] row_to_caseID_1, int[] row_to_caseID_2, int[] caseID_nrows_1, int[] caseID_nrows_2, matrix[] KF, vector[] T_onset, vector T_observed, int[] D, int UNCRT, int HMGNS, int USE_VAR_MASK, real[] vm_params, real[] alpha_idAge, real[] alpha_sharedAge, real[] alpha_diseaseAge, real[] alpha_continuous, real[] alpha_categAge, real[] alpha_categOffset, real[] lengthscale_idAge, real[] lengthscale_sharedAge, real[] lengthscale_diseaseAge, real[] lengthscale_continuous, real[] lengthscale_categAge, real[] warp_steepness, vector[] beta){
    int n1 = num_elements(X1[1]);
    int n2 = num_elements(X2[1]);
    real x1_age[n1] = to_array_1d(X1[2]);   // age covariate as an array
    real x2_age[n2] = to_array_1d(X2[2]);   // age covariate as an array
    int sum_D = sum(D);
    matrix[n1,n2] KX[sum_D];
    int r = 0;
    if(D[1]==1){
      real alp = alpha_idAge[1];
      real ell = lengthscale_idAge[1];
      r += 1;
      KX[r] = cov_exp_quad(x1_age, x2_age, alp, ell) .* KF[1];
    }
    if(D[2]==1){
      real alp = alpha_sharedAge[1];
      real ell = lengthscale_sharedAge[1];
      r += 1;
      KX[r] = cov_exp_quad(x1_age, x2_age, alp, ell);
    }
    if(D[3]==1){
      real alp = alpha_diseaseAge[1];
      real ell = lengthscale_diseaseAge[1];
      real stp = warp_steepness[1];
      vector[n1] x1_tilde;
      vector[n2] x2_tilde;
      real w1[n1];
      real w2[n2];
      r += 1;
    
      // Handle diseaseAge uncertainty
      if(UNCRT==0){
        x1_tilde = X1[3];
        x2_tilde = X2[3];
      }else{
        x1_tilde = STANFUNC_get_x_tilde(X1[3], T_onset[1], T_observed, caseID_to_rows_1, caseID_nrows_1);
        x2_tilde = STANFUNC_get_x_tilde(X2[3], T_onset[1], T_observed, caseID_to_rows_2, caseID_nrows_2); 
      }
      w1 = to_array_1d(STANFUNC_warp_input(x1_tilde, stp, 0.0, 1.0));
      w2 = to_array_1d(STANFUNC_warp_input(x2_tilde, stp, 0.0, 1.0));
    
      // Create disease effect kernel
      KX[r] = KF[2] .* cov_exp_quad(w1, w2, alp, ell);
      if(HMGNS==0){
        KX[r] = STANFUNC_K_beta(beta[1], row_to_caseID_1, row_to_caseID_2) .* KX[r];
      }
      if(USE_VAR_MASK==1){
        KX[r] = STANFUNC_K_var_mask(x1_tilde, x2_tilde, stp, vm_params) .* KX[r];
      }
    }
    for(j in 1:D[4]){
      real alp = alpha_continuous[j];
      real ell = lengthscale_continuous[j];
      int ix  = 2 + D[3] + j;
      r += 1;
      KX[r] = cov_exp_quad(to_array_1d(X1[ix]), to_array_1d(X2[ix]), alp, ell);
    }
    for(j in 1:D[5]){
      real alp = alpha_categAge[j];
      real ell = lengthscale_categAge[j];
      int ikf = 1 + D[3] + j;
      r += 1;
      KX[r] = cov_exp_quad(x1_age, x2_age, alp, ell) .* KF[ikf];
    }
    for(j in 1:D[6]){
      int ikf = 1 + D[3] + D[5] + j;
      real alp = alpha_categOffset[j];
      r += 1;
      KX[r] = square(alp) * KF[ikf];
    }
    return(KX);
  }

}
