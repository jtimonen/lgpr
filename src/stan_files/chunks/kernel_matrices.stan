// Compute the kernel matrix for each component
  matrix[n,n] Kx = rep_matrix(0,n,n);
  matrix[n,n] KX[sum_D];
  matrix[n,n] Ky;
  matrix[n,n] Ly;
  int r = 0;
  if(D[1]==1){
    real alp = alpha_idAge[1];
    real ell = lengthscale_idAge[1];
    r += 1;
    KX[r] = cov_exp_quad(x_age, alp, ell) .* KF[1];
  }
  if(D[2]==1){
    real alp = alpha_sharedAge[1];
    real ell = lengthscale_sharedAge[1];
    r += 1;
    KX[r] = cov_exp_quad(x_age, alp, ell);
  }
  if(D[3]==1){
    real alp = alpha_diseaseAge[1];
    real ell = lengthscale_diseaseAge[1];
    real stp = warp_steepness[1];
    vector[n] x_tilde;
    real w[n];
    r += 1;
    
    // Handle diseaseAge uncertainty
    if(UNCRT==0){
      x_tilde = X[3];
    }else{
      x_tilde = get_x_tilde(X[3], T_onset[1], T_observed, caseID_to_rows, caseID_nrows); 
    }
    w = to_array_1d(sigmoid(x_tilde, stp, 0.0, 1.0));
    
    // Homogeneous or heterogeneous disease effect
    if(HMGNS==1){
      KX[r] = cov_exp_quad(w, alp, ell).* KF[2];
    }else{
      KX[r] = K_beta(beta[1], row_to_caseID).* cov_exp_quad(w, alp, ell).* KF[2];
    }
  }
  for(j in 1:D[4]){
    real alp = alpha_continuous[j];
    real ell = lengthscale_continuous[j];
    int ix  = 2 + D[3] + j;
    r += 1;
    KX[r] = cov_exp_quad(to_array_1d(X[ix]), alp, ell);
  }
  for(j in 1:D[5]){
    real alp = alpha_categAge[j];
    real ell = lengthscale_categAge[j];
    int ikf = 1 + D[3] + j;
    r += 1;
    KX[r] = cov_exp_quad(x_age, alp, ell) .* KF[ikf];
  }
  for(j in 1:D[6]){
    int ikf = 1 + D[3] + D[5] + j;
    real alp = alpha_categOffset[j];
    r += 1;
    KX[r] = square(alp) * KF[ikf];
  }

  // Sum
  for(j in 1:sum_D){
    Kx += KX[j]; 
  }

// Add jitter to ensure that Kx is pos. def.
Kx += diag_matrix(rep_vector(DELTA, n));

// Cholesky decomposition
Ky = Kx + diag_matrix(rep_vector(square(sigma_n[1]), n));
Ly = cholesky_decompose(Ky);


