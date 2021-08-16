  // Magnitude parameters
  for(j in 1:J){
    target += STAN_log_prior(alpha[j], prior_alpha[j], hyper_alpha[j]);
  }
  
  // Lengthscale parameters
  for(j in 1:num_ell){
    target += STAN_log_prior(ell[j], prior_ell[j], hyper_ell[j]);
  }
  
  // Input warping parameters
  for(j in 1:num_wrp){
    target += STAN_log_prior(wrp[j], prior_wrp[j], hyper_wrp[j]);
  }
  
  // Heterogeneity parameters
  if(num_het > 0) {
    target += beta_lpdf(beta[1] | hyper_beta[1][1], hyper_beta[1][2]);
  }
  
  // Disease-related age uncertainty
  if(idx_unc > 0) {
    int ptype = prior_teff[1][1];
    int is_backwards = prior_teff[1][2];
    real direction = (-1.0)^(is_backwards);
    vector[num_teff] tx = direction * (teff[1] - teff_zero[1]);
    for(k in 1:num_teff){
      target += STAN_log_prior(tx[k], {ptype, 0}, hyper_teff[1]);
    }
  }
