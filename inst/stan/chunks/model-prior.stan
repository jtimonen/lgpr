// Magnitude parameters
for(j in 1:num_comps){
  target += STAN_log_prior(alpha[j], prior_alpha[j], hyper_alpha[j]);
}

// Lengthscale parameters
for(j in 1:num_ell){
  target += STAN_log_prior(ell[j], prior_ell[j], hyper_ell[j]);
}

// Input warping parameters
for(j in 1:num_ns){
  target += STAN_log_prior(wrp[j], prior_wrp[j], hyper_wrp[j]);
}

// Noise level parameters
if(obs_model==1){
  target += STAN_log_prior(sigma[1], prior_sigma[1], hyper_sigma[1]);
}else if(obs_model==3){
  target += STAN_log_prior(phi[1], prior_phi[1], hyper_phi[1]);
}

// Heterogeneity parameters
for(j in 1:num_heter){
  real a[2] = hyper_beta[j];
  target += beta_lpdf(beta[j] | a[1], a[2]);
}

// Disease-related age uncertainty
for(j in 1:num_uncrt){
  
  vector[num_cases] tx;
  int ptype = prior_teff[1][1];
  int is_backwards = prior_teff[1][2];
  int is_relative = prior_teff[1][3];
  
  if(is_relative){
    tx = teff[1] - teff_obs[1];
  }else if(is_backwards){
    tx = teff_obs[1] - teff[1];
  }else{
      tx = teff[1];
  }
  
  for(k in 1:num_cases){
    target += STAN_log_prior(tx[k], {ptype, 0}, hyper_teff[1]);
  }

}

// Isotropic normals for auxiliary variables when F is sampled
if(is_f_sampled){
  for(j in 1:num_comps){
    target += normal_lpdf(eta[1,j] | 0, 1);
  }
}
