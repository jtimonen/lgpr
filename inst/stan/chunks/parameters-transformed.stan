// Declare transformed parameters
vector[num_cases] teff[is_uncrt];
vector[num_obs] f_latent[is_f_sampled, num_comps];

// Transform raw effect times
if(is_uncrt){
  teff[1] = teff_lb[1] + (teff_ub[1] - teff_lb[1]) .* teff_raw[1];
}

// Transform isotropic latent function components
if(is_f_sampled){
  
  matrix[num_obs, num_obs] Delta = diag_matrix(rep_vector(delta, num_obs));
  matrix[num_obs, num_obs] KX[num_comps] = STAN_kernel_all(
      K_fixed, components, x_cont, x_cont, alpha, ell, wrp, beta, teff,
      vm_params, idx_expand, idx_expand, teff_obs);
      
  for(j in 1:num_comps){
    f_latent[1, j] = cholesky_decompose(KX[j] + Delta) * eta[1, j];
  }
  
}
