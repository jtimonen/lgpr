// Declare transformed parameters
vector[num_bt] teff[num_uncrt>0];
vector[num_obs] f_latent[is_f_sampled, num_comps];

// Transform raw effect times
for(j in 1:num_uncrt){
  teff[j] = teff_lb[j] + (teff_ub[j] - teff_lb[j]) .* teff_raw[j];
}

// Transform isotropic latent function components
if(is_f_sampled){
  
  matrix[num_obs, num_obs] Delta = diag_matrix(delta_vec);
  matrix[num_obs, num_obs] KX[num_comps] = STAN_kernel_all(num_obs, num_obs,
      K_const, components, x_cont, x_cont, x_cont_unnorm, x_cont_unnorm,
      alpha, ell, wrp, beta, teff, 
      vm_params, idx_expand, idx_expand, teff_zero);
      
  for(j in 1:num_comps){
    f_latent[1, j] = cholesky_decompose(KX[j] + Delta) * eta[1, j];
  }
  
}
