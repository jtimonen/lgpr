vector[num_obs] f_post[is_generated_done, 2*(num_comps+1)];

if(is_generated_done){
  
  // Compute all kernel matrices
  matrix[num_obs, num_obs] KX[num_comps] = STAN_kernel_all(
      K_fixed, components, x_cont, x_cont, alpha, ell, wrp, beta, teff,
      vm_params, idx_expand, idx_expand, teff_obs);
      
  // Compute component-wise and total function posteriors
  f_post[1] = STAN_gp_posterior(KX, y_cont[1], delta, sigma[1]);
}
