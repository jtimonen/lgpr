/* Component-wise and total posterior means and variances. Let D be the
  number of components (<num_comps>). Then the elemens of <f_post> are
    - [1:D] = component means
    - [D+1] = total mean
    - [(D+2):(D+2+D)] = component variances
    - [2*(D+1)] = total variance
*/
vector[num_obs] f_post[is_generated_done, 2*(num_comps+1)];

if(is_generated_done){
  
  // Compute all kernel matrices
  matrix[num_obs, num_obs] KX[num_comps] = STAN_kernel_all(num_obs, num_obs,
      K_const, components, x_cont, x_cont, x_cont_unnorm, x_cont_unnorm,
      alpha, ell, wrp, beta, teff,
      vm_params, idx_expand, idx_expand, teff_obs);
      
  // Compute component-wise and total function posteriors
  f_post[1] = STAN_gp_posterior(KX, y_cont[1], delta, sigma[1]);
}
