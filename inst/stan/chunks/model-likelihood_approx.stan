if(is_f_sampled){
  
  reject("NOT IMPLEMENTED!");

}else{

  // delta?
    
    //vector[R*M] D_diag = STAN_bfa_create_D(num_obs, num_obs,
    //  K_const, components, x_cont, x_cont, x_cont_unnorm, x_cont_unnorm,
    //  alpha, ell, wrp, beta, teff,
    //  vm_params, idx_expand, idx_expand, teff_zero);
      
  // Increase log prob
  //target += STAN_multi_normal_bf_lpdf(y, V, D_diag, sigma_n[1]);
  
  // Compute Phi and Lambda
  vector[num_basisfun] bfa_lambda[num_comps] = STAN_lambda_matrix(alpha, ell, 
      num_basisfun, width_basisfun);
  matrix[num_obs, num_basisfun] bfa_phi[num_comps] = STAN_phi_matrix(
      x_cont, num_basisfun, width_basisfun, components);
  vector[RM] bfa_D = STAN_D_matrix(alpha, bfa_lambda, bfa_delta, ranks); // beta?
  matrix[num_obs, RM] bfa_V = STAN_V_matrix(bfa_phi, bfa_theta, RM);
      
  target += STAN_multi_normal_bfa_lpdf(y_cont[1]| bfa_V, bfa_D, sigma[1]);

}
