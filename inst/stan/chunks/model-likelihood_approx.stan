if(is_f_sampled){
  
  reject("NOT IMPLEMENTED!");

}else{

  // delta?


    //matrix[n, R*M] V = STAN_bfa_create_V(num_obs, num_obs,
    //  K_const, components, x_cont, x_cont, x_cont_unnorm, x_cont_unnorm,
    //  alpha, ell, wrp, beta, teff,
    //  vm_params, idx_expand, idx_expand, teff_zero);
    
    
    //vector[R*M] D_diag = STAN_bfa_create_D(num_obs, num_obs,
    //  K_const, components, x_cont, x_cont, x_cont_unnorm, x_cont_unnorm,
    //  alpha, ell, wrp, beta, teff,
    //  vm_params, idx_expand, idx_expand, teff_zero);
      
  // Increase log prob
  //target += STAN_multi_normal_bf_lpdf(y, V, D_diag, sigma_n[1]);
  
  // Compute Phi and Lambda
  vector[num_basisfun] bfa_lambda[num_comps] = STAN_lambda_matrix(alpha, ell, 
      num_basisfun, width_basisfun);
  matrix[num_obs, num_basisfun] bfa_phi[num_comps] = STAN_phi_matrix(x_cont, 
      num_basisfun, width_basisfun);
  //print(PHI);
  print(bfa_lambda);

}
