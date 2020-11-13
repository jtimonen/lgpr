if(is_f_sampled){
  
  reject("NOT IMPLEMENTED!");

}else{
  print("num_basisfun = ", num_basisfun);
  print("width_basisfun = ", width_basisfun);
  // F NOT SAMPLED
  //vector[num_obs] sigma2_vec = rep_vector(square(sigma[1]), num_obs);
  //matrix[num_obs, num_obs] Ky = diag_matrix(num_comps * delta_vec);
  //matrix[num_obs, num_obs] KX[num_comps] = STAN_kernel_all(num_obs, num_obs,
  //    K_const, components, x_cont, x_cont, x_cont_unnorm, x_cont_unnorm,
  //    alpha, ell, wrp, beta, teff,
  //    vm_params, idx_expand, idx_expand, teff_zero);

  //for(j in 1:num_comps){
  //  Ky += KX[j];
  //}
  //Ky = Ky + diag_matrix(sigma2_vec);
  //y_cont[1] ~ multi_normal_cholesky(c_hat, cholesky_decompose(Ky));
}
