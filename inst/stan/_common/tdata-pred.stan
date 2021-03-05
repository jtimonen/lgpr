  // Precompute fixed kernel matrices (pred vs. data)
  matrix[num_pred, num_obs] K_const_s[num_comps] = STAN_kernel_const_all(
    num_pred, num_obs, x_cat_PRED, x_cat, x_cont_mask_PRED, x_cont_mask, 
    x_cat_num_levels, components
  );
  
  // Precompute fixed kernel matrices (pred vs. pred)
  matrix[num_pred, num_pred] K_const_ss[num_comps] = STAN_kernel_const_all(
    num_pred, num_pred, x_cat_PRED, x_cat_PRED,
    x_cont_mask_PRED, x_cont_mask_PRED, x_cat_num_levels, components
  );
