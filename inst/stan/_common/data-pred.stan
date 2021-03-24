  // Covariates
  int num_pred;
  vector[num_pred] x_cont_PRED[num_cov_cont];
  vector[num_pred] x_cont_unnorm_PRED[num_cov_cont];
  int x_cont_mask_PRED[num_cov_cont, num_pred];
  int x_cat_PRED[num_cov_cat, num_pred];
  int<lower=1, upper=num_bt+1> idx_expand_PRED[num_pred];
