
int nf = 1 + D[3] + D[5] + D[6];     // number of fixed kernel matrices
int sum_D = sum(D);                  // total number of covariates
matrix[n,n] KF[nf] = STANFUNC_compute_fixed_kernel_matrices(X, X, X_notnan, X_notnan, D, cat_interact_kernel);
vector[n] mu = rep_vector(C_hat, n); // GP mean

if(verbose_mode==1){
  print(" ")
  print("* Observation model = ", LH);
  print("* Number of data points = ", n);
  print("* Number of model components = ", sum_D);
  print("* Number of individuals = ", N_tot);
  print("* Additional model info:")
  if(LH==2 || LH==3){
    print("  - C_hat = ", C_hat);
  }
  print("  - D = ", D);
  print("  - F_is_sampled = ", F_is_sampled)
  print("  - cat_interact_kernel = ", cat_interact_kernel)
  if(D[3]==1){
    print("* Disease modeling info: ");
    print("  - Number of cases = ", N_cases);
    print("  - UNCRT = ", UNCRT);
    print("  - HMGNS = ", HMGNS);
    print("  - USE_VAR_MASK = ", USE_VAR_MASK);
    if(USE_VAR_MASK==1){
      print("      o vm_params = ", vm_params);
    }
  }
  print(" ")
}
