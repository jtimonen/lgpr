
int nf = 1 + D[3] + D[5] + D[6];     // number of fixed kernel matrices
int sum_D = sum(D);                  // total number of covariates
matrix[n,n] KF[nf] = STAN_compute_fixed_kernel_matrices(X, X_notnan, D, N_tot, N_cat);
vector[n] mu = rep_vector(C_hat, n); // GP mean
int row_to_caseID_plus1[n] = rep_array(0, n);
for(i in 1:n){
  row_to_caseID_plus1[i] = row_to_caseID[i] + 1;
}

if(VERBOSE==1){
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
  print("  - F_IS_SAMPLED = ", F_IS_SAMPLED)
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
