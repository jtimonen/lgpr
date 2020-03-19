
int nf = 1 + D[3] + D[5] + D[6];     // number of fixed kernel matrices
int sum_D = sum(D);                  // total number of covariates
matrix[n,n] KF[nf] = STAN_kernels_fixed(X, X_notnan, D, N_tot, N_cat);
int DO_GEN_QUANT = (1 - F_IS_SAMPLED) * (1 - SKIP_GQ);

if(VERBOSE==1){
  print(" ")
  print("* Number of data points = ", n);
  print("* Number of model components = ", sum_D);
  print("* Number of individuals = ", N_tot);
  print("* Additional model info:")
  print("  - LH = ", LH);
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
