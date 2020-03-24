if(VERBOSE==1){
  print(" ")
  print("* Number of data points = ", n);
  print("* Number of model components = ", sum_D);
  print("* Number of individuals = ", N_tot);
  print("* Components = ", D);
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
