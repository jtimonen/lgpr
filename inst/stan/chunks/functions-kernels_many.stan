/* 
  Compute all fixed kernel matrices. These do not depend on parameters and
  therefore this function never needs to be evaluated during sampling.
*/
matrix[] STAN_kernel_fixed_all(
  data int n1,
  data int n2,
  data int[,] x1, 
  data int[,] x2,
  data int[] num_levels,
  data int[,] components)
{
  int L1 = size(x1);
  int L2 = size(x2);
  int L3 = size(num_levels);
  int L4 = size(components);
  int num_comps = size(components[1]);
  matrix[n1, n2] K_fixed[num_comps];
  if(L1!=L2){
    reject("first dims of <x1> and <x2> must match! ",
           "found = (", L1, ", ", L2, ")");
  }
  if(L1!=L3){
    reject("size of <num_levels> must match first dim of <x1> and <x2>! ",
           "found = ", L3, ", should be = ", L1);
  }
  if(L4!=4){
    reject("first dimension of <components> must be 4!");
  }

  for(j in 1:num_comps){
    int ctype = components[1][j];
    if(ctype != 1){
      int ktype = components[2][j];
      int idx = components[3][j];
      int x1_j[n1] = x1[idx];
      int x2_j[n2] = x2[idx];
      int n_levels = num_levels[idx];
      K_fixed[j] = STAN_kernel_fixed(x1_j, x2_j, n_levels, ctype, ktype); 
    }else{
      K_fixed[j] = rep_matrix(0, n1, n2);
    }
  }
  return(K_fixed);
}

/* 
  Compute all kernel matrices. These depend on parameters and
  therefore this function needs to be evaluated repeatedly during sampling.
*/
matrix[] STAN_kernel_all(
  data int n1,
  data int n2,
  data matrix[] K_fixed,
  data int[,] components,
  data vector[] x1,
  data vector[] x2,
  real[] alpha,
  real[] ell,
  real[] wrp,
  vector[] beta,
  vector[] teff,
  data real[,] vm_params,
  data int[,] idx1_expand,
  data int[,] idx2_expand,
  data vector[] teff_obs)
{
  int ell_idx = 0;
  int num_comps = size(components[1]);
  matrix[n1, n2] KX[num_comps];
  int L1 = size(x1);
  int L2 = size(x2);
  int L3 = size(K_fixed);
  int L4 = size(components);

  if(L1 != L2){
    reject("first dims of <x1> and <x2> must match! ",
           "found = (", L1, ", ", L2, ")");
  }
  if(num_comps != L3){
    reject("second dim of <components> and first dim of <K_fixed> must match! ",
           "found = (", num_comps, ", ", L3, ")");
  }
  if(L4 != 4){
    reject("first dimension of <components> must be 4!");
  }
  
  for(j in 1:num_comps){
    
    // Component type + covariate
    int ctype = components[1][j];
    int idx = components[4][j];
    
    // Compute the kernel matrix
    if(ctype==0){
      KX[j] = square(alpha[j]) * K_fixed[j];
    }else{
      ell_idx += 1;
      if(ctype!=3){
        KX[j] = STAN_kernel_stationary(
          K_fixed[j], x1[idx], x2[idx], ctype, alpha[j], ell[ell_idx]);
      }else{
        KX[j] = STAN_kernel_disease(
          K_fixed[j], x1[idx], x2[idx], alpha[j], ell[ell_idx], wrp[1],
          beta, teff, vm_params, idx1_expand, idx2_expand, teff_obs);
      }
    }
  }
  
  return(KX);
}
