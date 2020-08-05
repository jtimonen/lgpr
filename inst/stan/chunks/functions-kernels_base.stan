// Categorical zero-sum kernel
matrix STAN_kernel_base_zerosum(
  data int[] x1,
  data int[] x2,
  data int num_cat)
{
  int n1 = size(x1);
  int n2 = size(x2);
  matrix[n1, n2] K;
  if(num_cat <= 1){
    reject("STAN_kernel_base_zerosum: <num_cat> must be at least 2!")
  }
  for(i in 1:n1){
    for(j in 1:n2){
      if(x1[i]==x2[j]){
        K[i,j] = 1;
      }else{
       K[i,j] = -inv(num_cat-1); 
      }
    }
  }
  return(K);
}

// Categorical kernel
matrix STAN_kernel_base_cat(
  data int[] x1,
  data int[] x2)
{
  int n1 = size(x1);
  int n2 = size(x2);
  matrix[n1,n2] K;
  for(i in 1:n1){
    for(j in 1:n2){
      K[i,j] = (x1[i]==x2[j]);
    }
  }
  return(K);
}

// BINARY MASK KERNEL
matrix STAN_kernel_base_bin(
  data int[] x1,
  data int[] x2,
  data int c)
{
  int n1 = size(x1);
  int n2 = size(x2);
  matrix[n1,n2] K;
  for(i in 1:n1){
    for(j in 1:n2){
      K[i,j] = (x1[i]==c)*(x2[j]==c);
    }
  }
  return(K);
}

// DISEASE MASK KERNEL
matrix STAN_kernel_base_disease_mask(
  data int[] x1,
  data int[] x2)
{
  int n1 = size(x1);
  int n2 = size(x2);
  matrix[n1, n2] K;
  for(i in 1:n1){
    for(j in 1:n2){
      K[i,j] = (x1[i]>0)*(x2[j]>0);
    }
  }
  return(K);
}

// Multiplier matrix to enable variance masking
matrix STAN_kernel_base_var_mask(
  data vector x1,
  data vector x2,
  real steepness,
  data real[] vm_params)
{
  int n1 = num_elements(x1);
  int n2 = num_elements(x2);
  real a = steepness * vm_params[2];
  real r = inv(a)*logit(vm_params[1]);
  matrix[n1, 1] s1 = to_matrix(STAN_var_mask(x1 - r, a));
  matrix[n2, 1] s2 = to_matrix(STAN_var_mask(x2 - r, a));
  matrix[n1, n2] K = s1 * transpose(s2);
  STAN_check_real_positive(steepness);
  STAN_check_real_positive(vm_params[2]);
  STAN_check_prob_positive(vm_params[1]);
  return(K);
}
