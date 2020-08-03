// CATEGORICAL ZERO-SUM KERNEL
matrix STAN_kernel_base_zerosum(
  data int[] x1,
  data int[] x2,
  data int num_cat)
{
  int n1 = size(x1);
  int n2 = size(x2);
  matrix[n1,n2] K;
  if(num_cat <= 1){
    reject("num_cat must be at least 2!")
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

// CATEGORICAL KERNEL
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
  real c)
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

// Compute a discrete kernel matrix
// - Does not depend on parameters and therefore this function
//   never needs to be evaluated during sampling
matrix STAN_kernel_discrete(
  data int[] x1,
  data int[] x2,
  data int kernel_type,
  data int num_cat)
{
  int n1 = num_elements(x1);
  int n2 = num_elements(x2);
  matrix[n1, n2] K;
  if(kernel_type==0){
    K = STAN_kernel_base_zerosum(x1, x2, num_cat);
  }else if(kernel_type==1){
    K = STAN_kernel_base_cat(x1, x2);
  }else if(kernel_type==2){
    K = STAN_kernel_base_bin(x1, x2, 1);
  }else{
    reject("invalid kernel type");
  }
  return(K);
}
