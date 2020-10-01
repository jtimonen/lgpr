// Categorical zero-sum kernel
matrix STAN_kernel_base_zerosum(
  data int[] x1,
  data int[] x2,
  data int num_cat)
{
  int n1 = size(x1);
  int n2 = size(x2);
  matrix[n1, n2] K;
  //if(num_cat <= 1){
  //  reject("STAN_kernel_base_zerosum: <num_cat> must be at least 2!")
  //}
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
matrix STAN_kernel_base_bin_mask(
  data int[] x1,
  data int[] x2)
{
  int n1 = size(x1);
  int n2 = size(x2);
  matrix[n1,n2] K;
  for(i in 1:n1){
    for(j in 1:n2){
      K[i,j] = (x1[i]==0)*(x2[j]==0);
    }
  }
  return(K);
}

/* 
  Compute one constant kernel matrix. Does not depend on parameters and
  therefore this function never needs to be evaluated during sampling.
*/
matrix STAN_kernel_const(
  data int[] x1,
  data int[] x2,
  data int kernel_type,
  data int num_cat)
{
  int n1 = num_elements(x1);
  int n2 = num_elements(x2);
  matrix[n1, n2] K;
  if(kernel_type==1){
    K = STAN_kernel_base_cat(x1, x2);
  }else if(kernel_type==2){
    K = STAN_kernel_base_bin_mask(x1, x2);
  }else{
    // kernel_type should be 0
    K = STAN_kernel_base_zerosum(x1, x2, num_cat);
  }
  return(K);
}

/* 
  Compute all constant kernel matrices. These do not depend on parameters and
  therefore this function never needs to be evaluated during sampling.
*/
matrix[] STAN_kernel_const_all(
  data int n1,
  data int n2,
  data int[,] x1,
  data int[,] x2,
  data int[,] x1_mask,
  data int[,] x2_mask,
  data int[] num_levels,
  data int[,] components)
{
  int num_comps = size(components);
  matrix[n1, n2] K_const[num_comps];

  for(j in 1:num_comps){
    matrix[n1, n2] K;
    int opts[9] = components[j];
    int ctype = opts[1];
    int ktype = opts[2];
    int idx_cat = opts[8];
    int idx_cont = opts[9];
    
    // Compute mask kernel for continuous covariate
    if(idx_cont != 0){
      K = STAN_kernel_const(x1_mask[idx_cont], x2_mask[idx_cont], 2, 0);
    }else{
      K = rep_matrix(1, n1, n2);
    }
    
    // Compute kernel for categorical covariate
    if(ctype == 0 || ctype == 2){
      int M = num_levels[idx_cat];
      K = K .* STAN_kernel_const(x1[idx_cat], x2[idx_cat], ktype, M);
    }
    K_const[j] = K;
  }
  return(K_const);
}
