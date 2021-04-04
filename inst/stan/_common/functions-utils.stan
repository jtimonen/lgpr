  // Compute total signal by summing components
  vector STAN_vectorsum(vector[] vecs, data int L){
    int num_vecs = size(vecs);
    vector[L] s = rep_vector(0, L);
    for (j in 1:num_vecs){
      s += vecs[j];
    }
    return(s);
  }
  
  // Sum an array of matrices
  matrix STAN_matrix_array_sum(matrix[] K){
    int n1 = rows(K[1]);
    int n2 = cols(K[1]);
    matrix[n1, n2] K_sum = K[1];
    for(j in 2:size(K)){
      K_sum += K[j];
    }
    return(K_sum);
  }
