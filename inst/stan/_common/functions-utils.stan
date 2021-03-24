  // Compute total signal by summing components
  vector STAN_vectorsum(vector[] vecs, data int L){
    int num_vecs = size(vecs);
    vector[L] s = rep_vector(0, L);
    for (j in 1:num_vecs){
      s += vecs[j];
    }
    return(s);
  }
