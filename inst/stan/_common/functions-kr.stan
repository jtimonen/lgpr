  // Kernel regression
  vector[] STAN_kr(
    vector[] f_latent,
    matrix[] KX,
    matrix[] KX_s,
    matrix[] KX_ss,
    real delta)
  {
    
    // Declare variables
    int num_comps = size(KX);
    int n = rows(KX[1]);
    int p = rows(KX_s[1]);
    int J = num_comps + 1;
    int inds[2];
    vector[p] F_KR[J];
    return(F_KR);
  
  }
