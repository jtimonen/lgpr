  // Log prior density to be added to target
  real STAN_log_prior(real x, data array[] int types, data array[] real p) {
    real log_prior = 0;
    real t = x;
  
    // Possible transform and log of its absolute derivative
    if (types[2]==1){
      log_prior += log(abs(2*x));
      t = square(x);
    }
    
    // Value of pdf
    if (types[1]==2){
      log_prior += normal_lpdf(t | p[1], p[2]); // 2 = normal
    }else if (types[1]==3){
      log_prior += student_t_lpdf(t | p[1], 0.0, 1.0); // 3 = student-t
    }else if (types[1]==4){
      log_prior += gamma_lpdf(t | p[1], p[2]); // 4 = gamma
    }else if (types[1]==5){
      log_prior += inv_gamma_lpdf(t | p[1], p[2]); // 5 = inv-gamma
    }else if (types[1]==6){
      log_prior += lognormal_lpdf(t | p[1], p[2]); // 6 = log-normal
    }
    
    return(log_prior);
  }
