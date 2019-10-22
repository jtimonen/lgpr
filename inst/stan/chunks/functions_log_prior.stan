// LOG PRIOR TO BE ADDED TO TARGET
real STANFUNC_log_prior(real x, int[] types, real[] hp){
  real lp;
  real a = hp[1]; // prior hyperparameter 1
  real b = hp[2]; // prior hyperparameter 2
  real c = hp[3]; // prior hyperparameter 3 (currently not used)
  real theta;
  
  // Possible transform and log of its absolute derivative
  if (types[2]==0){
    lp = 0;
    theta = x;
  }else if (types[2]==1){
    lp = log(fabs(2*x));
    theta = square(x);
  }else{
    reject("invalid value of types[2]!")
  }
  
  // Value of pdf
  if (types[1]==1){
    // do nothing
  }else if (types[1]==2){
    lp += normal_lpdf(theta|a,b);
  }else if (types[1]==3){
    lp += student_t_lpdf(theta|a,0,b);
  }else if (types[1]==4){
    lp += gamma_lpdf(theta|a,b);
  }else if (types[1]==5){
    lp += inv_gamma_lpdf(theta|a,b);
  }else if (types[1]==6){
    lp += lognormal_lpdf(theta|a,b);
  }else{
    reject("types[1] must be an integer between 1 and 6; found =", types[1]);
  }
  
  return(lp);
}
