// Log prior density to be added to target
real STAN_log_prior(
  real x, 
  data int[] types,
  data real[] hyper)
{
  real lp;
  real a = hyper[1]; // prior hyperparameter 1
  real b = hyper[2]; // prior hyperparameter 2 (not used in student-t)
  real c = hyper[3]; // prior hyperparameter 3 (currently not used)
  real theta;
  int L1 = size(types);
  int L2 = size(hyper);
  if(L1!=2){ reject("size of <types> must be 2, found = ", L1); }
  if(L2!=3){ reject("size of <hyper> must be 3, found = ", L2); }

  // Possible transform and log of its absolute derivative
  if (types[2]==0){
    lp = 0;
    theta = x;
  }else if (types[2]==1){
    lp = log(fabs(2*x));
    theta = square(x);
  }else{
    reject("types[2] must be either 0 or 1, found = ", types[2]);
  }
  
  // Value of pdf
  if (types[1]==1){
    // do nothing
  }else if (types[1]==2){
    lp += normal_lpdf(theta|a,b); // 2 = normal
  }else if (types[1]==3){
    lp += student_t_lpdf(theta|a,0.0,1.0); // 3 = student-t
  }else if (types[1]==4){
    lp += gamma_lpdf(theta|a,b); // 4 = gamma
  }else if (types[1]==5){
    lp += inv_gamma_lpdf(theta|a,b); // 5 = inv-gamma
  }else if (types[1]==6){
    lp += lognormal_lpdf(theta|a,b); // 6 = log-normal
  }else{
    reject("types[1] must be an integer between 1 and 6, found = ", types[1]);
  }
  
  return(lp);
}
