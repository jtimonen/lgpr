  int<lower=0> y_int[obs_model>1, N]; // response variable (int)
  real y_real[obs_model==1, N]; // response variable (real)
  int<lower=1> y_num_trials[obs_model>3, N]; // for Binom or BB model
  vector[N] c_hat;
