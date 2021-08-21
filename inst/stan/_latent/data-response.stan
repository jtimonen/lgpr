  int<lower=0> y_int[obs_model>1, N]; // response variable (int)
  vector[N] y[obs_model==1]; // response variable (real)
  vector[N] c_hat[obs_model>1];
  int<lower=1> y_num_trials[obs_model>3, N]; // for Binom or BB model
