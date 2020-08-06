library(lgpr)
library(rstan)
stanmodel <- lgpr::get_stan_model()
rstan::expose_stan_functions(stanmodel)

