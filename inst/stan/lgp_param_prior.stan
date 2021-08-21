#include _common/licence.stan

functions{
#include _common/functions-prior.stan
}

data {
#include _common/data-general.stan
#include _common/data-covariates.stan
#include _common/data-priors.stan
#include _latent/data-general.stan
}

parameters {
#include _common/params.stan
#include _latent/params.stan
}

transformed parameters {
#include _common/tparams.stan
}

model {
#include _common/priors.stan
#include _latent/priors.stan
}
