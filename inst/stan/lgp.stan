// *lgp.stan*
// This is the main Stan model file of the 'lgpr' package
// Author: Juho Timonen

/*
    lgpr is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    lgpr is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with lgpr.  If not, see <http://www.gnu.org/licenses/>.
*/

functions{
#include chunks/functions-utils.stan
#include chunks/functions-log_prior.stan
#include chunks/functions-kernels_discrete.stan
#include chunks/functions-kernels_continuous.stan
#include chunks/functions-kernels_single.stan
#include chunks/functions-kernels_many.stan
}

data {
#include chunks/data-data.stan
}

transformed data{
#include chunks/data-transformed.stan
}

parameters {
#include chunks/parameters-parameters.stan
}

transformed parameters {
#include chunks/parameters-transformed.stan
}

model {
// #include chunks/model-priors.stan
  if(is_likelihood_skipped){
  }else{
// #include chunks/model-likelihood.stan
  }
}

generated quantities {
// #include chunks/generated.stan
}
