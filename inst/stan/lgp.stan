/*
  This Stan file is part of the 'lgpr' package
  Author: Juho Timonen
    
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
#include chunks/functions-prior.stan
#include chunks/functions-posterior.stan
#include chunks/functions-kernels_const.stan
#include chunks/functions-kernels.stan
#include chunks/functions-basisfun.stan
}

data {
#include chunks/data.stan
}

transformed data{
#include chunks/data-transformed.stan
}

parameters {
#include chunks/parameters.stan
}

transformed parameters {
#include chunks/parameters-transformed.stan
}

model {
#include chunks/model-prior.stan
  if(is_likelihood_skipped){
  }else{
#include chunks/model-likelihood.stan
  }
}

generated quantities {
#include chunks/generated.stan
}
