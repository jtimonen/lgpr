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
#include chunks/functions_kernels_base.stan
#include chunks/functions_kernels.stan
#include chunks/functions_log_prior.stan
}


data {
#include chunks/data.stan
}


transformed data{
#include chunks/transformed_data.stan
}


parameters {
#include chunks/parameters.stan
}


transformed parameters {
#include chunks/transformed_parameters.stan
}


model {
#include chunks/model_priors.stan
  if(LH!=0){
#include chunks/model_likelihood.stan
  }
}


generated quantities {
#include chunks/generated_quantities.stan
}

