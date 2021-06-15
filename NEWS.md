# lgpr 1.1

## 1.1.1-1.1.2
  * Only small patches in documentation etc. in order to conform to CRAN policies
  
## 1.1.0
  * Adds `prior_pred()` for prior predictive sampling and `sample_param_prior()` for sampling from the parameter prior.
  * Adds `read_proteomics_data()` function.
  * Relax data type checking, to require that they only inherit from factor
  or numeric. Allow also `tibble`s and `data.table`s to be passed as data.
  * Adds more methods for `lgpfit` and `lgpmodel` objects, see their
  documentation.
  * Lot of improvements internally. Kernel computations
  in functions like `pred()` should take a lot less memory now. Two separete
  main Stan models now. One for latent GP (signal where f is sampled) and other for GP with marginalized f.
  * Improved documentation.
  * Improve verbose messages to user.
  
# lgpr 1.0

## 1.0.13
  * Fix bug that ignored the `group_by` argument in `get_teff_obs()` and
  caused at least `new_x()` to not work if the subject identifier variable
  was called something else than `id` (see [issue #22](https://github.com/jtimonen/lgpr/issues/22)).
  
## 1.0.12
  * Add more informative error message if trying to specify a model like
  `y ~ age + id | age`, which should be `y ~ age + age | id`, i.e. the
  continuous covariate on the left of `|` and categorical on the right.
  * New startup message that prints also `rstan` version
  * Update citation information
  
## 1.0.11
  * Add the `c_hat_pred` argument to `pred()`, to be used when `f` has been
  sampled and `c_hat` is not constant. Previously, `c_hat = 0` was used in
  all prediction points, which did not make sense in all cases.
  
## 1.0.10
  * Allow setting `group_by = NA` in `plot_pred()`,
  `plot_components()` and `new_x()` to avoid grouping in plots.
  * Allow setting `color_by` as the same factor as `group_by`.
  * Fix bug which caused an error when trying to define a separate prior
  for parameters of the same type.
    
## 1.0.9
  * Internal change for more effective computation of function (component)
    posterior variances.
  
## 1.0.8
  * Add option `do_yrng` which controls whether to do draws from the
  predictive distribution. This was previously always done if `sample_f`
  was `TRUE`. That is now considered a bug because it is unnecessary work if
  the `y_rng` draws are not needed. So the default is now `do_yrng = FALSE`,
  since `do_yrng = TRUE` can cause errors with the negative binomial model due
  to numerical problems (see [here](https://discourse.mc-stan.org/t/numerical-stability-of-gps-with-negative-binomial-likelihood/19343/5)). These problems should be addressed in a future
  release to allow more stable prior and posterior predictive sampling.
  
## 1.0.7
  * Small documentation update.
  
## 1.0.6
 * Fix bug in `get_pred()`, which was caused by not adding the GP mean to
 the sampled signal. This was causing postprocessing functions like
 `relevances()` and `plot_pred()` to give
 erroneous results if the GP mean was not a vector of zeros and 
 `sample_f = TRUE`.
 * Small edits in documentation and verbose information messages.
 
## 1.0.5
 * Make `plot_pred()` work with any response variable name (fixes 
 [issue #12](https://github.com/jtimonen/lgpr/issues/12)).
 * Avoid adding `ggplot2::color_scale_manual()` if number of colors > 5 
 (fixes [issue #11](https://github.com/jtimonen/lgpr/issues/11)).
 
## 1.0.4
Edit type checking to work more generally on all systems (fixes [issue #5](https://github.com/jtimonen/lgpr/issues/5)).

## 1.0.3
Fix CITATION to point to new preprint.

## 1.0.2
Added RcppParallel dependency explicitly.

## 1.0.1
Added warning if using default prior for input warping steepness.

## 1.0.0

### New features

* More general modeling options, allowing more mixing of different
  types of kernels/options
* Prior and posterior predictive checks using `ppc()`, which interfaces to
  [bayesplot](http://mc-stan.org/bayesplot/).

### Changes and improvements
* Formula syntax where `|` indicates interaction terms.
* Alternative advanced formula syntax with `gp()`, `gp_warp()`, `zerosum()` etc.
* Beta binomial observation model.
* Categorical covariates must now be specified as factors in data, and don't 
have to be numeric.
* Component relevance assessment is now separated from model fitting into the `relevances()` function and selection into `select()`.
* Easier prior specification with `normal()`, `log_normal()`, `student_t()` etc.
* Better prediction and plotting functionality with `get_pred()`, `pred()`, `plot_pred()`, and `plot_f()`.
* Extensive argument checking (see `check_positive_all()` etc.) to give
  users informative error messages

### Automated testing
* Thorough unit tests using [test_that](https://testthat.r-lib.org/).
* C++ versions of the Stan model functions are
  now exposed to package namespace and also tested.

# lgpr 0.33

* First release.

# lgpr < 0.33

* Earlier development versions.
