# lgpr 1.0

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
 [#12](https://github.com/jtimonen/lgpr/issues/12)).
 * Avoid adding `ggplot2::color_scale_manual()` if number of colors > 5 
 (fixes [#11](https://github.com/jtimonen/lgpr/issues/11)).
 
## 1.0.4
Edit type checking to work more generally on all systems (fixes [#5](https://github.com/jtimonen/lgpr/issues/5)).

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
* C++ versions of the stan model functions are
  now exposed to package namespace and also tested.

# lgpr 0.33

* First release.

# lgpr < 0.33

* Earlier development versions.
