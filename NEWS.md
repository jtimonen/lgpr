# lgpr 1.0

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
