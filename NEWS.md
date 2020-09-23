# lgpr 1.0

## New features

* More general modeling options.
* Formula syntax where `|` indicates interaction terms.
* Alternative advanced formula syntax with `gp()`, `gp_warp()`, `zerosum()` etc.
* Categorical covariates must now be specified as factors in data, and don't 
have to be numeric.
* Easier prior specification with `normal()`, `log_normal()`, `student_t()` etc.
* Beta binomial observation model.
* Prior and posterior predictive checks using `ppc()`, which interfaces to
  [bayesplot](http://mc-stan.org/bayesplot/).
* Component relevance assessment is now separated from model fitting into the `relevances()` function and selection into `select()`.
* Better prediction and plotting functionality with `get_pred()`, `pred()`, `plot_pred()`, and `plot_f()`.

## Automated testing
* Thorough unit tests using [test_that](https://testthat.r-lib.org/).
* C++ versions of the stan model functions are
  now exposed to package namespace and also tested.

# lgpr 0.33

* First release.

# lgpr < 0.33

* Earlier development versions.
