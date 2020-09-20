# lgpr 1.0

## New features

* More general modeling options.
* Formula syntax where `|` is used for indicates interaction terms.
* Alternative advanced formula syntax with `gp()`, `gp_warp()`, `zerosum()` etc.
* Categorical covariates must now be specified as factors in data, and don't 
have to be numeric.
* Easier prior specification with `normal()`, `log_normal()`, `student_t()` etc.
* Beta binomial observation model.
* Prior and posterior predictive checks using`ppc()`, which interfaces to
  [bayesplot](http://mc-stan.org/bayesplot/).
* Covariate relevance assessment is now separated from model fitting.
functionality into the `relevances()` function.
* More general visualization options with `plot_data()`, `plot_fit()`, 
 `plot_fit_components()`
* Unified API for model fit and component plots using 
  `plot_api_g()` (grouped) and `plot_api_c()` (componentwise).

# lgpr 0.33

* First release.

# lgpr < 0.33

* Earlier development versions.
