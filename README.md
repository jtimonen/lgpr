# lgpr
[![coverage](https://codecov.io/gh/jtimonen/lgpr/branch/master/graph/badge.svg)](https://codecov.io/gh/jtimonen/lgpr)
[![travis](https://travis-ci.org/jtimonen/lgpr.svg?branch=master)](https://travis-ci.org/github/jtimonen/lgpr)

R-package for interpretable nonparametric modeling of longitudinal data using additive Gaussian processes. Contains functionality for inferring covariate effects and assessing covariate relevances. Various models can be specified using a convenient formula syntax.

## Getting started
See overview, tutorials and documentation at https://jtimonen.github.io/lgpr-usage/index.html. 

## Requirements

* The package should work on all major operating systems. 
* R 3.4 or later is required, R 4.0.2 or later is recommended
* We don't have pre-compiled binaries distributed yet, so currently *lgpr* needs to be installed from source. This is why you need to have your toolchain setup properly. See the instructions for:
  - [Linux](https://github.com/stan-dev/rstan/wiki/Configuring-C-Toolchain-for-Linux)
  - [Windows](https://github.com/stan-dev/rstan/wiki/Configuring-C---Toolchain-for-Windows)
  - [Mac](https://github.com/stan-dev/rstan/wiki/Configuring-C---Toolchain-for-Mac)

## License
[GPL>=3](https://www.gnu.org/licenses/gpl-3.0.html)

## Installation
* Install *lgpr* by running
```r
install.packages('devtools') # if you don't have devtools already
devtools::install_github('jtimonen/lgpr', dependencies = TRUE)
```
* If you have trouble installing the dependency [rstan](https://mc-stan.org/rstan/), see [these instructions](https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started)

* **Note:** In some cases, the exact version 2.0.0 of *rstantools* is required due
to problems with some recently updated versions of *StanHeaders* or *rstantools* (see [this thread](https://github.com/stan-dev/rstantools/issues/76)). You can remove possible incompatible version of *rstantools* and install version 2.0.0 by
```r
remove.packages('rstantools')
devtools::install_version("rstantools", version = "2.0.0", repos = "http://cran.us.r-project.org")
```

## Real data and reproducing the experiments
For code to reproduce the experiments of our manuscript see https://github.com/jtimonen/lgpr-usage. Preprocessed longitudinal proteomics
data is also provided there. See also the built-in `read_proteomics_data()` function.
