# lgpr

[![travis](https://app.travis-ci.com/jtimonen/lgpr.svg?branch=master)](https://app.travis-ci.com/github/jtimonen/lgpr)
[![coverage](https://codecov.io/gh/jtimonen/lgpr/branch/master/graph/badge.svg)](https://app.codecov.io/gh/jtimonen/lgpr)
[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/lgpr)](https://cran.r-project.org/package=lgpr)
[![metacran downloads](https://cranlogs.r-pkg.org/badges/grand-total/lgpr)](https://cran.r-project.org/package=lgpr)
[![license](https://img.shields.io/badge/license-GPL>=3-lightgrey.svg)](https://www.gnu.org/licenses/gpl-3.0.html)

R-package for interpretable nonparametric modeling of longitudinal data using additive Gaussian processes. Contains functionality for inferring covariate effects and assessing covariate relevances. Various models can be specified using a convenient formula syntax.

## Getting started
See overview, tutorials, vignettes and documentation at https://jtimonen.github.io/lgpr-usage/index.html. 

## Requirements
* The package should work on all major operating systems. 
* R 3.4 or later is required, R 4.0.2 or later is recommended

## Installing from CRAN
* The latest released version that is available from CRAN can be installed simply via
```r
install.packages("lgpr")
```
Installing from CRAN is probably the easiest option since they might have binaries for your system (so no need to build the package from source yourself).

## Installing from source
* The latest released version (which might not be in CRAN yet) can be installed via
```r
install.packages('devtools') # if you don't have devtools already
devtools::install_github('jtimonen/lgpr', build_vignettes = TRUE)
```
* The latest development version can be installed via
```r
devtools::install_github('jtimonen/lgpr', ref = "develop")
``` 
Github installations are source installations (they require a C++ compiler).

* If you have trouble installing the dependency [rstan](https://mc-stan.org/rstan/), see [these instructions](https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started)
* Installing from source requires that you have your toolchain setup properly.
See the instructions for:
  - [Linux](https://github.com/stan-dev/rstan/wiki/Configuring-C-Toolchain-for-Linux)
  - [Windows](https://github.com/stan-dev/rstan/wiki/Configuring-C---Toolchain-for-Windows)
  - [Mac](https://github.com/stan-dev/rstan/wiki/Configuring-C---Toolchain-for-Mac)

## Real data and reproducing the experiments
For code to reproduce the experiments of our manuscript see https://github.com/jtimonen/lgpr-usage. Preprocessed longitudinal proteomics
data is also provided there. See also the built-in `read_proteomics_data()` function.
