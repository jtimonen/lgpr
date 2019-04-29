# lgpr
An R-package for **L**ongitudinal **G**aussian **P**rocess **R**egression and covariate selection.

* Nonparametric modeling of longitudinal data
* Selection of categorical and continuous covariates
* Disease effect modeling either homogeneously or heterogeneously across diagnosed patients
* Modeling uncertainty in the disease onset
* Efficient posterior inference using [Stan](https://mc-stan.org/)
* Visualization of inferred covariate effects

## Installation

### 1. Install RStan
To avoid problems, we recommend first installing the `rstan` package using these detailed instructions:
https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started

### 2. Install lgpr
The development version can be installed using the `devtools` package.
~~~r
require(devtools)
devtools::install_github('jtimonen/lgpr', build_vignettes = TRUE)
~~~

## Basic usage
See vignettes and demos. To be updated soon.
