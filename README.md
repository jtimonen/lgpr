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
Instructions:
* download or clone this repository
* put it in directory `lgpr-pkg` under parent directory `lgpr`
* go to `lgpr` and run `R CMD build lgpr-pkg`
* this will build the package into file `lgpr_X.Y.Z.tar.gz`
* install the package using `R CMD install lgpr_X.Y.Z.tar.gz`

## Basic usage
To be updated.
