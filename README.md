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
To avoid problems, we recommend first installing the *rstan* package using these detailed instructions:
https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started

### 2. Install lgpr
Instructions:
* download or clone this repository
* put it in directory _lgpr-pkg_ under parent directory _lgpr_
* go to _lgpr_ and run `R CMD build lgpr-pkg`
* this will build the package into file _lgpr_X.Y.Z.tar.gz_
* install the package using `R CMD INSTALL lgpr_X.Y.Z.tar.gz`

## Basic usage
To be updated.
