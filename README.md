# lgpr
An R-package for **L**ongitudinal **G**aussian **P**rocess **R**egression.

* Bayesian nonparametric modeling of longitudinal data using additive Gaussian process regression
* Selection of categorical and continuous covariates, and interpretable visualization of covariate effects
* Disease effect modeling either homogeneously or heterogeneously across diagnosed patients
* Modeling uncertainty in the disease effect time
* Gaussian, Poisson and Negative Binomial observation models
* Efficient posterior inference using [Stan](https://mc-stan.org/)

## 1. Installation

### 1.1 Install RStan
To avoid problems, we recommend first installing the **rstan** package using these detailed instructions:
https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started

### 1.2 Install lgpr

A) On command line or Rstudio terminal:
* download or clone this repository
* put it in directory **lgpr** under parent directory **parent-dir**
* go to **parent-dir** on command line or Rstudio terminal
* run `R CMD build lgpr`
* this will build the package into **lgpr_X.Y.Z.tar.gz**
* install the package using `R CMD INSTALL lgpr_X.Y.Z.tar.gz`

OR

B) Using `devtools`:
* Install the `devtools` package, for example by running `install.packages('devtools')`
* Run `devtools::install_github('jtimonen/lgpr', build_vignettes = TRUE)`

## 2. Tutorials
See the demos at https://github.com/jtimonen/lgpr-usage, and the package vignettes.

## 3. Experiments
To reproduce the experiments of our manuscript, see ...
