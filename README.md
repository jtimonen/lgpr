# lgpr
[![coverage](https://codecov.io/gh/jtimonen/lgpr/branch/development/graph/badge.svg)](https://codecov.io/gh/jtimonen/lgpr)
![travis](https://travis-ci.org/jtimonen/lgpr.svg?branch=development)
[![zenodo](https://zenodo.org/badge/DOI/10.5281/zenodo.3632542.svg)](https://doi.org/10.5281/zenodo.3632542)

R package for **L**ongitudinal **G**aussian **P**rocess **R**egression. 

## Getting started
See overview, tutorials and documentation at https://jtimonen.github.io/lgpr-usage/index.html. 

## Requirements

* Linux is the preferred operating system. However, the package works also on
  Mac and Windows. We dont have pre-compiled binaries distributed yet, so 
  currently *lgpr* needs to be installed from source. 
* To compile the Stan code included in the package on Windows or Mac, you need
  to have your toolchain setup properly. 
  -  On Windows, install Rtools as explained [here](https://github.com/stan-dev/rstan/wiki/Installing-RStan-from-source-on-Windows#configuration). You also need to complete the **Configuration** step, as described in the
  above link.
  - On Mac, see toolchain configuration [here](https://github.com/stan-dev/rstan/wiki/Installing-RStan-from-source-on-a-Mac)
  - R 3.4 or later is required, R 4.0.2 or later is recommended

## Installation
* Install the dependency [rstan](https://mc-stan.org/rstan/) using [these instructions](https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started)
* Install *lgpr* by running
```r
install.packages('devtools') # if you don't have devtools already
devtools::install_github('jtimonen/lgpr', dependencies = TRUE)
```

* **Note:** In some cases, the exact version 2.0.0 of *rstantools* is required due
to problems with some recently updated versions of *StanHeaders* or *rstantools* (see [this thread](https://github.com/stan-dev/rstantools/issues/76)). You can remove possible incompatible version of *rstantools* and install version 2.0.0 by
```r
remove.packages('rstantools')
devtools::install_version("rstantools", version = "2.0.0", repos = "http://cran.us.r-project.org")
```
