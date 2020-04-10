# lgpr
This repository contains the source code for the R-package **lgpr**. 

### Getting started
See tutorials, installation instructions and documentation [here](https://jtimonen.github.io/lgpr-usage/index.html).

### Static release
v0.33.0: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3632542.svg)](https://doi.org/10.5281/zenodo.3632542)

### Notes
Disease effect modeling can currently be ensured to work correctly only if the `data` rows are sorted so that the values of `id_variable` are increasing and belong in the integer set *{1, ..., N}*, where *N* is the total number of individuals. Note that for control individuals, you need to indicate missing values of `disAge_variable` with `NaN`.

### Citation
Juho Timonen, Henrik Mannerström, Aki Vehtari and Harri Lähdesmäki (2019). *An interpretable probabilistic
 machine learning method for heterogeneous longitudinal studies*. https://arxiv.org/abs/1912.03549
