# lgpr v0
This repository contains the source code for the R-package **lgpr**. This is an old version.

### Static release
v0.33.0: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3632542.svg)](https://doi.org/10.5281/zenodo.3632542)

### Notes
In this version, if you include a disease component in your model, you currently need to ensure that the `data` rows are sorted so that the values of `id_variable` are increasing and belong in the integer set *{1, ..., N}*, where *N* is the total number of individuals. Note that for control individuals, you need to indicate missing values of `disAge_variable` with `NaN`.

