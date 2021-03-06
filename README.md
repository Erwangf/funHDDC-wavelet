
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Functional Model-Based Discrimination and Clustering using wavelets

<!-- badges: start -->

[![Travis build
status](https://travis-ci.org/Erwangf/funHDDC-wavelet.svg?branch=master)](https://travis-ci.org/Erwangf/funHDDC-wavelet)
[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://www.tidyverse.org/lifecycle/#experimental)
<!-- badges: end -->

The R package `funHDDCwavelet` allows to perform time series clustering,
via discrete wavelet transform, and modeling wavelet coefficients with a
parcimonious gaussian mixture model. Parameter estimation is done with
an EM algorithm.

## Installation

You can install the development version from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("Erwangf/funHDDC-wavelet")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(funHDDCwavelet)
#>    __             _   _ ____  ____   ____                         _      _   
#>   / _|_   _ _ __ | | | |  _ \|  _ \ / ___|_      ____ ___   _____| | ___| |_ 
#>  | |_| | | | '_ \| |_| | | | | | | | |   \ \ /\ / / _` \ \ / / _ \ |/ _ \ __|
#>  |  _| |_| | | | |  _  | |_| | |_| | |___ \ V  V / (_| |\ V /  __/ |  __/ |_ 
#>  |_|  \__,_|_| |_|_| |_|____/|____/ \____| \_/\_/ \__,_| \_/ \___|_|\___|\__|
#> 
#> 
#> 
## basic example code
```

# Acknowledgements

This research benefited from the support of the FMJH ’Program Gaspard
Monge for optimization and operations research and their interactions
with data science’, and from the support from EDF and Thales.
