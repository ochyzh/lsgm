
# lsgm: Local Structure Graph Model Estimation

<!-- badges: start -->
<!-- badges: end -->

The goal of lsgm is to implement the methodologies described in Chyzh and Kaiser
(2019, DOI: 10.1017/pan.2019.8) and Nieman, Machain, Chyzh and Bell (2022, DOI:
10.1086/711716).

## Installation

You can install the development version of lsgm like so:

``` r
remotes::install_github("ochyzh/lsgm")
```

This package is not yet on CRAN.

## Example

This is a basic example which shows you how to estimate a local structure graph
model by using the package data `toy_data` that provides:

* Y: an 100x1 matrix that contains the values of the binary dependent variable.
* W: an 100x100 square matrix whose values represent connectivity between pairs
  of edges.
* X: an 100x1 matrix that contains 1 independent variable for 100 observations.

``` r
library(lsgm)
lsgm(toy_data$Y, toy_data$W, toy_data$X)
```

The documentation provides more examples and details on the functions and
arguments.
