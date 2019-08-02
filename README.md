# RcppCoDA

<!-- badges: start -->
[![Travis build status](https://travis-ci.org/jsilve24/RcppCoDA.svg?branch=master)](https://travis-ci.org/jsilve24/RcppCoDA)
[![Codecov test coverage](https://codecov.io/gh/jsilve24/RcppCoDA/branch/master/graph/badge.svg)](https://codecov.io/gh/jsilve24/RcppCoDA?branch=master)
<!-- badges: end -->

The goal of RcppCoDA is to create a blazing fast library for the analysis of 
compositional data with both C++ and R interfaces. 

Beyond standard compositional data transforms, RcppCoDA features the following:

* Support for compositional arrays (e.g., posterior samples of compositional matricies)
* Transformations of Covariace Matricies (and covariance arrays)
* Support for partially compositional arrays
* A library of functions for transforming between coordinate systems while staying in the simplex
* The Inter-Quartile Log Ratio Transform
* Proportionality Calculations (i.e., Phi Statistics)


## Installation

You can install the released version of RcppCoDA from GitHub with:

``` r
devtools::install_github("jsilve24/RcppCoDA")
```

## Example

I will fill this out more in the future. For the moment, here is an example
of how you can use RcppCoDA to transform a compositional array. 

``` r
> library(RcppCoDA)
>
> X = array(rnorm(20), dim=c(5, 2, 2))
> ilrInv(X)

, , 1

           [,1]       [,2]
[1,] 0.05249135 0.04709032
[2,] 0.06032787 0.31502648
[3,] 0.09483247 0.08592917
[4,] 0.23780862 0.04948804
[5,] 0.12890097 0.12073137
[6,] 0.42563871 0.38173463

, , 2

            [,1]       [,2]
[1,] 0.392597544 0.29713085
[2,] 0.312753433 0.19622643
[3,] 0.182176938 0.04474728
[4,] 0.046990696 0.33968326
[5,] 0.055951886 0.10022110
[6,] 0.009529502 0.02199108
```

