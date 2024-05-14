# SLSEdesign: Optimal designs using the second-order Least squares estimator 

*Chi-Kuang Yeh, Julie Zhou*

*May 13, 2024*

<!-- badges: start -->
[![R-CMD-check](https://github.com/chikuang/SLSEdesign/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/chikuang/SLSEdesign/actions/workflows/R-CMD-check.yaml)
[![Codecov test coverage](https://codecov.io/gh/chikuang/SLSEdesign/branch/master/graph/badge.svg)](https://app.codecov.io/gh/chikuang/SLSEdesign?branch=master)
<!-- badges: end -->

---

## Description

This is a package to compute the optimal regression design under the second-order Least squares estimator 

## Installation

SLSEdesign is available in Python, Julia and R. To install in R:

```r
devtools::install_github("chikuang/SLSEdesign")
```

## Examples

1. Calculate the D-optimal design for the Peleg model:

```r
peleg <- function(xi, theta){
  deno <- (theta[1] + xi * theta[2])^2
  rbind(-xi/deno, -xi^2/deno)
}
my_design <- Dopt(N = 31, u = seq(0, 180, length.out = 31), tt = 0, FUN = peleg,
                  theta = c(0.05, 0.5), num_iter = 500)
my_design$design
my_design$val
```

## TODO

+ [x] Add doi in the DESCRIPTION file
+ [x] Add files in inst/doc/ folder
+ [x] Add examples to all the exported functions to fulfill the requirement of CRAN
+ [x] Add information matrices in the examples within the calc_directional_derivatives functions
+ [ ] Upload the package to CRAN

## Reference 

1. Gao, Lucy L. and Zhou, Julie. (2017). [D-optimal designs based on the second-order least squares estimator](https://link.springer.com/article/10.1007/s00362-015-0688-9). *Statistical Papers*, 58, 77–94.
2. Gao, Lucy L. and Zhou, Julie. (2014). [New optimal design criteria for regression models with asymmetric errors](https://www.sciencedirect.com/science/article/pii/S037837581400007X). *Journal of Statistical Planning and Inference*, 149, 140-151.
3. Wang, Liqun and Leblanc, Alexandre. (2008). [Second-order nonlinear least squares estimation](https://link.springer.com/article/10.1007/s10463-007-0139-z). *Annals of the Institute of Statistical Mathematics*, 60, 883–900.
4. Yeh, Chi-Kuang and Zhou, Julie. (2021). [Properties of optimal regression designs under the second-order least squares estimator](https://link.springer.com/article/10.1007/s00362-018-01076-6). *Statistical Papers*, 62, 75–92.
5. Yin, Yue and Zhou, Julie. (2017). [Optimal designs for regression models using the second-order least squares estimator](https://www.jstor.org/stable/26384103). *Statistica Sinica*, 27, 1841-1856. 
