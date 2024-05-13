# SLSEdesign: Optimal designs using the second-order Least squares estimator 

*Chi-Kuang Yeh, Julie Zhou*

*May 13, 2024*

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

## Reference 

1. [Properties of optimal regression designs under the second-order least squares estimator](https://link.springer.com/article/10.1007/s00362-018-01076-6)
