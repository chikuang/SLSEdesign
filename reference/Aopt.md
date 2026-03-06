# Calculate the A-optimal design under the second-order Least squares estimator

Calculate the A-optimal design under the second-order Least squares
estimator

## Usage

``` r
Aopt(N, u, tt, FUN, theta, show_cvxr_status = FALSE)
```

## Arguments

- N:

  The number of sample points in the design space.

- u:

  The discretized design space.

- tt:

  The level of skewness between 0 to 1 (inclusive). When tt=0, it is
  equivalent to compute the A-optimal design under the ordinary least
  squares estimator.

- FUN:

  The function to calculate the derivative of the given model.

- theta:

  The parameter value of the model.

- show_cvxr_status:

  A boolean variable to indicate whether to show the status of the CVXR
  optimization. By default, it is set to FALSE.

## Value

A list that contains 1. Value of the objective function at solution. 2.
Status. 3. Optimal design

## Details

This function calculates the A-optimal design and the loss function
under the A-optimality. The loss function under A-optimality is defined
as the trace of the inverse of the Fisher information matrix

## Examples

``` r
poly3 <- function(xi, theta){
  matrix(c(1, xi, xi^2, xi^3), ncol = 1)
}
Npt <- 101
my_design <- Aopt(N = Npt, u = seq(-1, +1, length.out = Npt),
   tt = 0, FUN = poly3, theta = rep(0,4))
round(my_design$design, 3)
#>     location weight
#> 1      -1.00   0.15
#> 28     -0.46   0.35
#> 74      0.46   0.35
#> 101     1.00   0.15
my_design$val
#> [1] 37.52534
```
