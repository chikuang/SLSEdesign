# Calculate the c-optimal design under the SLSE with the given combination of the parameters

Calculate the c-optimal design under the SLSE with the given combination
of the parameters

## Usage

``` r
copt(N, u, tt, FUN, theta, cVec, show_cvxr_status = FALSE)
```

## Arguments

- N:

  The number of sample points in the design space.

- u:

  The discretized design space.

- tt:

  The level of skewness. When tt=0, it is equivalent to compute the
  c-optimal design under the ordinary least squares estimator.

- FUN:

  The function to calculate the derivative of the given model.

- theta:

  The parameter value of the model.

- cVec:

  c vector used to determine the combination of the parameters

- show_cvxr_status:

  A boolean variable to indicate whether to show the status of the CVXR
  optimization. By default, it is set to FALSE.

## Value

A list that contains 1. Value of the objective function at solution. 2.
Status. 3. Optimal design

## Details

This function calculates the c-optimal design and the loss function
under the c-optimality. The loss function under c-optimality is defined
as the log determinant of the inverse of the Fisher information matrix.

## Examples

``` r
poly3 <- function(xi, theta){
  matrix(c(1, xi, xi^2, xi^3), ncol = 1)
}
Npt <- 101
my_design <- copt(N = Npt, u = seq(-1, +1, length.out = Npt),
   tt = 0, FUN = poly3, theta = rep(0,4),
   cVec = c(0,1,1,1))
round(my_design$design, 3)
#>     location weight
#> 51         0    0.5
#> 101        1    0.5
my_design$val
#> [1] 4.000021
```
