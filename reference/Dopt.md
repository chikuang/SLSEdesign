# Calculate the D-optimal design under the SLSE

Calculate the D-optimal design under the SLSE

## Usage

``` r
Dopt(N, u, tt, FUN, theta, show_cvxr_status = FALSE)
```

## Arguments

- N:

  The number of sample points in the design space.

- u:

  The discretized design space.

- tt:

  The level of skewness. When tt=0, it is equivalent to compute the
  D-optimal design under the ordinary least squares estimator.

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

This function calculates the D-optimal design and the loss function
under the D-optimality. The loss function under D-optimality is defined
as the log determinant of the inverse of the Fisher information matrix.

## Examples

``` r
poly3 <- function(xi, theta){
  matrix(c(1, xi, xi^2, xi^3), ncol = 1)
}
Npt <- 101
my_design <- Dopt(N = Npt, u = seq(-1, +1, length.out = Npt),
   tt = 0, FUN = poly3, theta = rep(0,4))
round(my_design$design, 3)
#>     location weight
#> 1      -1.00  0.250
#> 28     -0.46  0.058
#> 29     -0.44  0.192
#> 73      0.44  0.192
#> 74      0.46  0.058
#> 101     1.00  0.250
my_design$val
#>           [,1]
#> [1,] -5.275115
```
