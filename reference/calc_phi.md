# Calculate the loss function of the A-, c- or D-optimal design

Calculate the loss function of the A-, c- or D-optimal design

## Usage

``` r
calc_phi(
  design,
  theta,
  FUN,
  tt,
  A,
  criterion = "D",
  cVec = rep(0, length(theta))
)
```

## Arguments

- design:

  The resulted design that contains the design points and the associated
  weights

- theta:

  The parameter value of the model

- FUN:

  The function to calculate the derivative of the given model.

- tt:

  The level of skewness

- A:

  The calculated covariance matrix

- criterion:

  The criterion to be used for the design, either "D" for D-optimality
  or "A" for A-optimality. Default is "D".

- cVec:

  c vector used to determine the combination of the parameters. This is
  only used in c-optimality

## Value

The loss of the model at each design points

## Details

This function calculates the loss function of the design problem under
the A- or D-optimality. The loss functions under A-, or D-optimality are
defined as the trace and log determinant of the inverse of the Fisher
information matrix

## Examples

``` r
my_design <- data.frame(location = c(0, 180), weight = c(1/2, 1/2))
theta <- c(0.05, 0.5)
peleg <- function(xi, theta){
   deno <- (theta[1] + xi * theta[2])^2
   rbind(-xi/deno, -xi^2/deno)
}
A <- matrix(c(1, 0, 0, 0, 0.2116, 1.3116, 0, 1.3116, 15.462521), byrow = TRUE, ncol = 3)
res <- calc_phi(my_design, theta, peleg, 0, A, criterion = "A")
res
#> [1]  0.0000 10.2395
```
