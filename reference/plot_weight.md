# Plot the weight distribution of the optimal design for univaraite regression model

Plot the weight distribution of the optimal design for univaraite
regression model

## Usage

``` r
plot_weight(design)
```

## Arguments

- design:

  The resulted design that contains the design points and the associated
  weights

## Value

The plot that shows the given optimal design

## Details

This functions produce a figure that contains the location and their
associated weights of the resulted optimal design measures.

## Examples

``` r
Des = list(location = c(-1, +1), weight = c(0.5, 0.5))
plot_weight(Des)

```
