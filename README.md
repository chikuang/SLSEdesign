SLSEDesign: Optimal designs using the second-order Least squares
estimator
================
*Chi-Kuang Yeh, Julie Zhou*  

*May 21, 2024*

<!-- badges: start -->

[![CRAN
status](https://www.r-pkg.org/badges/version/SLSEdesign)](https://CRAN.R-project.org/package=SLSEdesign)
[![R-CMD-check](https://github.com/chikuang/SLSEdesign/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/chikuang/SLSEdesign/actions/workflows/R-CMD-check.yaml)
[![CRAN
download](http://cranlogs.r-pkg.org/badges/grand-total/SLSEdesign?color=blue)](https://cran.r-project.org/package=SLSEdesign)
[![](https://img.shields.io/github/languages/code-size/chikuang/SLSEdesign.svg)](https://github.com/chikuang/SLSEdesign)
<!-- badges: end -->

## Description

This is a package to compute the optimal regression design under the
second-order Least squares estimator

## Installation

SLSEdesign is now available on [CRAN](https://cran.r-project.org/).
Hence you may install it by typing

``` r
install.packages("SLSEdesign")
```

or you may download the develop version by typing

``` r
devtools::install_github("chikuang/SLSEdesign") # or pak::pkg_install("chikuang/SLSEdesign")
```

## Details

Consider a general regression model,
$$y_i=\eta(\mathbf{x}_i, \mathbf{\theta})+ \epsilon_i, \quad i=1, \ldots, n,$$
where $y_i$ is the $i$-th observation of a response variable $y$ at
design point $\mathbf{x}_i \in S \subset \mathbb{R}^p$, $S$ is a design
space, $\mathbf{\theta} \in \mathbb{R}^q$ is the unknown regression
parameter vector, response function
$\eta(\mathbf{x}_i, \mathbf{\theta})$ can be a linear or nonlinear
function of $\mathbf{\theta}$, and the errors $\epsilon_i$ are assumed
to be uncorrelated with mean zero and finite variance $\sigma^2$.

Let $\hat{\mathbf{\theta}}$ be an estimator of $\mathbf{\theta}$, such
as the least squares estimator. Various optimal designs are defined by
minimizing $\phi \left( \mathbb{c}ov(\hat{\mathbf{\theta}}) \right)$
over the design points $\mathbf{x}_1, \ldots, \mathbf{x}_n$, where
function $\phi(\cdot)$ can be determinant, trace, or other scalar
functions. The resulting designs are called optimal exact designs
(OEDs), which depend on the response function $\eta(\cdot,\cdot)$, the
design space $S$, the estimator $\hat{\mathbf{\theta}}$, the scalar
function $\phi(\cdot)$, and the number of points $n$.

Second order least-squares estimator is defined as

``` math
(\boldsymbol{\hat{\theta}}^\top,\hat{\sigma}^2)^\top:=\underset{\boldsymbol{\theta},\sigma^2}{\mathrm{argmin}}\sum_{i=1}^n \begin{pmatrix}
y_i-\eta(\boldsymbol{x}_i;\boldsymbol{\theta})\\
y_i^2-\eta^2(\boldsymbol{x}_i;\boldsymbol{\theta})-\sigma^2
\end{pmatrix}^\top W(\boldsymbol{x_i}) \begin{pmatrix}
y_i-\eta(\boldsymbol{x_i};\boldsymbol{\theta})\\
y_i^2-\eta^2(\boldsymbol{x_i};\boldsymbol{\theta})-\sigma^2
\end{pmatrix}.
```

### Comparison between ordinary least-squares and second order least-squares estimators

Note that $`W(\boldsymbol{x_i})`$ is a $`2\times 2`$ non-negative
semi-definite matrix which may or may not depend on $\boldsymbol{x_i}$
(Wang and Leblanc, 2008). It is clear that SLSE is a natural extension
of the OLSE which is defined based on the first-order difference
function
(i.e. $`y_i-\mathbb{E}[y_i]=y_i-\eta(\boldsymbol{x_i};\boldsymbol{\theta})`$).
On the other hand, SLSE is defined using not only the first-order
difference function, but also second-order difference function
(i.e. $`y_i^2-\mathbb{E}[y_i^2]=y_i^2-(\eta^2(\boldsymbol{x_i};\boldsymbol{\theta})+\sigma^2))`$.
One might think about the downsides of SLSE after discussing its
advantages over OLSE. SLSE does have its disadvantages. It is not a
linear estimator, and there is no closed-form solution. It requires more
computational resources compared to OLSE due to its nonlinearity.
However, numerical results can be easily computed for SLSE nowadays. As
a result, SLSE is a powerful alternative estimator to be considered in
research studies and real-life applications.

In particular, if we set the skewness parameter $t$ to be zero, the
resulting optimal designs under SLSE and OLSE **will be the same**!

## Examples

1.  Calculate the D-optimal design for the 3rd order polynomial
    regression model:

$$
y_i = \beta_0 + \beta_1 x_i + \beta_2 x_i^2 + \beta_3 x_i^3 +\varepsilon_i
$$

``` r
poly3 <- function(xi,theta){
    matrix(c(1, xi, xi^2, xi^3), ncol = 1)
}
my_design <- Dopt(N = 31, u = seq(-1, 1, length.out = 31), 
     tt = 0, FUN = poly3, theta = rep(1, 4), num_iter = 500)
my_design$design
my_design$val
```

Add equivalence theorem plot for D-optimal design:

``` r
design = data.frame(location = c(-1, -0.447, 0.447, 1),
 weight = rep(0.25, 4))
u = seq(-1, 1, length.out = 201)
plot_direction_Dopt(u, design, tt=0, FUN = poly3,
  theta = rep(0, 4))
```

<img src="man/fig/README-demo-equivalence.png" width="50%" />

## TODO

- [x] Add doi in the DESCRIPTION file
- [x] Add files in inst/doc/ folder
- [x] Add examples to all the exported functions to fulfill the
  requirement of CRAN
- [x] Add information matrices in the examples within the
  calc_directional_derivatives functions
- [x] Upload the package to CRAN
- [x] Update description
- [x] Generate README.md file using RMarkdown
- [ ] Python and Julia version of the package, which are expected to be
  faster than in R
- [ ] Merge the functions that compute the directional derivatives.
  Maybe adding an extra argument to indicate the design criterion used.

## Reference

1.  Gao, Lucy L. and Zhou, Julie. (2017). [D-optimal designs based on
    the second-order least squares
    estimator](https://link.springer.com/article/10.1007/s00362-015-0688-9).
    *Statistical Papers*, 58, 77–94.
2.  Gao, Lucy L. and Zhou, Julie. (2014). [New optimal design criteria
    for regression models with asymmetric
    errors](https://www.sciencedirect.com/science/article/pii/S037837581400007X).
    *Journal of Statistical Planning and Inference*, 149, 140-151.
3.  Wang, Liqun and Leblanc, Alexandre. (2008). [Second-order nonlinear
    least squares
    estimation](https://link.springer.com/article/10.1007/s10463-007-0139-z).
    *Annals of the Institute of Statistical Mathematics*, 60, 883–900.
4.  Yeh, Chi-Kuang and Zhou, Julie. (2021). [Properties of optimal
    regression designs under the second-order least squares
    estimator](https://link.springer.com/article/10.1007/s00362-018-01076-6).
    *Statistical Papers*, 62, 75–92.
5.  Yin, Yue and Zhou, Julie. (2017). [Optimal designs for regression
    models using the second-order least squares
    estimator](https://www.jstor.org/stable/26384103). *Statistica
    Sinica*, 27, 1841-1856.
