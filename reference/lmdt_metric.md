# Supervised Metric Learning for Sparse Noise Dilution

Computes correlation-based feature weights to soft-threshold
uninformative ambient noise coordinates.

## Usage

``` r
lmdt_metric(X, u, lambda = 0.15)
```

## Arguments

- X:

  A numeric matrix of covariates (n x p).

- u:

  A numeric vector of log-dispersion scores.

- lambda:

  Shrinkage prior in (0, 1). Default is 0.15.

## Value

A rescaled covariate matrix with filtered noise dimensions.

## Examples

``` r
set.seed(42)
X <- matrix(rnorm(50 * 10), nrow = 50, ncol = 10)
u <- rnorm(50)
X_rescaled <- lmdt_metric(X, u, lambda = 0.15)
dim(X_rescaled)
#> [1] 50 10
```
