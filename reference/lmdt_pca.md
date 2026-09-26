# Gavish–Donoho Optimal Adaptive PCA Denoising

Filters high-dimensional correlated features by hard-thresholding
singular values at the Gavish–Donoho (2014) optimal threshold.

## Usage

``` r
lmdt_pca(X)
```

## Arguments

- X:

  A numeric matrix of covariates (n x p).

## Value

A matrix of retained principal components.

## Examples

``` r
set.seed(42)
X <- matrix(rnorm(100 * 20), nrow = 100, ncol = 20)
X_clean <- lmdt_pca(X)
ncol(X_clean)
#> [1] 2
```
