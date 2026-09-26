# Self-Tuning Manifold Graph Construction

Constructs an undirected, self-tuning k-nearest neighbors affinity graph
and its normalized Laplacian matrix from a covariate matrix.

## Usage

``` r
lmdt_graph(X, k = NULL)
```

## Arguments

- X:

  A numeric matrix or data frame of covariates (n x p).

- k:

  Number of nearest neighbors. Default is adaptive rule of thumb.

## Value

A list containing matrices W, S, L, local scales sigma, and vertex
degrees.

## Examples

``` r
set.seed(42)
X <- matrix(rnorm(50 * 2), ncol = 2)
g <- lmdt_graph(X, k = 5)
dim(g$S)
#> [1] 50 50
```
