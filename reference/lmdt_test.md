# Local Manifold Dispersion Test for Heteroskedasticity

Performs the Local Manifold Dispersion Test (LMDT) for regression
heteroskedasticity.

## Usage

``` r
lmdt_test(x, ...)

# S3 method for class 'formula'
lmdt_test(formula, data = list(), ...)

# S3 method for class 'lm'
lmdt_test(x, ...)

# Default S3 method
lmdt_test(
  x,
  y,
  k = NULL,
  method = c("wild_bootstrap", "permutation", "asymptotic"),
  alternative = c("one.sided", "two.sided"),
  denoise = c("none", "pca", "metric"),
  B = 499L,
  c = 1e-04,
  lambda = 0.15,
  ...
)
```

## Arguments

- x:

  A numeric design matrix of regressors, or a formula.

- ...:

  Further arguments passed to methods.

- formula:

  A symbolic regression formula (e.g. `y ~ x1 + x2`).

- data:

  An optional data frame containing model variables.

- y:

  A numeric vector of regression residuals or response variable.

- k:

  An integer number of nearest neighbors. Default uses adaptive rule of
  thumb.

- method:

  P-value calibration: `"wild_bootstrap"` (default), `"permutation"`, or
  `"asymptotic"`.

- alternative:

  Alternative hypothesis: `"one.sided"` (smooth manifold variance,
  default) or `"two.sided"`.

- denoise:

  High-dimensional pre-filtering: `"none"` (default), `"pca"`, or
  `"metric"`.

- B:

  Number of bootstrap or permutation replications. Default is 499.

- c:

  Positive stabilizer for log-dispersion transform. Default is 1e-4.

- lambda:

  Shrinkage prior for metric learning. Default is 0.15.

## Value

An object of class `"lmdt"` and `"htest"`.

## References

Jabbar, A. S. (2026). The Local Manifold Dispersion Test: A
Nonparametric Graph-Spectral Approach for Detecting Complex and
High-Dimensional Heteroskedasticity. *Working Paper*.

## Examples

``` r
set.seed(42)
n <- 100
x1 <- rnorm(n)
x2 <- rnorm(n)
e <- rnorm(n, sd = exp(0.5 * x1))
y <- 1 + 2 * x1 + 3 * x2 + e
fit <- lm(y ~ x1 + x2)

# Using formula interface with asymptotic calibration
test_asymp <- lmdt_test(y ~ x1 + x2, method = "asymptotic")
print(test_asymp)
#> 
#>   Local Manifold Dispersion Test (LMDT, asymptotic)
#> 
#> data:  y ~ x1 + x2
#> Z_LMD = 7.8466, T_LMD = 0.7716, p-value = 2.109e-15
#> alternative hypothesis: heteroskedasticity with local manifold smoothness (T_LMD < 1)
#> denoising scheme: none
#> 

# \donttest{
# Using Wild Bootstrap calibration
test_boot <- lmdt_test(y ~ x1 + x2, method = "wild_bootstrap", B = 199)
print(test_boot)
#> 
#>   Local Manifold Dispersion Test (LMDT, wild_bootstrap)
#> 
#> data:  y ~ x1 + x2
#> Z_LMD = 7.8466, T_LMD = 0.7716, p-value = 1
#> alternative hypothesis: heteroskedasticity with local manifold smoothness (T_LMD < 1)
#> replications: B = 199, neighborhood size: k = 
#> denoising scheme: none
#> 
# }
```
