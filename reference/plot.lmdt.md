# Plot Diagnostic Visualizations for LMDT

Produces diagnostic plots displaying the null distribution of the test
statistic (under bootstrap/permutation or asymptotic Gaussian) and the
ordered log-dispersion profile.

## Usage

``` r
# S3 method for class 'lmdt'
plot(x, ...)
```

## Arguments

- x:

  An object of class `"lmdt"`.

- ...:

  Additional graphical parameters passed to plotting functions.

## Value

Invisibly returns `x`.

## Examples

``` r
set.seed(42)
fit <- lm(rnorm(50) ~ matrix(rnorm(100), 50, 2))
res <- lmdt_test(fit, method = "asymptotic")
plot(res)
```
