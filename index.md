# LMDT: Local Manifold Dispersion Test for Heteroskedasticity

[![R-CMD-check](https://img.shields.io/badge/R--CMD--check-passing-brightgreen.svg)](https://github.com/amjed-droid/LMDT)
[![License: GPL
v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![CRAN
status](https://www.r-pkg.org/badges/version/LMDT)](https://CRAN.R-project.org/package=LMDT)

**LMDT** implements the **Local Manifold Dispersion Test** for detecting
complex, non-linear, and high-dimensional heteroskedasticity in linear
regression models.

Unlike classical tests (Breusch–Pagan, White) that rely on restrictive
linear specifications or suffer from degrees-of-freedom explosion when
$`p > 15`$, LMDT measures the **graph-spectral smoothness (Dirichlet
energy)** of robust log-dispersion residuals over a self-tuning
$`k`$-nearest neighbors normalized Laplacian matrix.

------------------------------------------------------------------------

## Key Features

- **Asymptotic Kurtosis Invariance:** The zero-diagonal structure of the
  affinity operator asymptotically eliminates fourth-order error moments
  from the limiting variance at rate $`O(1/n)`$, providing robust size
  calibration under heavy-tailed (e.g., Student-$`t`$) error
  distributions.
- **Fast Computational Complexity:** Evaluates in
  $`O(k \cdot n \log n)`$ time, scaling easily to thousands of
  observations.
- **High-Dimensional Scaling ($`p > 50`$):** Combines unsupervised
  Gavish–Donoho optimal singular-value thresholding and supervised
  metric learning to prevent noise dilution.
- **S3 Object-Oriented Interface:** Works seamlessly with
  [`lm()`](https://rdrr.io/r/stats/lm.html) fit objects, symbolic
  formulas (`y ~ x1 + x2`), and numeric matrices.
- **Diagnostic Visualizations:** Native S3
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) method
  displays the null distribution and log-dispersion profiles.

------------------------------------------------------------------------

## Installation

### From GitHub (Development Version)

``` r

# install.packages("remotes")
remotes::install_github("amjed-droid/LMDT")
```

### From CRAN (Upcoming)

``` r

install.packages("LMDT")
```

------------------------------------------------------------------------

## Quick Start

### 1. Basic Usage with `lm()`

``` r

library(LMDT)

# Simulate regression data with non-linear heteroskedasticity
set.seed(42)
n <- 200
x1 <- runif(n, -2, 2)
x2 <- runif(n, -2, 2)
r <- sqrt(x1^2 + x2^2)
sigma <- 0.3 + 3.0 * exp(-(r - 1.2)^2 / (2 * 0.3^2))
y <- 2 * x1 - x2 + rnorm(n, sd = sigma)

# Fit linear model
fit <- lm(y ~ x1 + x2)

# Perform LMDT test
res <- lmdt_test(fit)
print(res)

# Plot diagnostic visualizations
plot(res)
```

### 2. High-Dimensional Regression ($`p = 60`$)

``` r

# Regressors with ambient noise dimensions
X <- matrix(rnorm(n * 60), n, 60)
y <- X[, 1] + rnorm(n, sd = 1 + abs(X[, 1]))

# Use supervised metric learning or adaptive PCA denoising
res_hd <- lmdt_test(X, residuals(lm(y ~ X)), denoise = "metric")
summary(res_hd)
```

------------------------------------------------------------------------

## Calibration Methods

- `method = "wild_bootstrap"` (Default): Recommended for finite samples
  and non-Gaussian errors.
- `method = "permutation"`: Exact non-parametric permutation calibration
  under $`H_0`$.
- `method = "asymptotic"`: Fast Gaussian approximation via de Jong’s
  Central Limit Theorem for degenerate quadratic forms.

------------------------------------------------------------------------

## Citation

To cite `LMDT` in publications, please use:

``` bibtex
@article{jabbar2026lmdt,
  title={The Local Manifold Dispersion Test: A Nonparametric Graph-Spectral Approach for Detecting Complex and High-Dimensional Heteroskedasticity},
  author={Jabbar, Ahmed Sattar},
  journal={Working Paper, Department of Statistics, Mustansiriyah University},
  year={2026}
}
```

------------------------------------------------------------------------

## Author & Maintainer

**Ahmed Sattar Jabbar**  
Department of Statistics, College of Administration and Economics,  
Mustansiriyah University, Baghdad, Iraq.  
Email: <ahmed.state.me@gmail.com>

------------------------------------------------------------------------

## License

GPL (\>= 3) © Ahmed Sattar Jabbar.
