# Local Manifold Dispersion Test for Heteroskedasticity

Implements the Local Manifold Dispersion Test (LMDT) for detecting
complex, non-linear, and high-dimensional heteroskedasticity in linear
regression models. LMDT measures the graph-spectral Dirichlet energy of
robust log-dispersion residuals over a self-tuning k-nearest neighbors
normalized graph Laplacian.

## Details

The primary user function is
[`lmdt_test`](https://amjed-droid.github.io/LMDT/reference/lmdt_test.md),
which supports both formula and matrix interfaces.

## Author

Ahmed Sattar Jabbar <ahmed.state.me@gmail.com>

## References

Jabbar, A. S. (2026). The Local Manifold Dispersion Test: A
Nonparametric Graph-Spectral Approach for Detecting Complex and
High-Dimensional Heteroskedasticity. *Working Paper*.
