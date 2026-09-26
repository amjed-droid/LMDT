\name{LMDT-package}
\alias{LMDT-package}
\alias{LMDT}
\docType{package}
\title{Local Manifold Dispersion Test for Heteroskedasticity}
\description{
Implements the Local Manifold Dispersion Test (LMDT) for detecting complex,
non-linear, and high-dimensional heteroskedasticity in linear regression models.
LMDT measures the graph-spectral Dirichlet energy of robust log-dispersion residuals
over a self-tuning k-nearest neighbors normalized graph Laplacian.
}
\details{
The primary user function is \code{\link{lmdt_test}}, which supports both
formula and matrix interfaces.
}
\author{
Ahmed Sattar Jabbar \email{ahmed.state.me@gmail.com}
}
\references{
Jabbar, A. S. (2026). The Local Manifold Dispersion Test: A Nonparametric
Graph-Spectral Approach for Detecting Complex and High-Dimensional
Heteroskedasticity. \emph{Working Paper}.
}
\keyword{package}
