\name{lmdt_metric}
\alias{lmdt_metric}
\title{Supervised Metric Learning for Sparse Noise Dilution}
\usage{
lmdt_metric(X, u, lambda = 0.15)
}
\arguments{
\item{X}{A numeric matrix of covariates (n x p).}

\item{u}{A numeric vector of log-dispersion scores.}

\item{lambda}{Shrinkage prior in (0, 1). Default is 0.15.}
}
\value{
A rescaled covariate matrix with filtered noise dimensions.
}
\description{
Computes correlation-based feature weights to soft-threshold uninformative
ambient noise coordinates.
}
\examples{
set.seed(42)
X <- matrix(rnorm(50 * 10), nrow = 50, ncol = 10)
u <- rnorm(50)
X_rescaled <- lmdt_metric(X, u, lambda = 0.15)
dim(X_rescaled)
}
