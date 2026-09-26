\name{lmdt_pca}
\alias{lmdt_pca}
\title{Gavish--Donoho Optimal Adaptive PCA Denoising}
\usage{
lmdt_pca(X)
}
\arguments{
\item{X}{A numeric matrix of covariates (n x p).}
}
\value{
A matrix of retained principal components.
}
\description{
Filters high-dimensional correlated features by hard-thresholding singular values
at the Gavish--Donoho (2014) optimal threshold.
}
\examples{
set.seed(42)
X <- matrix(rnorm(100 * 20), nrow = 100, ncol = 20)
X_clean <- lmdt_pca(X)
ncol(X_clean)
}
