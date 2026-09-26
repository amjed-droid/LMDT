\name{lmdt_graph}
\alias{lmdt_graph}
\title{Self-Tuning Manifold Graph Construction}
\usage{
lmdt_graph(X, k = NULL)
}
\arguments{
\item{X}{A numeric matrix or data frame of covariates (n x p).}

\item{k}{Number of nearest neighbors. Default is adaptive rule of thumb.}
}
\value{
A list containing matrices W, S, L, local scales sigma, and vertex degrees.
}
\description{
Constructs an undirected, self-tuning k-nearest neighbors affinity graph
and its normalized Laplacian matrix from a covariate matrix.
}
\examples{
set.seed(42)
X <- matrix(rnorm(50 * 2), ncol = 2)
g <- lmdt_graph(X, k = 5)
dim(g$S)
}
