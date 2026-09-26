\name{lmdt_test}
\alias{lmdt_test}
\alias{lmdt_test.default}
\alias{lmdt_test.formula}
\alias{lmdt_test.lm}
\title{Local Manifold Dispersion Test for Heteroskedasticity}
\usage{
lmdt_test(x, ...)

\method{lmdt_test}{formula}(formula, data = list(), ...)

\method{lmdt_test}{lm}(x, ...)

\method{lmdt_test}{default}(
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
}
\arguments{
\item{x}{A numeric design matrix of regressors, or a formula.}

\item{...}{Further arguments passed to methods.}

\item{formula}{A symbolic regression formula (e.g. \code{y ~ x1 + x2}).}

\item{data}{An optional data frame containing model variables.}

\item{y}{A numeric vector of regression residuals or response variable.}

\item{k}{An integer number of nearest neighbors. Default uses adaptive rule of thumb.}

\item{method}{P-value calibration: \code{"wild_bootstrap"} (default), \code{"permutation"}, or \code{"asymptotic"}.}

\item{alternative}{Alternative hypothesis: \code{"one.sided"} (smooth manifold variance, default) or \code{"two.sided"}.}

\item{denoise}{High-dimensional pre-filtering: \code{"none"} (default), \code{"pca"}, or \code{"metric"}.}

\item{B}{Number of bootstrap or permutation replications. Default is 499.}

\item{c}{Positive stabilizer for log-dispersion transform. Default is 1e-4.}

\item{lambda}{Shrinkage prior for metric learning. Default is 0.15.}
}
\value{
An object of class \code{"lmdt"} and \code{"htest"}.
}
\description{
Performs the Local Manifold Dispersion Test (LMDT) for regression heteroskedasticity.
}
\examples{
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

\donttest{
# Using Wild Bootstrap calibration
test_boot <- lmdt_test(y ~ x1 + x2, method = "wild_bootstrap", B = 199)
print(test_boot)
}
}
\references{
Jabbar, A. S. (2026). The Local Manifold Dispersion Test: A Nonparametric
Graph-Spectral Approach for Detecting Complex and High-Dimensional
Heteroskedasticity. \emph{Working Paper}.
}
