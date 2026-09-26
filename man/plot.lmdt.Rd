\name{plot.lmdt}
\alias{plot.lmdt}
\title{Plot Diagnostic Visualizations for LMDT}
\usage{
\method{plot}{lmdt}(x, ...)
}
\arguments{
\item{x}{An object of class \code{"lmdt"}.}

\item{...}{Additional graphical parameters passed to plotting functions.}
}
\value{
Invisibly returns \code{x}.
}
\description{
Produces diagnostic plots displaying the null distribution of the test statistic
(under bootstrap/permutation or asymptotic Gaussian) and the ordered log-dispersion profile.
}
\examples{
set.seed(42)
fit <- lm(rnorm(50) ~ matrix(rnorm(100), 50, 2))
res <- lmdt_test(fit, method = "asymptotic")
plot(res)
}
