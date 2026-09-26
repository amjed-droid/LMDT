#' Gavish--Donoho Optimal Adaptive PCA Denoising
#'
#' Filters high-dimensional correlated features by hard-thresholding singular values
#' at the Gavish--Donoho (2014) optimal threshold.
#'
#' @param X A numeric matrix or data frame of covariates (n x p).
#'
#' @return A matrix of retained principal component projections.
#'
#' @references
#' Gavish, M., & Donoho, D. L. (2014). The optimal hard threshold for singular
#' values is 4/sqrt(3). \emph{IEEE Transactions on Information Theory}, 60(8), 5040-5053.
#'
#' @export
#' @examples
#' set.seed(42)
#' X <- matrix(rnorm(200 * 30), nrow = 200, ncol = 30)
#' X_clean <- lmdt_pca(X)
#' ncol(X_clean)
lmdt_pca <- function(X) {
  X <- as.matrix(X)
  n <- nrow(X)
  p <- ncol(X)
  
  X_scaled <- scale(X, center = TRUE, scale = TRUE)
  X_scaled[is.na(X_scaled)] <- 0
  
  svd_decomp <- svd(X_scaled)
  s <- svd_decomp$d
  
  beta <- min(n, p) / max(n, p)
  omega <- 0.56 * beta^3 - 0.95 * beta^2 + 1.82 * beta + 1.43
  tau_gd <- omega * median(s)
  
  k_opt <- sum(s > tau_gd)
  k_opt <- max(2, min(k_opt, min(n, p) - 1))
  
  svd_decomp$u[, seq_len(k_opt), drop = FALSE] %*% diag(s[seq_len(k_opt)], nrow = k_opt)
}

#' Supervised Metric Learning for Sparse Noise Dilution
#'
#' Computes correlation-based feature weights to soft-threshold uninformative
#' ambient noise coordinates, preventing distance concentration in regimes where p > 50.
#'
#' @param X A numeric matrix or data frame of covariates (n x p).
#' @param u A numeric vector of robust log-dispersion scores.
#' @param lambda A numeric value in (0, 1) specifying the uniform shrinkage prior.
#'   Default is 0.15.
#'
#' @return A rescaled covariate matrix with filtered noise dimensions.
#'
#' @export
#' @examples
#' set.seed(42)
#' X <- matrix(rnorm(100 * 20), nrow = 100, ncol = 20)
#' u <- rnorm(100)
#' X_weighted <- lmdt_metric(X, u, lambda = 0.15)
#' dim(X_weighted)
lmdt_metric <- function(X, u, lambda = 0.15) {
  X <- as.matrix(X)
  n <- nrow(X)
  p <- ncol(X)
  
  tau_n <- 2 / sqrt(n)
  w <- numeric(p)
  
  for (j in seq_len(p)) {
    xj <- X[, j]
    cor_lin <- abs(cor(u, xj))
    cor_quad <- abs(cor(u, (xj - mean(xj))^2))
    rho_j <- max(cor_lin, cor_quad, na.rm = TRUE)
    if (is.na(rho_j)) rho_j <- 0
    w[j] <- (max(0, rho_j - tau_n))^2
  }
  
  sum_w <- sum(w)
  if (sum_w > 0) {
    w_norm <- w / sum_w
    w_final <- (1 - lambda) * w_norm + lambda * (1 / p)
  } else {
    w_final <- rep(1 / p, p)
  }
  
  X * rep(sqrt(w_final * p), each = n)
}
