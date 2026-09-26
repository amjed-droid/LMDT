#' Self-Tuning Manifold Graph Construction
#'
#' Constructs an undirected, self-tuning k-nearest neighbors affinity graph
#' and its normalized Laplacian matrix from a covariate matrix.
#'
#' @param X A numeric matrix or data frame of covariates (n x p).
#' @param k An integer specifying the number of nearest neighbors. Default is
#'   \code{max(5, min(floor(2 * sqrt(nrow(X))), nrow(X) - 2))}.
#'
#' @return A list containing:
#'   \item{W}{The symmetric affinity matrix with zero diagonal.}
#'   \item{S}{The normalized adjacency matrix \eqn{S = D^{-1/2} W D^{-1/2}}.}
#'   \item{L}{The normalized graph Laplacian \eqn{L = I - S}.}
#'   \item{sigma}{A numeric vector of local bandwidth scales.}
#'   \item{degrees}{A numeric vector of graph vertex degrees.}
#'
#' @references
#' Zelnik-Manor, L., & Perona, P. (2004). Self-tuning spectral clustering.
#' \emph{Advances in Neural Information Processing Systems}, 17, 1601-1608.
#'
#' @export
#' @examples
#' set.seed(42)
#' X <- matrix(rnorm(100), ncol = 2)
#' g <- lmdt_graph(X, k = 5)
#' dim(g$S)
lmdt_graph <- function(X, k = NULL) {
  X <- as.matrix(X)
  n <- nrow(X)
  if (n < 6) stop("Sample size n must be at least 6.")
  
  if (is.null(k)) {
    k <- max(5, min(floor(2 * sqrt(n)), n - 2))
  } else {
    k <- as.integer(k)
    if (k < 2 || k >= n) stop("k must be between 2 and n - 1.")
  }
  
  # Pairwise squared Euclidean distances
  # Using cross-product decomposition: ||xi - xj||^2 = ||xi||^2 + ||xj||^2 - 2 <xi, xj>
  norms <- rowSums(X^2)
  D2 <- outer(norms, norms, "+") - 2 * (X %*% t(X))
  D2[D2 < 0] <- 0
  diag(D2) <- 0
  D_mat <- sqrt(D2)
  
  # Identify k-NN and local bandwidth sigma_i
  knn_idx <- matrix(0L, nrow = n, ncol = k)
  sigma <- numeric(n)
  
  for (i in seq_len(n)) {
    dists <- D_mat[i, ]
    dists[i] <- Inf
    ord <- order(dists)[seq_len(k)]
    knn_idx[i, ] <- ord
    sigma[i] <- dists[ord[k]]
  }
  
  sigma[sigma <= .Machine$double.eps] <- quantile(sigma[sigma > 0], 0.05, na.rm = TRUE)
  if (is.na(sigma[1]) || sigma[1] <= 0) sigma[] <- 1.0
  
  # Compute self-tuning affinity matrix W
  W <- matrix(0, nrow = n, ncol = n)
  for (i in seq_len(n)) {
    for (j in knn_idx[i, ]) {
      denom <- 2 * sigma[i] * sigma[j]
      if (denom > 0) {
        val <- exp(-D2[i, j] / denom)
        W[i, j] <- val
        W[j, i] <- val
      }
    }
  }
  diag(W) <- 0
  
  # Degree and normalized matrices
  deg <- rowSums(W)
  deg[deg <= .Machine$double.eps] <- 1e-6
  deg_inv_sqrt <- 1 / sqrt(deg)
  
  S <- deg_inv_sqrt * W * rep(deg_inv_sqrt, each = n)
  diag(S) <- 0
  L <- diag(n) - S
  
  list(W = W, S = S, L = L, sigma = sigma, degrees = deg)
}
