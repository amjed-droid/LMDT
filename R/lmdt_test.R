#' Local Manifold Dispersion Test for Heteroskedasticity
#'
#' Performs the Local Manifold Dispersion Test (LMDT) to detect complex, non-linear,
#' and manifold-structured heteroskedasticity in linear regression models.
#'
#' @param x Either a formula (e.g. \code{y ~ x1 + x2}) or a numeric design matrix of regressors.
#' @param ... Additional arguments passed to specific methods.
#'
#' @return An object of class \code{"lmdt"} and \code{"htest"} containing:
#'   \item{statistic}{The standardized LMDT test statistic (Z-score).}
#'   \item{p.value}{The p-value corresponding to the chosen inference calibration.}
#'   \item{tlmd}{The raw Rayleigh quotient statistic \eqn{T_{LMD}}.}
#'   \item{method_name}{A character string indicating the calibration method used.}
#'   \item{null.value}{The expected value under the null (1.000).}
#'   \item{alternative}{A character string describing the alternative hypothesis.}
#'   \item{data.name}{A character string giving the name(s) of the data.}
#'   \item{call}{The matched call.}
#'
#' @references
#' Jabbar, A. S. (2026). The Local Manifold Dispersion Test: A Nonparametric
#' Graph-Spectral Approach for Detecting Complex and High-Dimensional
#' Heteroskedasticity. \emph{Working Paper}.
#'
#' @export
lmdt_test <- function(x, ...) {
  UseMethod("lmdt_test")
}

#' @rdname lmdt_test
#' @method lmdt_test formula
#' @param formula A symbolic formula representing the regression model.
#' @param data An optional data frame or environment containing the model variables.
#' @export
lmdt_test.formula <- function(formula, data = list(), ...) {
  cl <- match.call()
  mf <- stats::model.frame(formula, data = data)
  y <- stats::model.response(mf)
  X <- stats::model.matrix(formula, data = data)
  
  # Remove intercept column from manifold graph construction
  intercept_col <- which(colnames(X) == "(Intercept)")
  if (length(intercept_col) > 0) {
    X_cov <- X[, -intercept_col, drop = FALSE]
  } else {
    X_cov <- X
  }
  
  fit <- stats::lm.fit(X, y)
  resids <- fit$residuals
  
  res <- lmdt_test.default(x = X_cov, y = resids, ...)
  res$data.name <- paste(deparse(formula))
  res$call <- cl
  res
}

#' @rdname lmdt_test
#' @method lmdt_test lm
#' @export
lmdt_test.lm <- function(x, ...) {
  cl <- match.call()
  X <- stats::model.matrix(x)
  intercept_col <- which(colnames(X) == "(Intercept)")
  if (length(intercept_col) > 0) {
    X_cov <- X[, -intercept_col, drop = FALSE]
  } else {
    X_cov <- X
  }
  resids <- stats::residuals(x)
  res <- lmdt_test.default(x = X_cov, y = resids, ...)
  res$data.name <- deparse(substitute(x))
  res$call <- cl
  res
}

#' @rdname lmdt_test
#' @method lmdt_test default
#' @param y A numeric vector of regression residuals (or response vector if x is regressors).
#' @param k An integer neighborhood size. If \code{NULL}, default rule of thumb is used.
#' @param method A character string specifying the p-value calibration method:
#'   \code{"wild_bootstrap"} (default), \code{"permutation"}, or \code{"asymptotic"}.
#' @param alternative A character string specifying the alternative hypothesis:
#'   \code{"one.sided"} (smooth manifold variance, default) or \code{"two.sided"} (general non-smooth).
#' @param denoise A character string specifying high-dimensional pre-filtering:
#'   \code{"none"} (default), \code{"pca"} (Gavish--Donoho thresholding), or
#'   \code{"metric"} (supervised metric learning).
#' @param B An integer specifying the number of bootstrap/permutation replications. Default is 499.
#' @param c A positive numeric stabilizer for the log-dispersion transform. Default is 1e-4.
#' @param lambda A numeric shrinkage prior for metric learning. Default is 0.15.
#' @export
lmdt_test.default <- function(x, y,
                              k = NULL,
                              method = c("wild_bootstrap", "permutation", "asymptotic"),
                              alternative = c("one.sided", "two.sided"),
                              denoise = c("none", "pca", "metric"),
                              B = 499L,
                              c = 1e-4,
                              lambda = 0.15,
                              ...) {
  cl <- match.call()
  X <- as.matrix(x)
  e <- as.numeric(y)
  n <- nrow(X)
  p <- ncol(X)
  
  if (length(e) != n) stop("Length of residuals y must match nrow(x).")
  method <- match.arg(method)
  alternative <- match.arg(alternative)
  denoise <- match.arg(denoise)
  
  # Robust log-dispersion transformation
  u <- log(abs(e) + c)
  
  # High-dimensional pre-filtering
  if (denoise == "pca") {
    X_eff <- lmdt_pca(X)
  } else if (denoise == "metric") {
    X_eff <- lmdt_metric(X, u, lambda = lambda)
  } else {
    X_eff <- X
  }
  
  # Graph Laplacian
  graph_obj <- lmdt_graph(X_eff, k = k)
  L <- graph_obj$L
  S <- graph_obj$S
  
  # Centered dispersion signal
  z_raw <- u - mean(u)
  denom <- sum(z_raw^2)
  if (denom <= .Machine$double.eps) denom <- 1e-10
  
  # Test statistic T_LMD = (z^T L z) / (z^T z)
  T_LMD_obs <- as.numeric((t(z_raw) %*% L %*% z_raw) / denom)
  
  # Asymptotic standard error sigma_0 = sqrt(2 * Tr(S^2)) / n
  tr_S2 <- sum(S^2)
  sigma_0 <- sqrt(2 * tr_S2) / n
  
  # Standardized test statistic Z_LMD = (1 - T_LMD) / sigma_0
  Z_stat <- (1 - T_LMD_obs) / sigma_0
  
  # Calibration
  if (method == "asymptotic") {
    if (alternative == "one.sided") {
      p_val <- 1 - stats::pnorm(Z_stat)
    } else {
      p_val <- 2 * (1 - stats::pnorm(abs(Z_stat)))
    }
  } else if (method == "permutation") {
    T_perm <- numeric(B)
    for (b in seq_len(B)) {
      u_perm <- sample(u)
      if (denoise == "metric") {
        X_b <- lmdt_metric(X, u_perm, lambda = lambda)
        g_b <- lmdt_graph(X_b, k = k)
        L_b <- g_b$L
      } else {
        L_b <- L
      }
      z_b <- u_perm - mean(u_perm)
      T_perm[b] <- as.numeric((t(z_b) %*% L_b %*% z_b) / sum(z_b^2))
    }
    if (alternative == "one.sided") {
      p_val <- (1 + sum(T_perm <= T_LMD_obs)) / (B + 1)
    } else {
      p_val <- (1 + sum(abs(T_perm - 1) >= abs(T_LMD_obs - 1))) / (B + 1)
    }
  } else if (method == "wild_bootstrap") {
    T_boot <- numeric(B)
    for (b in seq_len(B)) {
      v <- sample(c(-1, 1), size = n, replace = TRUE)
      e_star <- e * v
      u_star <- log(abs(e_star) + c)
      if (denoise == "metric") {
        X_b <- lmdt_metric(X, u_star, lambda = lambda)
        g_b <- lmdt_graph(X_b, k = k)
        L_b <- g_b$L
      } else {
        L_b <- L
      }
      z_b <- u_star - mean(u_star)
      T_boot[b] <- as.numeric((t(z_b) %*% L_b %*% z_b) / sum(z_b^2))
    }
    if (alternative == "one.sided") {
      p_val <- (1 + sum(T_boot <= T_LMD_obs)) / (B + 1)
    } else {
      p_val <- (1 + sum(abs(T_boot - 1) >= abs(T_LMD_obs - 1))) / (B + 1)
    }
  }
  
  names(Z_stat) <- "Z_LMD"
  names(T_LMD_obs) <- "T_LMD"
  null_val <- 1.0
  names(null_val) <- "Dirichlet energy ratio"
  
  alt_desc <- if (alternative == "one.sided") {
    "heteroskedasticity with local manifold smoothness (T_LMD < 1)"
  } else {
    "general non-smooth heteroskedasticity (|T_LMD - 1| > 0)"
  }
  
  res <- list(
    statistic = Z_stat,
    p.value = p_val,
    tlmd = T_LMD_obs,
    sigma0 = sigma_0,
    tr_s2 = tr_S2,
    method_name = paste0("Local Manifold Dispersion Test (LMDT, ", method, ")"),
    null.value = null_val,
    alternative = alt_desc,
    method = paste0("Local Manifold Dispersion Test (", method, ")"),
    data.name = deparse(cl$x),
    k = k,
    B = if (method != "asymptotic") B else NULL,
    denoise = denoise,
    u = u,
    residuals = e,
    resamples = if (method == "wild_bootstrap") T_boot else if (method == "permutation") T_perm else NULL,
    call = cl
  )
  
  class(res) <- c("lmdt", "htest")
  res
}
