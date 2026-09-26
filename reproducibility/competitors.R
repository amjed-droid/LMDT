# ==============================================================================
# Benchmark Competing Heteroskedasticity Tests for Simulation Studies
# Part of the Reproducibility Suite for:
# "The Local Manifold Dispersion Test: A Nonparametric Graph-Spectral Approach"
# ==============================================================================

#' Breusch-Pagan / Koenker Studentized Test
bp_test <- function(X, residuals, studentize = TRUE) {
  n <- length(residuals)
  e2 <- residuals^2
  mean_e2 <- mean(e2)
  Z <- cbind(1, as.matrix(X))
  fit <- stats::lm.fit(Z, e2)
  res_ss <- sum(fit$residuals^2)
  tot_ss <- sum((e2 - mean_e2)^2)
  r_sq <- if (tot_ss > 0) 1 - (res_ss / tot_ss) else 0
  df <- ncol(X)
  
  if (studentize) {
    v <- mean((e2 - mean_e2)^2)
    if (v < 1e-12) return(list(test = "Breusch-Pagan", statistic = 0, p_value = 1))
    stat <- n * r_sq * (mean_e2^2) / v
  } else {
    stat <- 0.5 * n * r_sq
  }
  
  p_val <- 1 - stats::pchisq(stat, df = df)
  list(test = "Breusch-Pagan", statistic = stat, df = df, p_value = p_val)
}

#' White's General Heteroskedasticity Test
white_test <- function(X, residuals) {
  n <- length(residuals)
  p <- ncol(X)
  e2 <- residuals^2
  
  terms_list <- list(as.matrix(X), as.matrix(X)^2)
  if (p > 1) {
    combos <- utils::combn(p, 2)
    cross_terms <- matrix(0, nrow = n, ncol = ncol(combos))
    for (c_idx in seq_len(ncol(combos))) {
      cross_terms[, c_idx] <- X[, combos[1, c_idx]] * X[, combos[2, c_idx]]
    }
    terms_list[[3]] <- cross_terms
  }
  
  W_design <- do.call(cbind, terms_list)
  if (ncol(W_design) >= n - 1) {
    return(list(test = "White", statistic = NA, df = ncol(W_design), p_value = NA, error = "DF_EXHAUSTED"))
  }
  
  qr_decomp <- qr(cbind(1, W_design))
  rank_W <- qr_decomp$rank - 1
  
  fit <- stats::lm.fit(cbind(1, W_design), e2)
  tot_ss <- sum((e2 - mean(e2))^2)
  r_sq <- if (tot_ss > 0) 1 - (sum(fit$residuals^2) / tot_ss) else 0
  stat <- n * r_sq
  df <- rank_W
  
  p_val <- 1 - stats::pchisq(stat, df = df)
  list(test = "White", statistic = stat, df = df, p_value = p_val)
}

#' Goldfeld-Quandt Test
gq_test <- function(X, residuals, order_by = 1, split_fraction = 0.2) {
  n <- length(residuals)
  p <- ncol(X)
  
  if (order_by == "norm") {
    ord_val <- sqrt(rowSums(X^2))
  } else {
    ord_val <- X[, order_by]
  }
  
  idx_order <- order(ord_val)
  n_omit <- round(n * split_fraction)
  n_sub <- floor((n - n_omit) / 2)
  
  idx_low <- idx_order[seq_len(n_sub)]
  idx_high <- idx_order[(n - n_sub + 1):n]
  
  X_low <- X[idx_low, , drop = FALSE]
  e_low <- residuals[idx_low]
  X_high <- X[idx_high, , drop = FALSE]
  e_high <- residuals[idx_high]
  
  fit_low <- stats::lm.fit(cbind(1, X_low), e_low)
  fit_high <- stats::lm.fit(cbind(1, X_high), e_high)
  
  s_low <- sum(fit_low$residuals^2)
  s_high <- sum(fit_high$residuals^2)
  
  df_num <- n_sub - p - 1
  if (df_num <= 1) {
    return(list(test = "Goldfeld-Quandt", statistic = NA, p_value = NA))
  }
  
  f_stat <- (s_high / df_num) / (s_low / df_num)
  p_val <- 1 - stats::pf(f_stat, df1 = df_num, df2 = df_num)
  
  list(test = "Goldfeld-Quandt", statistic = f_stat, df1 = df_num, df2 = df_num, p_value = p_val)
}

#' Distance Covariance (dCor) Test for Independence between X and |e|
dcor_test <- function(X, residuals, R = 199) {
  n <- length(residuals)
  X_mat <- as.matrix(X)
  Y_vec <- abs(residuals)
  
  # Pairwise distance matrices
  dist_X <- as.matrix(stats::dist(X_mat))
  dist_Y <- as.matrix(stats::dist(Y_vec))
  
  # U-centering
  A <- dist_X - rowMeans(dist_X) - matrix(colMeans(dist_X), n, n, byrow = TRUE) + mean(dist_X)
  B <- dist_Y - rowMeans(dist_Y) - matrix(colMeans(dist_Y), n, n, byrow = TRUE) + mean(dist_Y)
  
  dcov2 <- mean(A * B)
  dvarX <- mean(A * A)
  dvarY <- mean(B * B)
  
  if (dvarX * dvarY <= 0) return(list(test = "dCor", statistic = 0, p_value = 1))
  dcor_stat <- sqrt(max(0, dcov2) / sqrt(dvarX * dvarY))
  
  # Permutation test
  perm_stats <- numeric(R)
  for (r in seq_len(R)) {
    Y_perm <- sample(Y_vec)
    dist_Yp <- as.matrix(stats::dist(Y_perm))
    Bp <- dist_Yp - rowMeans(dist_Yp) - matrix(colMeans(dist_Yp), n, n, byrow = TRUE) + mean(dist_Yp)
    perm_stats[r] <- sqrt(max(0, mean(A * Bp)) / sqrt(dvarX * mean(Bp * Bp)))
  }
  
  p_val <- (1 + sum(perm_stats >= dcor_stat)) / (R + 1)
  list(test = "dCor", statistic = dcor_stat, p_value = p_val)
}
