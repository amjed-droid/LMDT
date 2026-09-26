# ==============================================================================
# Monte Carlo Simulation Suite: Reproducing Benchmark Simulation Tables
# "The Local Manifold Dispersion Test: A Nonparametric Graph-Spectral Approach"
# Author: Ahmed Sattar Jabbar (Mustansiriyah University)
# ==============================================================================

if (!requireNamespace("LMDT", quietly = TRUE)) {
  stop("Package 'LMDT' must be installed. Run devtools::install('.') or remotes::install_github('amjed-droid/LMDT').")
}
library(LMDT)

source("competitors.R")

generate_dgp <- function(n = 200, p = 2, dgp = "H0", error_dist = "normal") {
  beta <- rep(1, p)
  X <- matrix(stats::runif(n * p, min = -2, max = 2), nrow = n, ncol = p)
  sigma <- rep(1, n)
  
  if (dgp == "H0") {
    sigma <- rep(1, n)
  } else if (dgp == "DGP1_Linear") {
    sigma <- exp(0.6 * X[, 1])
  } else if (dgp == "DGP2_Ring") {
    radius <- sqrt(rowSums(X[, 1:min(2, p), drop = FALSE]^2))
    sigma <- 0.3 + 3.0 * exp(- (radius - 1.2)^2 / (2 * 0.3^2))
  } else if (dgp == "DGP3_Cluster") {
    cluster_assign <- sample(c(1, 2), size = n, replace = TRUE)
    X[cluster_assign == 1, 1:2] <- matrix(stats::rnorm(sum(cluster_assign == 1) * 2, mean = -1.2, sd = 0.5), ncol = 2)
    X[cluster_assign == 2, 1:2] <- matrix(stats::rnorm(sum(cluster_assign == 2) * 2, mean = 1.2, sd = 0.5), ncol = 2)
    sigma <- ifelse(cluster_assign == 1, 3.0, 0.4)
  } else if (dgp == "DGP4_HighDim_Manifold") {
    k_fac <- 3
    F_lat <- matrix(stats::rnorm(n * k_fac), nrow = n, ncol = k_fac)
    Loadings <- matrix(stats::rnorm(30 * k_fac), nrow = 30, ncol = k_fac)
    Noise_E <- matrix(stats::rnorm(n * 30, sd = 0.3), nrow = n, ncol = 30)
    X <- F_lat %*% t(Loadings) + Noise_E
    beta <- rep(0.5, 30)
    latent_r <- sqrt(F_lat[, 1]^2 + F_lat[, 2]^2)
    sigma <- 0.4 + 2.8 * exp(- (latent_r - 1.2)^2 / (2 * 0.35^2))
  } else if (dgp == "DGP5_Sparse_Noise") {
    p_sparse <- 50
    X <- matrix(stats::runif(n * p_sparse, min = -2, max = 2), nrow = n, ncol = p_sparse)
    beta <- rep(0.3, p_sparse)
    r <- sqrt(X[, 1]^2 + X[, 2]^2)
    sigma <- 0.3 + 3.0 * exp(- (r - 1.2)^2 / (2 * 0.3^2))
  }
  
  if (error_dist == "normal") {
    eps <- stats::rnorm(n)
  } else if (error_dist == "t3") {
    eps <- stats::rt(n, df = 3) / sqrt(3)
  }
  
  errors <- sigma * eps
  y <- as.numeric(X %*% beta) + errors
  list(X = X, y = y, residuals = errors)
}

run_simulation_benchmark <- function(n_reps = 100, n = 200, alpha = 0.05) {
  scenarios <- list(
    list(name = "H0 (Normal)", p = 2, dgp = "H0", dist = "normal", denoise = "none"),
    list(name = "H0 (Heavy-Tail t3)", p = 2, dgp = "H0", dist = "t3", denoise = "none"),
    list(name = "DGP 1: Linear", p = 2, dgp = "DGP1_Linear", dist = "normal", denoise = "none"),
    list(name = "DGP 2: Ring Manifold", p = 2, dgp = "DGP2_Ring", dist = "normal", denoise = "none"),
    list(name = "DGP 3: Clustered", p = 2, dgp = "DGP3_Cluster", dist = "normal", denoise = "none"),
    list(name = "DGP 4: Latent Manifold (p=30)", p = 30, dgp = "DGP4_HighDim_Manifold", dist = "normal", denoise = "pca"),
    list(name = "DGP 5: Sparse Noise (p=50)", p = 50, dgp = "DGP5_Sparse_Noise", dist = "normal", denoise = "metric")
  )
  
  results <- data.frame(
    Scenario = character(),
    p = integer(),
    LMDT_Asymp = numeric(),
    LMDT_WildBoot = numeric(),
    Breusch_Pagan = numeric(),
    White = character(),
    GQ = numeric(),
    stringsAsFactors = FALSE
  )
  
  cat("=============================================================================\n")
  cat(sprintf("Running Monte Carlo Benchmark (%d replications, n = %d, alpha = %.2f)\n", n_reps, n, alpha))
  cat("=============================================================================\n")
  
  set.seed(42)
  for (sc in scenarios) {
    cat(sprintf("\nProcessing %s ... ", sc$name))
    rej_lmdt_a <- 0
    rej_lmdt_b <- 0
    rej_bp     <- 0
    rej_white  <- 0
    white_nan  <- FALSE
    rej_gq     <- 0
    
    for (i in seq_len(n_reps)) {
      d <- generate_dgp(n = n, p = sc$p, dgp = sc$dgp, error_dist = sc$dist)
      fit <- stats::lm.fit(cbind(1, d$X), d$y)
      e_hat <- fit$residuals
      
      # LMDT Asymptotic
      res_a <- LMDT::lmdt_test(d$X, e_hat, method = "asymptotic", denoise = sc$denoise)
      if (res_a$p.value < alpha) rej_lmdt_a <- rej_lmdt_a + 1
      
      # LMDT Wild Bootstrap
      res_b <- LMDT::lmdt_test(d$X, e_hat, method = "wild_bootstrap", B = 199L, denoise = sc$denoise)
      if (res_b$p.value < alpha) rej_lmdt_b <- rej_lmdt_b + 1
      
      # Competitor: Breusch-Pagan
      bp_res <- bp_test(d$X, e_hat)
      if (!is.na(bp_res$p_value) && bp_res$p_value < alpha) rej_bp <- rej_bp + 1
      
      # Competitor: White
      wh_res <- white_test(d$X, e_hat)
      if (is.na(wh_res$p_value)) {
        white_nan <- TRUE
      } else if (wh_res$p_value < alpha) {
        rej_white <- rej_white + 1
      }
      
      # Competitor: Goldfeld-Quandt
      gq_res <- gq_test(d$X, e_hat, order_by = 1)
      if (!is.na(gq_res$p_value) && gq_res$p_value < alpha) rej_gq <- rej_gq + 1
    }
    
    white_display <- if (white_nan) "NaN (Collapses)" else sprintf("%.3f", rej_white / n_reps)
    
    row_data <- data.frame(
      Scenario = sc$name,
      p = sc$p,
      LMDT_Asymp = rej_lmdt_a / n_reps,
      LMDT_WildBoot = rej_lmdt_b / n_reps,
      Breusch_Pagan = rej_bp / n_reps,
      White = white_display,
      GQ = rej_gq / n_reps,
      stringsAsFactors = FALSE
    )
    results <- rbind(results, row_data)
    cat("Done.\n")
  }
  
  cat("\n=============================================================================\n")
  cat("Simulation Results (Replication of Table 3):\n")
  cat("=============================================================================\n")
  print(results)
  cat("=============================================================================\n")
  invisible(results)
}

if (!interactive()) {
  run_simulation_benchmark(n_reps = 100)
}
