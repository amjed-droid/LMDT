# ==============================================================================
# Generate Publication-Quality Figures for LMDT Research
# Reproducing Figure 1: Geometric Behavior of LMDT vs Breusch-Pagan
# Author: Ahmed Sattar Jabbar (Mustansiriyah University)
# ==============================================================================

set.seed(42)
n <- 300

# 1. Generate Ring Manifold Data (DGP 2)
X <- matrix(stats::runif(n * 2, -2, 2), nrow = n, ncol = 2)
radius <- sqrt(X[, 1]^2 + X[, 2]^2)
sigma <- 0.3 + 3.0 * exp(-(radius - 1.2)^2 / (2 * 0.3^2))
eps <- stats::rnorm(n, 0, 1)
y <- X[, 1] + X[, 2] + sigma * eps
fit <- lm(y ~ X[, 1] + X[, 2])
resids <- residuals(fit)
abs_res <- abs(resids)
u <- log(abs_res + 1e-4)

# 2. Output High-Resolution PNG
png("figure1_manifold_visualization.png", width = 3600, height = 1100, res = 300)

par(mfrow = c(1, 3), mar = c(4.5, 4.5, 3.5, 1.5), oma = c(0, 0, 1, 0))

# Subplot (A): Spatial distribution of |e_i| over 2D feature manifold
palette_cols <- hcl.colors(100, "Viridis")
res_scaled <- (abs_res - min(abs_res)) / (max(abs_res) - min(abs_res))
col_idx <- ceiling(res_scaled * 99) + 1
col_points <- palette_cols[col_idx]

plot(X[, 1], X[, 2], col = col_points, pch = 19, cex = 1.3,
     main = "(A) Residual Dispersion |e_i| on Manifold",
     xlab = expression(X[1]), ylab = expression(X[2]),
     cex.main = 1.2, font.main = 2)
grid()

# Subplot (B): Breusch-Pagan Linear Projection (Flat line, completely blind)
plot(X[, 1], resids^2, col = "firebrick3", pch = 19, cex = 1.1,
     main = "(B) Breusch-Pagan Linear Trend (Blind)",
     xlab = expression(X[1]), ylab = expression(e[i]^2),
     cex.main = 1.2, font.main = 2)
abline(h = mean(resids^2), col = "royalblue3", lwd = 2.5, lty = 2)
legend("topright", legend = "BP Flat OLS Fit", col = "royalblue3", lwd = 2.5, lty = 2, bty = "n")
grid()

# Subplot (C): LMDT Spectral Signal along Manifold Radius
ord <- order(radius)
plot(radius[ord], u[ord], col = "purple3", pch = 19, cex = 1.0,
     main = "(C) LMDT Spectral Signal along Manifold Radius",
     xlab = expression("Manifold Radius" ~ "||" * X[i] * "||"[2]),
     ylab = expression("Log-Dispersion" ~ u[i] == ln(abs(e[i]) + c)),
     cex.main = 1.2, font.main = 2)

# Smoothed trend
smooth_trend <- stats::lowess(radius[ord], u[ord], f = 0.2)
lines(smooth_trend, col = "gold3", lwd = 3.5)
legend("topright", legend = c("u_i Observation", "Manifold Signal Trend"),
       col = c("purple3", "gold3"), pch = c(19, NA), lwd = c(NA, 3.5), bty = "n")
grid()

dev.off()
cat("Figure 1 generated successfully as 'figure1_manifold_visualization.png' at 300 DPI.\n")
