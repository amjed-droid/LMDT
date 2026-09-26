#' Print Method for LMDT Results
#'
#' Formats and prints the results of a Local Manifold Dispersion Test.
#'
#' @param x An object of class \code{"lmdt"}.
#' @param digits Number of significant digits to print.
#' @param ... Additional arguments.
#'
#' @return Invisibly returns \code{x}.
#' @export
print.lmdt <- function(x, digits = 4, ...) {
  cat("\n")
  cat("  ", x$method_name, "\n", sep = "")
  cat("\n")
  cat("data:  ", x$data.name, "\n", sep = "")
  cat("Z_LMD = ", format(round(x$statistic, digits)), 
      ", T_LMD = ", format(round(x$tlmd, digits)),
      ", p-value = ", format.pval(x$p.value, digits = digits), "\n", sep = "")
  cat("alternative hypothesis: ", x$alternative, "\n", sep = "")
  if (!is.null(x$B)) {
    cat("replications: B = ", x$B, ", neighborhood size: k = ", x$k, "\n", sep = "")
  }
  cat("denoising scheme: ", x$denoise, "\n", sep = "")
  cat("\n")
  invisible(x)
}

#' Summary Method for LMDT Results
#'
#' Provides a detailed summary of the LMDT test including spectral properties.
#'
#' @param object An object of class \code{"lmdt"}.
#' @param ... Additional arguments.
#'
#' @return Invisibly returns \code{object}.
#' @export
summary.lmdt <- function(object, ...) {
  print.lmdt(object, ...)
  cat("Spectral Graph Properties:\n")
  cat("  Trace of squared normalized adjacency Tr(S^2): ", round(object$tr_s2, 4), "\n", sep = "")
  cat("  Asymptotic null standard deviation sigma_0:     ", round(object$sigma0, 6), "\n", sep = "")
  cat("\n")
  invisible(object)
}

#' Plot Diagnostic Visualizations for LMDT
#'
#' Produces diagnostic plots for an object of class \code{"lmdt"}.
#'
#' @param x An object of class \code{"lmdt"}.
#' @param ... Additional graphical parameters.
#'
#' @return Invisibly returns \code{x}.
#' @export
#' @examples
#' set.seed(42)
#' fit <- lm(rnorm(50) ~ matrix(rnorm(100), 50, 2))
#' res <- lmdt_test(fit, method = "asymptotic")
#' plot(res)
plot.lmdt <- function(x, ...) {
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par))
  
  has_resamples <- !is.null(x$resamples) && length(x$resamples) > 0
  graphics::par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3, 1))
  
  if (has_resamples) {
    stats <- x$resamples
    xlims <- range(c(stats, x$tlmd))
    xlims[1] <- xlims[1] - 0.05 * diff(xlims)
    xlims[2] <- xlims[2] + 0.05 * diff(xlims)
    
    graphics::hist(stats, breaks = 20, col = "skyblue2", border = "white",
                   main = "Null Distribution of T_LMD",
                   xlab = "T_LMD (Resamples under H0)", xlim = xlims, ...)
    graphics::abline(v = x$tlmd, col = "red3", lwd = 2.5, lty = 1)
    graphics::abline(v = 1.0, col = "darkgreen", lwd = 1.5, lty = 2)
    graphics::legend("topright",
                     legend = c(paste0("Observed T = ", round(x$tlmd, 3)), "Null = 1.0"),
                     col = c("red3", "darkgreen"), lwd = c(2.5, 1.5), lty = c(1, 2),
                     bty = "n", cex = 0.8)
  } else {
    z_seq <- seq(-4, 4, length.out = 200)
    y_seq <- stats::dnorm(z_seq)
    graphics::plot(z_seq, y_seq, type = "l", lwd = 2, col = "steelblue",
                   main = "Asymptotic Distribution (Z_LMD)",
                   xlab = "Standardized Z", ylab = "Density", ...)
    graphics::abline(v = x$statistic, col = "red3", lwd = 2.5)
    graphics::legend("topright",
                     legend = c(paste0("Z_obs = ", round(x$statistic, 3)),
                                paste0("p = ", format.pval(x$p.value, digits = 3))),
                     col = c("red3", "transparent"), lwd = c(2.5, 0), bty = "n", cex = 0.8)
  }
  
  ord <- order(x$u)
  graphics::plot(x$u[ord], pch = 20, col = "darkblue",
                 main = "Ordered Log-Dispersion",
                 xlab = "Rank Index", ylab = expression(ln(abs(e[i]) + c)), ...)
  graphics::grid()
  
  invisible(x)
}

