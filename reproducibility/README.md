# Replication Suite: Benchmarks and Simulations

This directory contains self-contained replication scripts to reproduce all empirical simulations, benchmark tables, and figures presented in:

> **The Local Manifold Dispersion Test: A Nonparametric Graph-Spectral Approach for Detecting Complex and High-Dimensional Heteroskedasticity**  
> *Ahmed Sattar Jabbar (Mustansiriyah University)*  

---

## Directory Contents

| File | Description |
| :--- | :--- |
| `simulation_experiments.R` | Master Monte Carlo simulation suite replicating benchmark simulations across DGPs 0 through 5. |
| `competitors.R` | Pure Base-R benchmark implementations of competing tests (Breusch–Pagan, White, Distance Covariance, and Goldfeld–Quandt). |
| `visualize_lmdt.R` | High-resolution script generating the three-panel geometric comparison figure (LMDT vs. Breusch–Pagan). |

---

## Prerequisites

Ensure the companion `LMDT` package is installed:

```r
# From parent directory
devtools::install(".")

# Or directly from GitHub
remotes::install_github("amjed-droid/LMDT")
```

---

## Replication Guide

### 1. Reproducing Simulation Experiments

Open an R session inside this directory and run:

```r
source("simulation_experiments.R")
```

This will run the Monte Carlo simulation across all scenarios:
- **$H_0$ Normal & $t_3$:** Validating exact kurtosis invariance and size calibration.
- **DGP 1 (Linear):** Classical sweet spot.
- **DGP 2 (Ring Manifold):** Demonstrating the total blindness of Breusch–Pagan (power = 0.000) versus LMDT (power > 0.98).
- **DGP 3 (Clustered):** Localized cluster variance.
- **DGP 4 (Latent Manifold, $p=30$):** Demonstrating White test collapse (`NaN`) and adaptive PCA denoising.
- **DGP 5 (Sparse Noise, $p=50$):** Demonstrating supervised metric learning in the presence of 48 irrelevant ambient noise dimensions.

### 2. Reproducing Figure 1 (Visualizations)

To generate the three-panel geometric figure:

```r
source("visualize_lmdt.R")
```

This generates `figure1_manifold_visualization.png` at 300 DPI, exactly matching Figure 1 in the JASA manuscript.

---

## Computational Environment

- **R Version:** $\ge$ 4.0.0
- **Dependencies:** Strictly Base R packages (`stats`, `graphics`, `grDevices`, `utils`) plus `LMDT`.
