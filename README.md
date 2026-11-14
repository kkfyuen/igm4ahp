# igm4ahp: Inverse Gram Matrix Methods for Analytic Hierarchy Process

igm4ahp is an R package implementing exact, closed-form Inverse Gram Matrix (IGM) prioritization methods for Pairwise Reciprocal Matrices (PRMs) in the Analytic Hierarchy Process (AHP).

Based on the theoretical framework established in Yuen (2024), the package derives exact global minimizers for the Weighted Least Squares (WLS) optimization problem in $O(n^3)$ algebraic operations—completely bypassing the need for non-linear iterative optimization solvers.

## Key Features

- **Exact Closed-Form Prioritization**: Solves global Weighted Least Squares (WLS) priority vectors directly without iteration.
- **Three IGM Algorithms**:
  - `PIGM()` / `pseudoIGM()`: Pseudo-Inverse Gram Matrix formulation via design matrix $D$.
  - `eNIGM()`: Element-form Normalized Inverse Gram Matrix computation.
  - `lagrangeIGM()`: Augmented Lagrangian-form Inverse Gram Matrix solving KKT block matrices.
- **Saaty Benchmark Utilities**:
  - `SaatyEigen()`: Principal eigenvector prioritization (Saaty's EVM).
  - `ConsistencyIndex()` & `ConsistencyRatio()`: Automated consistency evaluations (Saaty, 1980).
  - `WLSE()`: Objective function calculation for Weighted Least Squares Error.
- **Interactive Web Documentation Portal**: Explore live package vignettes at https://kkfyuen.github.io/igm4ahp.

## Installation

You can install the development version of igm4ahp directly from GitHub using `remotes` or `devtools`:

```R
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
remotes::install_github("kkfyuen/igm4ahp", build_vignettes = TRUE)
```

## Quick Start Example

```R
library(igm4ahp)

# Define a 3x3 Pairwise Reciprocal Matrix (PRM)
A <- matrix(c(
  1,   2, 5,
  1/2, 1, 3,
  1/5, 1/3, 1
), nrow = 3, byrow = TRUE)

# 1. Evaluate Matrix Consistency
cr <- ConsistencyRatio(A)
cat("Consistency Ratio (CR):", round(cr, 4), "\n")

# 2. Derive Priority Weights via Pseudo-IGM (PIGM)
w_pigm <- PIGM(A)
cat("PIGM Priority Weights:\n")
print(round(w_pigm, 4))

# 3. Derive Priority Weights via Element-form Normalized IGM (eNIGM)
w_enigm <- eNIGM(A)
cat("eNIGM Priority Weights:\n")
print(round(w_enigm, 4))

# 4. Compare with Saaty's Eigenvalue Method (EVM)
w_evm <- SaatyEigen(A)
cat("Saaty EVM Priority Weights:\n")
print(round(w_evm, 4))

# 5. Calculate Weighted Least Squares Error (WLSE)
cat("PIGM WLS Error:", round(WLSE(w_pigm, A), 6), "\n")
cat("EVM WLS Error: ", round(WLSE(w_evm, A), 6), "\n")
```

## Vignettes and Documentation Portal

View full online documentation and step-by-step numerical demonstrations at the interactive portal:

- Documentation Website: https://kkfyuen.github.io/igm4ahp

Or access compiled vignettes locally within R:

```R
# Browse all available package vignettes
vignette(package = "igm4ahp")

# Open Section 4 numerical examples vignette
vignette("section4_igm_examples", package = "igm4ahp")
```

## Citation

If you use `igm4ahp` in your research, please cite both the software release and the foundational research paper:

> Yuen, K. K. F. (2026). *igm4ahp: Inverse Gram Matrix Methods for Analytic Hierarchy Process* (Version v0.1.0) [Computer software]. Zenodo. https://doi.org/10.5281/zenodo.24223557

> Yuen, K. K. F. (2024). Inverse Gram Matrix Methods for Prioritization in Analytic Hierarchy Process: Explainability of Weighted Least Squares Optimization Method. *arXiv preprint arXiv:2401.01190*. https://doi.org/10.48550/arXiv.2401.01190

## License

This package is licensed under the GPL (>= 3) License.
