# Simulation Verification Comparing Pseudo-IGM and Non-Linear WLS Optimization

Executes a Monte Carlo simulation benchmarking the closed-form priority
weight vectors derived by Pseudo-IGM
([`pseudoIGM`](https://kkfyuen.github.io/igm4ahp/reference/pseudoIGM.md))
against the numerical global optimization solver
([`WLSopt`](https://kkfyuen.github.io/igm4ahp/reference/WLSopt.md)) via
`nloptr`, verifying that closed-form IGM formulas consistently achieve
the exact global minimum Weighted Least Squares Error (WLSE) without
relying on iterative numerical solvers (Section 4.2; Yuen, 2024).

## Usage

``` r
runWLSEOpt(
  size = 10,
  seedNo = 100,
  sigDigits = 4,
  errorFile = "wlseOptData.txt"
)
```

## Arguments

- size:

  Integer. Total number of random PRM simulation instances to evaluate
  (default: `10`).

- seedNo:

  Integer. Seed number for pseudo-random number generator initialization
  to ensure reproducibility (default: `100`).

- sigDigits:

  Integer. Number of significant decimal places used for rounding
  priority weight vectors during equality verification (default: `4`).

- errorFile:

  Character string. Output text file path for logging simulation results
  and discrepancy metrics (default: `"wlseOptData.txt"`).

## Value

Logical. Returns `FALSE` if all simulation iterations pass without
discrepancy, or `TRUE` if a discrepancy is detected.

## References

Yuen, K. K. F. (2024). Inverse Gram Matrix Methods for Prioritization in
Analytic Hierarchy Process: Explainability of Weighted Least Squares
Optimization Method. *arXiv preprint arXiv:2401.01190*.
<https://doi.org/10.48550/arXiv.2401.01190>

## Author

Kevin Kam Fung Yuen

## Examples

``` r
if (FALSE) { # \dontrun{
  # Run simulation comparing Pseudo-IGM against numerical solver for 5 random PRMs
  runWLSEOpt(size = 5, seedNo = 123, errorFile = tempfile(fileext = ".txt"))
} # }
```
