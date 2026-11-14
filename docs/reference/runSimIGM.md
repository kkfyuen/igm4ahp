# Simulation Verification for Inverse Gram Matrix (IGM) Methods

Executes a Monte Carlo simulation over `size` randomly generated
Pairwise Reciprocal Matrices (PRMs) to verify exact numerical weight
vector equivalence across closed-form Inverse Gram Matrix
variants—Pseudo-IGM
([`pseudoIGM`](https://kkfyuen.github.io/igm4ahp/reference/pseudoIGM.md)),
Element-form NIGM
([`eNIGM`](https://kkfyuen.github.io/igm4ahp/reference/eNIGM.md)), and
Lagrangian-IGM
([`lagrangeIGM`](https://kkfyuen.github.io/igm4ahp/reference/lagrangeIGM.md))—logging
verification metrics and Weighted Least Squares Error (WLSE) values to a
CSV file (Section 4.1; Yuen, 2024).

## Usage

``` r
runSimIGM(size, seedNo = 100, errorFile = "igmData.csv")
```

## Arguments

- size:

  Integer. Total number of random PRM simulation iterations to execute.

- seedNo:

  Integer. Seed number for the pseudo-random number generator to ensure
  reproducibility (default: `100`).

- errorFile:

  Character string. Output CSV file path for logging simulation results
  and discrepancies (default: `"igmData.csv"`).

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
  # Run a small simulation of 10 random PRM instances
  runSimIGM(size = 10, seedNo = 123, errorFile = tempfile(fileext = ".csv"))
} # }
```
