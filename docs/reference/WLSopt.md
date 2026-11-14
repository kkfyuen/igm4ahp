# Numerical Weighted Least Squares (WLS) Optimization Solver for AHP

Solves the non-linear Weighted Least Squares (WLS) optimization problem
for a Pairwise Reciprocal Matrix (PRM) in the Analytic Hierarchy Process
(AHP), originally formulated by Chu et al. (1979). Uses stochastic or
derivative-free algorithms via the `nloptr` package (default:
`"NLOPT_GN_ISRES"`) to find the global minimizer of the Weighted Least
Squares Error (WLSE), serving as a numerical benchmark against exact
closed-form Inverse Gram Matrix (IGM) methods.

## Usage

``` r
WLSopt(A, maxeval = 5e+05, xtol_rel = 1e-06, algorithm = "NLOPT_GN_ISRES")
```

## Arguments

- A:

  A numeric square matrix of dimension `n x n` representing a Pairwise
  Reciprocal Matrix (PRM).

- maxeval:

  Integer. Maximum number of function evaluations allowed for the
  optimization solver (default: `500000`).

- xtol_rel:

  Numeric. Relative convergence stopping tolerance for the objective
  function (default: `1.0e-6`).

- algorithm:

  Character string. The algorithm code recognized by `nloptr` (default:
  `"NLOPT_GN_ISRES"`).

## Value

An object of class `nloptr` containing the optimization results,
including the optimal priority weight vector (`$solution`) and the
objective function value (`$objective`).

## References

Chu, A. T. W., Kalaba, R. E., & Spingarn, K. (1979). A comparison of two
methods for determining the weights of belonging to fuzzy sets. *Journal
of Optimization Theory and Applications*, 27(4), 531-538.

Yuen, K. K. F. (2024). Inverse Gram Matrix Methods for Prioritization in
Analytic Hierarchy Process: Explainability of Weighted Least Squares
Optimization Method. *arXiv preprint arXiv:2401.01190*.
<https://doi.org/10.48550/arXiv.2401.01190>

## Author

Kevin Kam Fung Yuen

## Examples

``` r
if (FALSE) { # \dontrun{
  # Define a sample 3x3 Pairwise Reciprocal Matrix
  A <- matrix(c(1, 2, 5,
                1/2, 1, 3,
                1/5, 1/3, 1), nrow = 3, byrow = TRUE)

  # Run WLS numerical optimization
  res <- WLSopt(A)
  print(res$solution)
} # }
```
