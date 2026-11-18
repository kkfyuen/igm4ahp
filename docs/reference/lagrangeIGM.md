# Lagrangian Inverse Gram Matrix (LIGM) Prioritization

Computes the priority weight vector \\W\\ and exact Lagrange multiplier
\\\lambda\\ from a Pairwise Reciprocal Matrix (PRM) `A` using the
augmented Lagrangian Inverse Gram Matrix (LIGM) closed-form KKT system
solution \\\[W, \lambda\] = q \bar{G}^{-1}\\ (Theorem 4; Yuen, 2024),
eliminating the requirement for post-hoc vector normalization.

## Usage

``` r
lagrangeIGM(A, r = 0, detail = FALSE)

LIGM(A, r = 0, detail = FALSE)
```

## Arguments

- A:

  Numeric square matrix of dimension `n x n` representing a Pairwise
  Reciprocal Matrix (PRM).

- r:

  Numeric scalar parameter for the Gram matrix shift (default: `0`).

- detail:

  Logical. If `TRUE`, returns a list containing intermediate matrices
  (`W`, `G`, `G_inv`) (default: `FALSE`).

## Value

If `detail = FALSE`, a numeric vector of length `n + 1` containing the
priority weights \\W\\ and the Lagrange multiplier \\\lambda\\. If
`detail = TRUE`, a list containing intermediate calculation steps and
matrices.

## References

Yuen, K. K. F. (2024). Inverse Gram Matrix Methods for Prioritization in
Analytic Hierarchy Process: Explainability of Weighted Least Squares
Optimization Method. *arXiv preprint arXiv:2401.01190*.
<https://doi.org/10.48550/arXiv.2401.01190>

## Author

Kevin Kam Fung Yuen

## Examples

``` r
# Define a 3x3 Pairwise Reciprocal Matrix
A <- matrix(c(1, 2, 5,
              1/2, 1, 3,
              1/5, 1/3, 1), nrow = 3, byrow = TRUE)

# Compute priority weights via lagrangeIGM (r = 0)
W <- lagrangeIGM(A, r = 0)
print(W)
#> [1] 0.5858925 0.3014962 0.1126112 0.9976250
```
