# Pseudo-Normalized Inverse Gram Matrix (PNIGM) Prioritization

Computes the priority weight vector \\W\\ from a Pairwise Reciprocal
Matrix (PRM) `A` using the Pseudo Inverse Gram Matrix (PIGM) closed-form
solution for Weighted Least Squares (WLS) prioritization (Theorem 1 &
Lemma 1; Yuen, 2024).

## Usage

``` r
PNIGM(A, detail = FALSE)
```

## Arguments

- A:

  Numeric square matrix of dimension `n x n` representing a Pairwise
  Reciprocal Matrix (PRM).

- detail:

  Logical. If `TRUE`, returns a list containing intermediate matrices
  (`D`, `G`, `Ginv`), unnormalized row sum vector `V`, and final
  priority vector `W` (default: `FALSE`).

## Value

If `detail = FALSE`, a numeric vector of normalized priority weights
`W`. If `detail = TRUE`, a list of intermediate calculation steps and
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

# Compute priority weights via PNIGM
W <- PNIGM(A)
print(W)
#>        w1        w2        w3 
#> 0.5858925 0.3014962 0.1126112 
```
