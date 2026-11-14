# Element-form Normalized Inverse Gram Matrix (eNIGM) Prioritization

Computes the priority weight vector \\W\\ from a Pairwise Reciprocal
Matrix (PRM) `A` using the element-wise Normalized Inverse Gram Matrix
(NIGM) closed-form solution (Theorem 3 & Remark 1; Yuen, 2024).
Parameter `r != 0` (e.g., `r = 1`) guarantees invertibility even for
perfectly consistent PRMs.

## Usage

``` r
eNIGM(A, r = 1, detail = FALSE)
```

## Arguments

- A:

  Numeric square matrix of dimension `n x n` representing a Pairwise
  Reciprocal Matrix (PRM).

- r:

  Numeric scalar parameter for the Gram matrix shift (default: `1`).

- detail:

  Logical. If `TRUE`, returns a list containing intermediate results
  (`W`, `r`, `G`, `G_inv`, `rs_G_inv`) (default: `FALSE`).

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

# Compute priority weights via eNIGM (r = 1)
W <- eNIGM(A, r = 1)
print(W)
#> [1] 0.5858925 0.3014962 0.1126112
```
