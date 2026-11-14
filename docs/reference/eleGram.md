# Element-wise Gram Matrix Construction

Constructs the `n x n` Gram matrix \\G = (g\_{ij})\\ directly from
element-wise formulas derived from the Weighted Least Squares (WLS)
objective function (Theorem 2; Yuen, 2024), bypassing the construction
of design matrix `D`, with an optional scalar shift parameter `r`.

## Usage

``` r
eleGram(A, r = 0)
```

## Arguments

- A:

  Numeric square matrix of dimension `n x n` representing a Pairwise
  Reciprocal Matrix (PRM).

- r:

  Numeric scalar shift parameter added to all matrix entries (default:
  `0`).

## Value

A numeric square matrix of dimension `n x n` representing the
constructed Gram matrix.

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

# Compute element-form Gram matrix (r = 0)
G <- eleGram(A, r = 0)
print(G)
#>       [,1]      [,2]      [,3]
#> [1,]  2.29 -2.500000 -5.200000
#> [2,] -2.50  6.111111 -3.333333
#> [3,] -5.20 -3.333333 36.000000
```
