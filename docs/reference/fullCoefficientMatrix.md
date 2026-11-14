# Construct Full Design (Coefficient) Matrix for Pairwise Comparisons

Converts an `n x n` Pairwise Reciprocal Matrix (PRM) `A` into the full
design coefficient matrix `D` of dimension `(n(n-1) + 1) x n` via
Algorithm 1 (Yuen, 2024). Stacks all `n(n-1)` pairwise weighted
differential equations (\\w_i - a\_{ij} w_j = 0\\) and appends the unit
simplex normalization row vector (\\e = \[1, \dots, 1\]\\, representing
\\\sum\_{i=1}^n w_i = 1\\).

## Usage

``` r
fullCoefficientMatrix(A)
```

## Arguments

- A:

  Numeric square matrix of dimension `n x n` representing a Pairwise
  Reciprocal Matrix (PRM).

## Value

A numeric design matrix of dimension `(n(n-1) + 1) x n` with column
names `w1, ..., wn`.

## References

Yuen, K. K. F. (2024). Inverse Gram Matrix Methods for Prioritization in
Analytic Hierarchy Process: Explainability of Weighted Least Squares
Optimization Method. *arXiv preprint arXiv:2401.01190*.
<https://doi.org/10.48550/arXiv.2401.01190>

## Author

Kevin Kam Fung Yuen

## Examples

``` r
A <- matrix(c(1, 2, 1/2, 1), nrow = 2)
D <- fullCoefficientMatrix(A)
print(D)
#>      w1   w2
#> [1,]  1 -0.5
#> [2,] -2  1.0
#> [3,]  1  1.0
```
