# Construct Parametrized Design Matrix for NIGM

Converts an `n x n` Pairwise Reciprocal Matrix (PRM) `A` into a design
coefficient matrix `D` of dimension `(n(n-1) + 1) x n`, incorporating
pairwise comparison differential equations and appending a parametrized
shift row vector based on scalar parameter `r` (Section 3.2, Theorem 3;
Yuen, 2024).

## Usage

``` r
rDesignMatrix(A, r)
```

## Arguments

- A:

  Numeric square matrix of dimension `n x n` representing a Pairwise
  Reciprocal Matrix (PRM).

- r:

  Numeric scalar parameter determining the weight shift factor appended
  to the design matrix.

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
D <- rDesignMatrix(A, r = 1)
print(D)
#>      w1   w2
#> [1,]  1 -0.5
#> [2,] -2  1.0
#> [3,]  1  1.0
```
