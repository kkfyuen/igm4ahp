# Calculate Consistency Index (CI) for Pairwise Reciprocal Matrix

Computes Saaty's Consistency Index \\CI = (\lambda\_{\max} - n) / (n -
1)\\ for an `n x n` Pairwise Reciprocal Matrix (PRM), where
\\\lambda\_{\max}\\ is the maximum real eigenvalue and \\n\\ is the
matrix order (Saaty, 1977, 1980).

## Usage

``` r
ConsistencyIndex(PRM)
```

## Arguments

- PRM:

  A numeric square matrix of dimension `n x n` representing a Pairwise
  Reciprocal Matrix.

## Value

A numeric scalar representing the Consistency Index (CI), rounded to 10
decimal places.

## References

Saaty, T. L. (1977). A scaling method for priorities in hierarchical
structures. *Journal of Mathematical Psychology*, 15(3), 234-281.

Saaty, T. L. (1980). *The Analytic Hierarchy Process: Planning, Priority
Setting, Resource Allocation*. McGraw-Hill, New York.

Yuen, K. K. F. (2024). Inverse Gram Matrix Methods for Prioritization in
Analytic Hierarchy Process: Explainability of Weighted Least Squares
Optimization Method. *arXiv preprint arXiv:2401.01190*.
<https://doi.org/10.48550/arXiv.2401.01190>

## Author

Kevin Kam Fung Yuen

## Examples

``` r
A <- matrix(c(1, 2, 5,
              1/2, 1, 3,
              1/5, 1/3, 1), nrow = 3, byrow = TRUE)
ci <- ConsistencyIndex(A)
print(ci)
#> [1] 0.001847299
```
