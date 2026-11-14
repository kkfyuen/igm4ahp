# Principal Eigenvector Priority Derivation (Saaty's EVM)

Computes the priority weight vector \\W\\ from a Pairwise Reciprocal
Matrix (PRM) using Saaty's Eigenvalue Method (EVM) by extracting and
normalizing the principal eigenvector corresponding to the maximum real
eigenvalue \\\lambda\_{\max}\\ (Saaty, 1977, 1980). Automatically
converts an upper-triangular ratings vector into a full `n x n` PRM if
supplied.

## Usage

``` r
SaatyEigen(PRM)
```

## Arguments

- PRM:

  A numeric square matrix of dimension `n x n` representing a Pairwise
  Reciprocal Matrix, or an upper-triangular pairwise comparison ratings
  vector.

## Value

A numeric vector of normalized priority weights `W` summing to unity.

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
# Define a 3x3 Pairwise Reciprocal Matrix
A <- matrix(c(1, 2, 5,
              1/2, 1, 3,
              1/5, 1/3, 1), nrow = 3, byrow = TRUE)

# Compute EVM priority weight vector
w_evm <- SaatyEigen(A)
print(w_evm)
#> [1] 0.5815521 0.3089956 0.1094523
```
