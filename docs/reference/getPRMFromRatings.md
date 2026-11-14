# Construct Pairwise Reciprocal Matrix from Ratings Vector

Reconstructs an `n x n` Pairwise Reciprocal Matrix (PRM) from an
upper-triangular pairwise ratings vector `V` of length `n(n-1)/2`.

## Usage

``` r
getPRMFromRatings(V, n = (1 + sqrt(1 + 8 * length(V)))/2)
```

## Arguments

- V:

  Numeric vector containing upper-triangular pairwise comparison ratings
  of length `n(n-1)/2`.

- n:

  Numeric scalar representing the dimension of the square PRM.
  Calculated automatically from the length of `V` by solving
  `n(n-1)/2 = length(V)`.

## Value

A numeric square matrix of dimension `n x n` representing the Pairwise
Reciprocal Matrix (PRM), or `NA_real_` with a warning if the input
vector length is invalid.

## References

Yuen, K. K. F. (2024). Inverse Gram Matrix Methods for Prioritization in
Analytic Hierarchy Process: Explainability of Weighted Least Squares
Optimization Method. *arXiv preprint arXiv:2401.01190*.
<https://doi.org/10.48550/arXiv.2401.01190>

## Author

Kevin Kam Fung Yuen

## Examples

``` r
# Vector of 3 upper-triangular ratings for a 3x3 matrix
V <- c(2, 5, 3)
A <- getPRMFromRatings(V)
print(A)
#>      [,1]      [,2] [,3]
#> [1,]  1.0 2.0000000    5
#> [2,]  0.5 1.0000000    3
#> [3,]  0.2 0.3333333    1
```
