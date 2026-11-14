# Weighted Least Squares Error (WLSE) Objective Function

Computes the total Weighted Least Squares Error (WLSE) objective
function value \\f(W) = \sum\_{i=1}^n \sum\_{j=1}^n (w_i - a\_{ij}
w_j)^2\\ for a given candidate priority weight vector `W` and Pairwise
Reciprocal Matrix `A` (Chu et al., 1979; Yuen, 2024).

## Usage

``` r
WLSE(W, A)
```

## Arguments

- W:

  Numeric vector of length `n` representing candidate priority weights.

- A:

  Numeric square matrix of dimension `n x n` representing a Pairwise
  Reciprocal Matrix (PRM).

## Value

Numeric scalar. The total Weighted Least Squares Error value.

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
# Define a 3x3 Pairwise Reciprocal Matrix
A <- matrix(c(1, 2, 5,
              1/2, 1, 3,
              1/5, 1/3, 1), nrow = 3, byrow = TRUE)

# Candidate priority weight vector
W <- c(0.58, 0.28, 0.14)

# Calculate WLSE objective error
err <- WLSE(W, A)
print(err)
#> [1] 0.03725378
```
