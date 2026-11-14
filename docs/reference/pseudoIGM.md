# Pseudo Inverse Gram Matrix (PIGM) Prioritization

Computes the priority weight vector \\W\\ from a Pairwise Reciprocal
Matrix (PRM) `A` using the closed-form Pseudo Inverse Gram Matrix (PIGM)
formula \\W = \frac{e G^{-1}}{e G^{-1} e^T}\\ where \\G = D^T D\\
(Theorem 1; Yuen, 2024).

## Usage

``` r
pseudoIGM(A)
```

## Arguments

- A:

  Numeric square matrix of dimension `n x n` representing a Pairwise
  Reciprocal Matrix (PRM).

## Value

A numeric vector of normalized priority weights `W`.

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

# Compute priority weights via pseudoIGM
W <- pseudoIGM(A)
print(W)
#> [1] 0.5858925 0.3014962 0.1126112
```
