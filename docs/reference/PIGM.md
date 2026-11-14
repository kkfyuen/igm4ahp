# Pseudo Inverse Gram Matrix (PIGM) Prioritization

Computes the priority weight vector \\W\\ from a Pairwise Reciprocal
Matrix (PRM) `A` by calling
[`pseudoIGM`](https://kkfyuen.github.io/igm4ahp/reference/pseudoIGM.md).

## Usage

``` r
PIGM(A)
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
A <- matrix(c(1, 2, 5,
              1/2, 1, 3,
              1/5, 1/3, 1), nrow = 3, byrow = TRUE)
W <- PIGM(A)
print(W)
#> [1] 0.5858925 0.3014962 0.1126112
```
