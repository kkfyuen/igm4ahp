# Calculate Consistency Ratio (CR) for Pairwise Reciprocal Matrix

Computes Saaty's Consistency Ratio \\CR = CI / RI\\ for a Pairwise
Reciprocal Matrix (PRM), where \\CI\\ is the Consistency Index and
\\RI\\ is the average Random Index for matrix size \\n \le 15\\ (Saaty,
1980, p. 21).

## Usage

``` r
ConsistencyRatio(PRM, randomIndex = NULL)
```

## Arguments

- PRM:

  A numeric square matrix of dimension `n x n` representing a Pairwise
  Reciprocal Matrix.

- randomIndex:

  Optional numeric vector providing custom Random Index (RI) values for
  matrix dimensions `1` to `15`. Defaults to Saaty's standard RI table.

## Value

A numeric scalar representing the calculated Consistency Ratio (CR), or
`NA_real_` with a warning if `n > 15` without custom RI values.

## References

Saaty, T. L. (1980). *The Analytic Hierarchy Process: Planning, Priority
Setting, Resource Allocation*. McGraw-Hill, New York (p. 21).

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
cr <- ConsistencyRatio(A)
print(cr)
#> [1] 0.003184998
```
