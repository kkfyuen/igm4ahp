# Generate Random Pairwise Comparison Ratings Vector

Generates a random upper-triangular pairwise comparison ratings vector
sampled from Saaty's 1/9 to 9 ratio scale (or a custom evaluation scale)
for a Pairwise Reciprocal Matrix (PRM) of dimension `n` (Section 4;
Yuen, 2024).

## Usage

``` r
randomPRR(n = sample(3:15, 1), scale = c(1/(9:1), 2:9))
```

## Arguments

- n:

  Integer. Dimension of the Pairwise Reciprocal Matrix (PRM),
  representing the number of evaluation criteria. Defaults to a random
  integer sampled between 3 and 15.

- scale:

  Numeric vector. Set of discrete rating scale values from which
  pairwise comparisons are sampled. Defaults to Saaty's standard AHP
  scale `c(1/9, 1/8, ..., 1/2, 1, 2, ..., 9)`.

## Value

A numeric vector of length `n(n-1)/2` containing randomly sampled
pairwise comparison ratings.

## References

Yuen, K. K. F. (2024). Inverse Gram Matrix Methods for Prioritization in
Analytic Hierarchy Process: Explainability of Weighted Least Squares
Optimization Method. *arXiv preprint arXiv:2401.01190*.
<https://doi.org/10.48550/arXiv.2401.01190>

## Author

Kevin Kam Fung Yuen

## Examples

``` r
# Generate a random ratings vector for a 4x4 matrix
ratings <- randomPRR(n = 4)
print(ratings)
#> [1] 5.00 4.00 0.20 0.25 7.00 7.00
```
