# Check if Input Object is an Upper-Triangular Matrix (UTM) Ratings Vector

Evaluates whether an input object `X` is a 1D atomic vector whose length
`m` corresponds to the exact upper-triangular element count \\m =
n(n-1)/2\\ of a square Pairwise Reciprocal Matrix (PRM) of integer
dimension `n`.

## Usage

``` r
is.SizeOfUTM(X)
```

## Arguments

- X:

  A vector, matrix, or array candidate object to be evaluated.

## Value

A logical scalar: `TRUE` if `X` is a vector matching a valid
upper-triangular pairwise comparison ratings length for an integer
matrix order `n`, or `FALSE` otherwise.

## References

Yuen, K. K. F. (2024). Inverse Gram Matrix Methods for Prioritization in
Analytic Hierarchy Process: Explainability of Weighted Least Squares
Optimization Method. *arXiv preprint arXiv:2401.01190*.
<https://doi.org/10.48550/arXiv.2401.01190>

## Author

Kevin Kam Fung Yuen

## Examples

``` r
# Length 3 vector corresponds to a 3x3 PRM upper triangle: 3*(3-1)/2 = 3
is.SizeOfUTM(c(2, 5, 3))     # TRUE
#> [1] TRUE

# Length 4 vector does not correspond to an integer matrix order n
is.SizeOfUTM(c(1, 2, 3, 4))   # FALSE
#> [1] FALSE
```
