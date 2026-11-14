# Check if a Value is a Whole Number

Tests whether numerical values in `x` are whole numbers (integers)
within a specified floating-point machine precision tolerance threshold.
Used to verify that the matrix order `n` derived from upper-triangular
pairwise ratings count \\m = n(n-1)/2\\ is a valid integer.

## Usage

``` r
isWholeNo(x, tolerance = sqrt(.Machine$double.eps))
```

## Arguments

- x:

  Numeric vector or scalar to be evaluated.

- tolerance:

  Numeric scalar specifying the numerical tolerance threshold (defaults
  to machine precision `sqrt(.Machine$double.eps)`).

## Value

A logical vector or scalar of the same dimensions as `x`, indicating
`TRUE` for whole numbers and `FALSE` otherwise.

## References

Yuen, K. K. F. (2024). Inverse Gram Matrix Methods for Prioritization in
Analytic Hierarchy Process: Explainability of Weighted Least Squares
Optimization Method. *arXiv preprint arXiv:2401.01190*.
<https://doi.org/10.48550/arXiv.2401.01190>

## Author

Kevin Kam Fung Yuen

## Examples

``` r
isWholeNo(4.0)  # TRUE
#> [1] TRUE
isWholeNo(4.1)  # FALSE
#> [1] FALSE
```
