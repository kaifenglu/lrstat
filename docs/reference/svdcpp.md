# Singular Value Decomposition of a Matrix

Computes the singular-value decomposition of a rectangular matrix.

## Usage

``` r
svdcpp(X, outtransform = TRUE, decreasing = TRUE)
```

## Arguments

- X:

  A numeric matrix whose SVD decomposition is to be computed.

- outtransform:

  Whether the orthogonal matrices composing of the left and right
  singular vectors are to be computed.

- decreasing:

  Whether the singular values should be sorted in decreasing order and
  the corresponding singular vectors rearranged accordingly.

## Value

A list with the following components:

- `d`: A vector containing the singular values of \\X\\.

- `U`: A matrix whose columns contain the left singular vectors of
  \\X\\.

- `V`: A matrix whose columns contain the right singular vectors of
  \\X\\.

## Details

Given \\A \in R^{m\times n} (m \geq n)\\, the following algorithm
overwrites \\A\\ with \\U^T A V = D\\, where \\U\in R^{m\times m}\\ is
orthogonal, \\V \in R^{n\times n}\\ is orthogonal, and \\D \in
R^{m\times n}\\ is diagonal.

## References

Gene N. Golub and Charles F. Van Loan. Matrix Computations, second
edition. Baltimore, Maryland: The John Hopkins University Press, 1989,
p.434.

## Author

Kaifeng Lu, <kaifenglu@gmail.com>

## Examples

``` r
A <- matrix(c(1,0,0,0, 1,2,0,0, 0,1,3,0, 0,0,1,4), 4, 4)
svdcpp(A)
#> $d
#> [1] 4.2600067 3.1073486 2.1117846 0.8585417
#> 
#> $U
#>             [,1]       [,2]       [,3]        [,4]
#> [1,] -0.01354261 -0.1354349  0.5426380  0.82886552
#> [2,] -0.10934067 -0.5184187  0.6673454 -0.52338971
#> [3,] -0.47016281 -0.7142294 -0.4819622  0.19114341
#> [4,] -0.87567582  0.4503064  0.1670527 -0.05009358
#> 
#> $V
#>             [,1]        [,2]        [,3]        [,4]
#> [1,] -0.00317901 -0.04358535  0.25695706  0.96543425
#> [2,] -0.05451258 -0.37725804  0.88897741 -0.25381867
#> [3,] -0.35676684 -0.85639145 -0.36866505  0.05828548
#> [4,] -0.93259621  0.34981478  0.08819481 -0.01075186
#> 
```
