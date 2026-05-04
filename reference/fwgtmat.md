# Weight Matrix for All Intersection Hypotheses

Obtains the weight matrix for all intersection hypotheses.

## Usage

``` r
fwgtmat(w, G)
```

## Arguments

- w:

  The vector of weights for elementary hypotheses.

- G:

  The transition matrix.

## Value

The weight matrix starting with the global null hypothesis.

## Author

Kaifeng Lu, <kaifenglu@gmail.com>

## Examples

``` r
w <- c(0.5,0.5,0,0)
g <- matrix(c(0,0,1,0, 0,0,0,1, 0,1,0,0, 1,0,0,0),
            nrow=4, ncol=4, byrow=TRUE)
(wgtmat <- fwgtmat(w,g))
#>       [,1] [,2] [,3] [,4]
#>  [1,]  0.5  0.5  0.0  0.0
#>  [2,]  0.5  0.5  0.0  0.0
#>  [3,]  0.5  0.5  0.0  0.0
#>  [4,]  0.5  0.5  0.0  0.0
#>  [5,]  0.5  0.0  0.0  0.5
#>  [6,]  1.0  0.0  0.0  0.0
#>  [7,]  0.5  0.0  0.0  0.5
#>  [8,]  1.0  0.0  0.0  0.0
#>  [9,]  0.0  0.5  0.5  0.0
#> [10,]  0.0  0.5  0.5  0.0
#> [11,]  0.0  1.0  0.0  0.0
#> [12,]  0.0  1.0  0.0  0.0
#> [13,]  0.0  0.0  0.5  0.5
#> [14,]  0.0  0.0  1.0  0.0
#> [15,]  0.0  0.0  0.0  1.0
```
