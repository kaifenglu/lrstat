# Adjusted p-Values for Bonferroni-Based Graphical Approaches

Obtains the adjusted p-values for graphical approaches using weighted
Bonferroni tests.

## Usage

``` r
fadjpbon(w, G, p)
```

## Arguments

- w:

  The vector of initial weights for elementary hypotheses.

- G:

  The initial transition matrix.

- p:

  The raw p-values for elementary hypotheses.

## Value

A matrix of adjusted p-values.

## References

Frank Bretz, Willi Maurer, Werner Brannath and Martin Posch. A graphical
approach to sequentially rejective multiple test procedures. Statistics
in Medicine. 2009; 28:586-604.

## Author

Kaifeng Lu, <kaifenglu@gmail.com>

## Examples

``` r
pvalues <- matrix(c(0.01,0.005,0.015,0.022, 0.02,0.015,0.010,0.023),
                  nrow=2, ncol=4, byrow=TRUE)
w <- c(0.5,0.5,0,0)
g <- matrix(c(0,0,1,0,0,0,0,1,0,1,0,0,1,0,0,0),
            nrow=4, ncol=4, byrow=TRUE)
fadjpbon(w, g, pvalues)
#>      [,1] [,2] [,3] [,4]
#> [1,] 0.02 0.01 0.03 0.03
#> [2,] 0.04 0.03 0.04 0.04
```
