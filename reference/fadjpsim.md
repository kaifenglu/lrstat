# Adjusted p-Values for Simes-Based Graphical Approaches

Obtains the adjusted p-values for graphical approaches using weighted
Simes tests.

## Usage

``` r
fadjpsim(wgtmat, p, family = NULL)
```

## Arguments

- wgtmat:

  The weight matrix for intersection hypotheses.

- p:

  The raw p-values for elementary hypotheses.

- family:

  The matrix of family indicators for elementary hypotheses.

## Value

A matrix of adjusted p-values.

## References

Frank Bretz, Martin Posch, Ekkehard Glimm, Florian Klinglmueller, Willi
Maurer, and Kornelius Rohmeyer. Graphical approach for multiple
comparison procedures using weighted Bonferroni, Simes, or parameter
tests. Biometrical Journal. 2011; 53:894-913.

Kaifeng Lu. Graphical approaches using a Bonferroni mixture of weighted
Simes tests. Statistics in Medicine. 2016; 35:4041-4055.

## Author

Kaifeng Lu, <kaifenglu@gmail.com>

## Examples

``` r
pvalues <- matrix(c(0.01,0.005,0.015,0.022, 0.02,0.015,0.010,0.023),
                  nrow=2, ncol=4, byrow=TRUE)
w <- c(0.5,0.5,0,0)
g <- matrix(c(0,0,1,0,0,0,0,1,0,1,0,0,1,0,0,0),
            nrow=4, ncol=4, byrow=TRUE)
wgtmat <- fwgtmat(w,g)

family <- matrix(c(1,1,0,0,0,0,1,1), nrow=2, ncol=4, byrow=TRUE)
fadjpsim(wgtmat, pvalues, family)
#>      [,1] [,2]  [,3]  [,4]
#> [1,] 0.02 0.01 0.022 0.022
#> [2,] 0.04 0.02 0.040 0.040
```
