# Adjusted p-Values for Holm, Hochberg, and Hommel Procedures

Obtains the adjusted p-values for possibly truncated Holm, Hochberg, and
Hommel procedures.

## Usage

``` r
ftrunc(p, test = "hommel", gamma = 1)
```

## Arguments

- p:

  The raw p-values for elementary hypotheses.

- test:

  The test to use, e.g., "holm", "hochberg", or "hommel" (default).

- gamma:

  The value of the truncation parameter. Defaults to 1 for the regular
  Holm, Hochberg, or Hommel procedure.

## Value

A matrix of adjusted p-values.

## References

Alex Dmitrienko, Ajit C. Tamhane, and Brian L. Wiens. General multistage
gatekeeping procedures. Biometrical Journal. 2008; 5:667-677.

## Author

Kaifeng Lu, <kaifenglu@gmail.com>

## Examples

``` r
pvalues <- matrix(c(0.01,0.005,0.015,0.022, 0.02,0.015,0.010,0.023),
                  nrow=2, ncol=4, byrow=TRUE)
ftrunc(pvalues, "hochberg")
#>       [,1]  [,2]  [,3]  [,4]
#> [1,] 0.022 0.020 0.022 0.022
#> [2,] 0.023 0.023 0.023 0.023
```
