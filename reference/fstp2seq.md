# Adjusted p-Values for Stepwise Testing Procedures for Two Sequences

Obtains the adjusted p-values for the stepwise gatekeeping procedures
for multiplicity problems involving two sequences of hypotheses.

## Usage

``` r
fstp2seq(p, gamma, test = "hochberg", retest = TRUE)
```

## Arguments

- p:

  The raw p-values for elementary hypotheses.

- gamma:

  The truncation parameters for each family. The truncation parameter
  for the last family is automatically set to 1.

- test:

  The component multiple testing procedure. It is either "Holm" or
  "Hochberg", and it defaults to "Hochberg".

- retest:

  Whether to allow retesting. It defaults to `TRUE`.

## Value

A matrix of adjusted p-values.

## Author

Kaifeng Lu, <kaifenglu@gmail.com>

## Examples

``` r
p <- c(0.0194, 0.0068, 0.0271, 0.0088, 0.0370, 0.0018, 0.0814, 0.0066)
gamma <- c(0.6, 0.6, 0.6, 1)
fstp2seq(p, gamma, test="hochberg", retest=1)
#> [1] 0.02425 0.01360 0.03300 0.02425 0.03700 0.02425 0.08140 0.03300
```
