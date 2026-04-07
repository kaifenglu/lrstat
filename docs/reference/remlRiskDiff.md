# REML Estimates of Individual Proportions With Specified Risk difference

Obtains the restricted maximum likelihood estimates of individual
proportions with specified risk difference.

## Usage

``` r
remlRiskDiff(n1, y1, n2, y2, riskDiffH0 = 0)
```

## Arguments

- n1:

  The sample size for the active treatment group.

- y1:

  The number of responses for the active treatment group.

- n2:

  The sample size for the control group.

- y2:

  The number of responses for the control group.

- riskDiffH0:

  The specified risk difference.

## Value

A vector of the restricted maximum likelihood estimates of the response
probabilities for the two treatment groups.

## Author

Kaifeng Lu, <kaifenglu@gmail.com>

## Examples

``` r
remlRiskDiff(n1 = 10, y1 = 4, n2 = 20, y2 = 0, riskDiffH0 = 0.1)
#> [1] 0.14332071 0.04332071
```
