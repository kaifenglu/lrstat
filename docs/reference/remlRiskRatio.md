# REML Estimates of Individual Proportions With Specified Risk Ratio

Obtains the restricted maximum likelihood estimates of individual
proportions with specified risk ratio.

## Usage

``` r
remlRiskRatio(n1, y1, n2, y2, riskRatioH0 = 1)
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

- riskRatioH0:

  The specified risk ratio.

## Value

A vector of the restricted maximum likelihood estimates of the response
probabilities for the two treatment groups.

## Author

Kaifeng Lu, <kaifenglu@gmail.com>

## Examples

``` r
remlRiskRatio(n1 = 10, y1 = 4, n2 = 20, y2 = 2, riskRatioH0 = 1.2)
#> [1] 0.2281748 0.1901457
```
