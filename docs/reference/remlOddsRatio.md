# REML Estimates of Individual Proportions With Specified Odds Ratio

Obtains the restricted maximum likelihood estimates of individual
proportions with specified odds ratio.

## Usage

``` r
remlOddsRatio(n1, y1, n2, y2, oddsRatioH0 = 1)
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

- oddsRatioH0:

  The specified odds ratio.

## Value

A vector of the restricted maximum likelihood estimates of the response
probabilities for the two treatment groups.

## Author

Kaifeng Lu, <kaifenglu@gmail.com>

## Examples

``` r
remlOddsRatio(n1 = 10, y1 = 4, n2 = 20, y2 = 2, oddsRatioH0 = 1.25)
#> [1] 0.2242871 0.1878564
```
