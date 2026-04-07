# REML Estimates of Individual Rates With Specified Rate Ratio

Obtains the restricted maximum likelihood estimates of individual
proportions with specified rate ratio.

## Usage

``` r
remlRateRatio(t1, y1, t2, y2, rateRatioH0 = 1)
```

## Arguments

- t1:

  The exposure for the active treatment group.

- y1:

  The number of events for the active treatment group.

- t2:

  The exposure for the control group.

- y2:

  The number of events for the control group.

- rateRatioH0:

  The specified rate ratio.

## Value

A vector of the restricted maximum likelihood estimates of the incidence
rates for the two treatment groups.

## Author

Kaifeng Lu, <kaifenglu@gmail.com>

## Examples

``` r
remlRateRatio(t1 = 10, y1 = 4, t2 = 20, y2 = 2, rateRatioH0 = 1.1)
#> [1] 0.2129032 0.1935484
```
