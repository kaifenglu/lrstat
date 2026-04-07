# REML Estimates of Individual Rates With Specified Rate Difference

Obtains the restricted maximum likelihood estimates of individual
proportions with specified rate difference.

## Usage

``` r
remlRateDiff(t1, y1, t2, y2, rateDiffH0 = 0)
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

- rateDiffH0:

  The specified rate difference.

## Value

A vector of the restricted maximum likelihood estimates of the incidence
rates for the two treatment groups.

## Author

Kaifeng Lu, <kaifenglu@gmail.com>

## Examples

``` r
remlRateDiff(t1 = 10, y1 = 4, t2 = 20, y2 = 2, rateDiffH0 = 0.1)
#> [1] 0.2457427 0.1457427
```
