# Mean and Variance of Truncated Piecewise Exponential Distribution

Obtains the mean and variance from a truncated piecewise exponential
distribution.

## Usage

``` r
mtpwexp(piecewiseSurvivalTime = 0L, lambda = NA_real_, lowerBound = 0)
```

## Arguments

- piecewiseSurvivalTime:

  A vector that specifies the starting time of piecewise exponential
  survival time intervals. Must start with 0, e.g., `c(0, 6)` breaks the
  time axis into 2 event intervals: \\\[0, 6)\\ and \\\[6, \infty)\\.
  Defaults to 0 for exponential distribution.

- lambda:

  A vector of hazard rates for the event. One for each analysis time
  interval.

- lowerBound:

  The left truncation time point for the survival time. Defaults to 0
  for no truncation.

## Value

A list with two components, one for the mean, and the other for the
variance of the truncated piecewise exponential distribution.

## Author

Kaifeng Lu, <kaifenglu@gmail.com>

## Examples

``` r
mtpwexp(piecewiseSurvivalTime = c(0, 6, 9, 15),
        lambda = c(0.025, 0.04, 0.015, 0.007))
#> $mean
#> [1] 112.0532
#> 
#> $variance
#> [1] 19084.4
#> 
```
