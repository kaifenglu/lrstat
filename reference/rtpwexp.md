# Random Number Generation Function of Truncated Piecewise Exponential Distribution

Obtains random samples from a truncated piecewise exponential
distribution.

## Usage

``` r
rtpwexp(n, piecewiseSurvivalTime = 0, lambda = 0.0578, lowerBound = 0)
```

## Arguments

- n:

  The number of observations.

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

The random numbers from truncated piecewise exponential distribution.

## Author

Kaifeng Lu, <kaifenglu@gmail.com>

## Examples

``` r
rtpwexp(n = 10, piecewiseSurvivalTime = c(0, 6, 9, 15),
        lambda = c(0.025, 0.04, 0.015, 0.007))
#>  [1]   3.3678923 220.3965363  94.7421084   6.5258903   0.2970781  53.2995140
#>  [7]  61.9588298  13.8108359 152.1521031 175.0999428
```
