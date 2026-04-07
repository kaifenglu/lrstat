# Density Function of Truncated Piecewise Exponential Distribution

Obtains the density of a truncated piecewise exponential distribution.

## Usage

``` r
dtpwexp(
  q,
  piecewiseSurvivalTime = 0,
  lambda = 0.0578,
  lowerBound = 0,
  log.d = FALSE
)
```

## Arguments

- q:

  The vector of quantiles.

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

- log.d:

  Logical; if TRUE, densities d are given as log(d).

## Value

The density d such that d = lambda(q) \* P(X \> q \| X \> lowerBound).

## Author

Kaifeng Lu, <kaifenglu@gmail.com>

## Examples

``` r
dtpwexp(q = c(8, 18), piecewiseSurvivalTime = c(0, 6, 9, 15),
        lambda = c(0.025, 0.04, 0.015, 0.007))
#> [1] 0.031781344 0.004782245
```
