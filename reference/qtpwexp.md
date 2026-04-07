# Quantile Function of Truncated Piecewise Exponential Distribution

Obtains the quantile of a truncated piecewise exponential distribution.

## Usage

``` r
qtpwexp(
  p,
  piecewiseSurvivalTime = 0,
  lambda = 0.0578,
  lowerBound = 0,
  lower.tail = TRUE,
  log.p = FALSE
)
```

## Arguments

- p:

  The vector of probabilities.

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

- lower.tail:

  Logical; if TRUE (default), probabilities are P(X \<= x), otherwise,
  P(X \> x).

- log.p:

  Logical; if TRUE, probabilities p are given as log(p).

## Value

The quantile q such that P(X \> q \| X \> lowerBound) = 1 - p.

## Author

Kaifeng Lu, <kaifenglu@gmail.com>

## Examples

``` r
qtpwexp(p = c(0.205, 0.317), piecewiseSurvivalTime = c(0, 6, 9, 15),
        lambda = c(0.025, 0.04, 0.015, 0.007))
#> [1]  7.985329 18.037203
```
