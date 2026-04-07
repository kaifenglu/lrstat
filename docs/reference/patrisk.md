# Probability of Being at Risk

Obtains the probability of being at risk at given analysis times.

## Usage

``` r
patrisk(time, piecewiseSurvivalTime, lambda, gamma)
```

## Arguments

- time:

  A vector of analysis times at which to calculate the probability of
  being at risk.

- piecewiseSurvivalTime:

  A vector that specifies the starting time of piecewise exponential
  survival time intervals. Must start with 0, e.g., `c(0, 6)` breaks the
  time axis into 2 event intervals: \\\[0, 6)\\ and \\\[6, \infty)\\.
  Defaults to 0 for exponential distribution.

- lambda:

  A vector of hazard rates for the event. One for each analysis time
  interval.

- gamma:

  The hazard rate for exponential dropout, or a vector of hazard rates
  for piecewise exponential dropout.

## Value

A vector of probabilities of being at risk at the specified analysis
times after enrollment for a patient in a treatment group with specified
piecewise exponential survival and dropout distributions.

## Author

Kaifeng Lu, <kaifenglu@gmail.com>

## Examples

``` r
# Piecewise exponential survival with hazard 0.0533 in the first 6
# months, and hazard 0.0309 thereafter, and 5% dropout by the end of
# 1 year.

patrisk(time = c(3, 9), piecewiseSurvivalTime = c(0, 6),
        lambda = c(0.0533, 0.0309), gamma = -log(1-0.05)/12)
#> [1] 0.8413704 0.6370100
```
