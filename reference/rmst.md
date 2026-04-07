# Restricted Mean Survival Time

Obtains the restricted mean survival time over an interval.

## Usage

``` r
rmst(t1 = 0, t2 = NA_real_, piecewiseSurvivalTime = 0L, lambda = NA_real_)
```

## Arguments

- t1:

  Lower bound of the analysis time interval.

- t2:

  Upper bound of the analysis time interval.

- piecewiseSurvivalTime:

  A vector that specifies the starting time of piecewise exponential
  survival time intervals. Must start with 0, e.g., `c(0, 6)` breaks the
  time axis into 2 event intervals: \\\[0, 6)\\ and \\\[6, \infty)\\.
  Defaults to 0 for exponential distribution.

- lambda:

  A vector of hazard rates for the event. One for each analysis time
  interval.

## Value

The integral of the survival function from `t1` to `t2`

## Author

Kaifeng Lu, <kaifenglu@gmail.com>

## Examples

``` r
rmst(t1 = 0, t2 = 7, piecewiseSurvivalTime = c(0, 6),
     lambda = c(0.0533, 0.0309))
#> [1] 5.850379
```
