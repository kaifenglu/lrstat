# Covariance Between Restricted Mean Survival Times

Obtains the covariance between restricted mean survival times at two
different time points.

## Usage

``` r
covrmst(
  t2 = NA_real_,
  tau1 = NA_real_,
  tau2 = NA_real_,
  allocationRatioPlanned = 1,
  accrualTime = 0L,
  accrualIntensity = NA_real_,
  piecewiseSurvivalTime = 0L,
  lambda1 = NA_real_,
  lambda2 = NA_real_,
  gamma1 = 0L,
  gamma2 = 0L,
  accrualDuration = NA_real_,
  maxFollowupTime = NA_real_
)
```

## Arguments

- t2:

  The calendar time for analysis 2.

- tau1:

  The milestone time for analysis 1.

- tau2:

  The milestone time for analysis 2.

- allocationRatioPlanned:

  Allocation ratio for the active treatment versus control. Defaults to
  1 for equal randomization.

- accrualTime:

  A vector that specifies the starting time of piecewise Poisson
  enrollment time intervals. Must start with 0, e.g., `c(0, 3)` breaks
  the time axis into 2 accrual intervals: \\\[0, 3)\\ and \\\[3,
  \infty)\\.

- accrualIntensity:

  A vector of accrual intensities. One for each accrual time interval.

- piecewiseSurvivalTime:

  A vector that specifies the starting time of piecewise exponential
  survival time intervals. Must start with 0, e.g., `c(0, 6)` breaks the
  time axis into 2 event intervals: \\\[0, 6)\\ and \\\[6, \infty)\\.
  Defaults to 0 for exponential distribution.

- lambda1:

  A vector of hazard rates for the event for the active treatment group.
  One for each analysis time interval.

- lambda2:

  A vector of hazard rates for the event for the control group. One for
  each analysis time interval.

- gamma1:

  The hazard rate for exponential dropout, or a vector of hazard rates
  for piecewise exponential dropout for the active treatment group.

- gamma2:

  The hazard rate for exponential dropout, or a vector of hazard rates
  for piecewise exponential dropout for the control group.

- accrualDuration:

  Duration of the enrollment period.

- maxFollowupTime:

  Follow-up time for the first enrolled subject. For fixed follow-up,
  `maxFollowupTime = minFollowupTime`. For variable follow-up,
  `maxFollowupTime = accrualDuration + minFollowupTime`.

## Value

The covariance between the restricted mean survival times for each
treatment group.

## Author

Kaifeng Lu, <kaifenglu@gmail.com>

## Examples

``` r
covrmst(t2 = 25, tau1 = 16, tau2 = 18, allocationRatioPlanned = 1,
        accrualTime = c(0, 3), accrualIntensity = c(10, 20),
        piecewiseSurvivalTime = c(0, 6),
        lambda1 = c(0.0533, 0.0309), lambda2 = c(0.0533, 0.0533),
        gamma1 = -log(1-0.05)/12, gamma2 = -log(1-0.05)/12,
        accrualDuration = 12, maxFollowupTime = 30)
#> [1] 0.3739127 0.3532650
```
