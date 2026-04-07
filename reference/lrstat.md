# Number of Subjects Having an Event and Log-Rank Statistics

Obtains the number of subjects accrued, number of events, number of
dropouts, and number of subjects reaching the maximum follow-up in each
group, mean and variance of weighted log-rank score statistic, estimated
hazard ratio from weighted Cox regression and variance of log hazard
ratio estimate at given calendar times.

## Usage

``` r
lrstat(
  time = NA_real_,
  hazardRatioH0 = 1,
  allocationRatioPlanned = 1,
  accrualTime = 0L,
  accrualIntensity = NA_real_,
  piecewiseSurvivalTime = 0L,
  stratumFraction = 1L,
  lambda1 = NA_real_,
  lambda2 = NA_real_,
  gamma1 = 0L,
  gamma2 = 0L,
  accrualDuration = NA_real_,
  followupTime = NA_real_,
  fixedFollowup = FALSE,
  rho1 = 0,
  rho2 = 0,
  predictTarget = 2L
)
```

## Arguments

- time:

  A vector of calendar times at which to calculate the number of events
  and the mean and variance of log-rank test score statistic.

- hazardRatioH0:

  Hazard ratio under the null hypothesis for the active treatment versus
  control. Defaults to 1 for superiority test.

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

- stratumFraction:

  A vector of stratum fractions that sum to 1. Defaults to 1 for no
  stratification.

- lambda1:

  A vector of hazard rates for the event in each analysis time interval
  by stratum for the active treatment group.

- lambda2:

  A vector of hazard rates for the event in each analysis time interval
  by stratum for the control group.

- gamma1:

  The hazard rate for exponential dropout, a vector of hazard rates for
  piecewise exponential dropout applicable for all strata, or a vector
  of hazard rates for dropout in each analysis time interval by stratum
  for the active treatment group.

- gamma2:

  The hazard rate for exponential dropout, a vector of hazard rates for
  piecewise exponential dropout applicable for all strata, or a vector
  of hazard rates for dropout in each analysis time interval by stratum
  for the control group.

- accrualDuration:

  Duration of the enrollment period.

- followupTime:

  Follow-up time for the last enrolled subject.

- fixedFollowup:

  Whether a fixed follow-up design is used. Defaults to `FALSE` for
  variable follow-up.

- rho1:

  The first parameter of the Fleming-Harrington family of weighted
  log-rank test. Defaults to 0 for conventional log-rank test.

- rho2:

  The second parameter of the Fleming-Harrington family of weighted
  log-rank test. Defaults to 0 for conventional log-rank test.

- predictTarget:

  The target of prediction. Set `predictTarget = 1` to predict the
  number of events only. Set `predictTarget = 2` (default) to predict
  the number of events and log-rank score statistic mean and variance.
  Set `predictTarget = 3` to predict the number of events, log-rank
  score statistic mean and variance, and hazard ratio and variance of
  log hazard ratio.

## Value

A data frame containing the following variables if `predictTarget = 1`:

- `time`: The analysis time since trial start.

- `subjects`: The number of enrolled subjects.

- `nevents`: The total number of events.

- `nevents1`: The number of events in the active treatment group.

- `nevents2`: The number of events in the control group.

- `ndropouts`: The total number of dropouts.

- `ndropouts1`: The number of dropouts in the active treatment group.

- `ndropouts2`: The number of dropouts in the control group.

- `nfmax`: The total number of subjects reaching maximum follow-up.

- `nfmax1`: The number of subjects reaching maximum follow-up in the
  active treatment group.

- `nfmax2`: The number of subjects reaching maximum follow-up in the
  control group.

If `predictTarget = 2`, the following variables will also be included:

- `uscore`: The numerator of the log-rank test statistic.

- `vscore`: The variance of the log-rank score test statistic.

- `logRankZ`: The log-rank test statistic on the Z-scale.

- `hazardRatioH0`: The hazard ratio under the null hypothesis.

Furthermore, if `predictTarget = 3`, the following additional variables
will also be included:

- `HR`: The average hazard ratio from weighted Cox regression.

- `vlogHR`: The variance of log hazard ratio.

- `zlogHR`: The Z-statistic for log hazard ratio.

## Author

Kaifeng Lu, <kaifenglu@gmail.com>

## Examples

``` r
# Piecewise accrual, piecewise exponential survivals, and 5% dropout by
# the end of 1 year.

lrstat(time = c(22, 40), allocationRatioPlanned = 1,
       accrualTime = seq(0, 8),
       accrualIntensity = 26/9*seq(1, 9),
       piecewiseSurvivalTime = c(0, 6),
       lambda1 = c(0.0533, 0.0309),
       lambda2 = c(0.0533, 0.0533),
       gamma1 = -log(1-0.05)/12,
       gamma2 = -log(1-0.05)/12,
       accrualDuration = 22,
       followupTime = 18, fixedFollowup = FALSE)
#>   time subjects  nevents  nevents1  nevents2 ndropouts ndropouts1 ndropouts2
#> 1   22      468 154.4040  71.83183  82.57218  13.45783   6.835885    6.62195
#> 2   40      468 307.5078 138.43661 169.07120  29.03037  15.471551   13.55882
#>   nfmax nfmax1 nfmax2     uscore   vscore  logRankZ hazardRatioH0
#> 1     0      0      0  -6.404097 38.56497 -1.031244             1
#> 2     0      0      0 -24.351805 76.18939 -2.789870             1
```
