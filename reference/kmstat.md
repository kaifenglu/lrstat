# Stratified Difference in Milestone Survival Probabilities

Obtains the stratified milestone survival probabilities and difference
in milestone survival probabilities at given calendar times.

## Usage

``` r
kmstat(
  time = NA_real_,
  milestone = NA_real_,
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
  fixedFollowup = FALSE
)
```

## Arguments

- time:

  A vector of calendar times for data cut.

- milestone:

  The milestone time at which to calculate the survival probability.

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

## Value

A data frame containing the following variables:

- `time`: The calendar time since trial start.

- `subjects`: The number of enrolled subjects.

- `nevents`: The total number of events.

- `nevents1`: The number of events in the active treatment group.

- `nevents2`: The number of events in the control group.

- `ndropouts`: The total number of dropouts.

- `ndropouts1`: The number of dropouts in the active treatment group.

- `ndropouts2`: The number of dropouts in the control group.

- `milestone`: The milestone time relative to randomization.

- `nmilestone`: The total number of subjects reaching milestone.

- `nmilestone1`: The number of subjects reaching milestone in the active
  treatment group.

- `nmiletone2`: The number of subjects reaching milestone in the control
  group.

- `surv1`: The milestone survival probability for the treatment group.

- `surv2`: The milestone survival probability for the control group.

- `survDiff`: The difference in milestone survival probabilities, i.e.,
  `surv1 - surv2`.

- `vsurv1`: The variance for `surv1`.

- `vsurv2`: The variance for `surv2`.

- `vsurvDiff`: The variance for `survDiff`.

- `information`: The information for `survDiff`, equal to `1/vsurvDiff`.

- `survDiffZ`: The Z-statistic value, i.e., `survDiff/sqrt(vsurvDiff)`.

## Author

Kaifeng Lu, <kaifenglu@gmail.com>

## Examples

``` r
# Piecewise accrual, piecewise exponential survivals, and 5% dropout by
# the end of 1 year.

kmstat(time = c(22, 40),
       milestone = 18,
       allocationRatioPlanned = 1,
       accrualTime = seq(0, 8),
       accrualIntensity = 26/9*seq(1, 9),
       piecewiseSurvivalTime = c(0, 6),
       stratumFraction = c(0.2, 0.8),
       lambda1 = c(0.0533, 0.0309, 1.5*0.0533, 1.5*0.0309),
       lambda2 = c(0.0533, 0.0533, 1.5*0.0533, 1.5*0.0533),
       gamma1 = -log(1-0.05)/12,
       gamma2 = -log(1-0.05)/12,
       accrualDuration = 22,
       followupTime = 18, fixedFollowup = FALSE)
#>   time subjects  nevents  nevents1 nevents2 ndropouts ndropouts1 ndropouts2
#> 1   22      468 195.2491  91.89764 103.3514  12.16925   6.202213   5.967036
#> 2   40      468 356.3965 164.83437 191.5622  24.15559  13.010003  11.145587
#>   milestone nmilestone nmilestone1 nmilestone2     surv1     surv2  survDiff
#> 1        18   8.700526    5.138325    3.562201 0.3841805 0.2663374 0.1178431
#> 2        18 140.948521   83.240865   57.707656 0.3841805 0.2663374 0.1178431
#>        vsurv1       vsurv2   vsurvDiff information survDiffZ
#> 1 0.003856622 0.0038944814 0.007751103    129.0139  1.338512
#> 2 0.001037023 0.0008598885 0.001896912    527.1726  2.705706
```
