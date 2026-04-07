# Range of Accrual Duration for Target Number of Events

Obtains a range of accrual duration to reach the target number of
events.

## Usage

``` r
getDurationFromNevents(
  nevents = NA_real_,
  allocationRatioPlanned = 1,
  accrualTime = 0L,
  accrualIntensity = NA_real_,
  piecewiseSurvivalTime = 0L,
  stratumFraction = 1L,
  lambda1 = NA_real_,
  lambda2 = NA_real_,
  gamma1 = 0L,
  gamma2 = 0L,
  followupTime = NA_real_,
  fixedFollowup = FALSE,
  npoints = 23L
)
```

## Arguments

- nevents:

  The target number of events.

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

- followupTime:

  Follow-up time for the last enrolled subjects. Must be provided for
  fixed follow-up design.

- fixedFollowup:

  Whether a fixed follow-up design is used. Defaults to `FALSE` for
  variable follow-up.

- npoints:

  The number of accrual duration time points. Defaults to 23.

## Value

A data frame of the following variables:

- `nevents`: The target number of events.

- `fixedFollowup`: Whether a fixed follow-up design is used.

- `accrualDuration`: The accrual duration.

- `subjects`: The total number of subjects.

- `followupTime`: The follow-up time for the last enrolled subject.

- `studyDuration`: The study duration.

## Author

Kaifeng Lu, <kaifenglu@gmail.com>

## Examples

``` r
# Piecewise accrual, piecewise exponential survivals, and 5% dropout by
# the end of 1 year.

getDurationFromNevents(
  nevents = 80, allocationRatioPlanned = 1,
  accrualTime = seq(0, 8),
  accrualIntensity = 26/9*seq(1, 9),
  piecewiseSurvivalTime = c(0, 6),
  lambda1 = c(0.0533, 0.0309),
  lambda2 = c(0.0533, 0.0533),
  gamma1 = -log(1-0.05)/12,
  gamma2 = -log(1-0.05)/12,
  fixedFollowup = FALSE)
#>    nevents fixedFollowup accrualDuration  subjects followupTime studyDuration
#> 1       80         FALSE        7.307975  88.00653 1000.0000000    1007.30797
#> 2       80         FALSE        7.706836  97.22466   48.9720745      56.67891
#> 3       80         FALSE        8.105698 106.74814   34.2115019      42.31720
#> 4       80         FALSE        8.504559 117.11854   26.0766251      34.58118
#> 5       80         FALSE        8.903421 127.48894   20.9657948      29.86922
#> 6       80         FALSE        9.302282 137.85933   17.3556166      26.65790
#> 7       80         FALSE        9.701143 148.22973   14.6238062      24.32495
#> 8       80         FALSE       10.100005 158.60013   12.4589821      22.55899
#> 9       80         FALSE       10.498866 168.97052   10.6849566      21.18382
#> 10      80         FALSE       10.897728 179.34092    9.1933642      20.09109
#> 11      80         FALSE       11.296589 189.71132    7.9133804      19.20997
#> 12      80         FALSE       11.695451 200.08171    6.7965187      18.49197
#> 13      80         FALSE       12.094312 210.45211    5.8090767      17.90339
#> 14      80         FALSE       12.493173 220.82251    4.9448558      17.43803
#> 15      80         FALSE       12.892035 231.19290    4.1831506      17.07519
#> 16      80         FALSE       13.290896 241.56330    3.5008370      16.79173
#> 17      80         FALSE       13.689758 251.93370    2.8817430      16.57150
#> 18      80         FALSE       14.088619 262.30410    2.3141156      16.40273
#> 19      80         FALSE       14.487480 272.67449    1.7891415      16.27662
#> 20      80         FALSE       14.886342 283.04489    1.3000476      16.18639
#> 21      80         FALSE       15.285203 293.41529    0.8415218      16.12673
#> 22      80         FALSE       15.684065 303.78568    0.4093187      16.09338
#> 23      80         FALSE       16.082926 314.15608    0.0000000      16.08293
```
