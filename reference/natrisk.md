# Number of Subjects at Risk

Obtains the number of subjects at risk at given analysis times for each
treatment group.

## Usage

``` r
natrisk(
  t = NA_real_,
  allocationRatioPlanned = 1,
  accrualTime = 0L,
  accrualIntensity = NA_real_,
  piecewiseSurvivalTime = 0L,
  lambda1 = NA_real_,
  lambda2 = NA_real_,
  gamma1 = 0L,
  gamma2 = 0L,
  accrualDuration = NA_real_,
  maxFollowupTime = NA_real_,
  time = NA_real_
)
```

## Arguments

- t:

  A vector of analysis times at which to calculate the number of
  patients at risk.

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

- time:

  Calendar time for the analysis.

## Value

A matrix of the number of patients at risk at the specified analysis
times (row) for each treatment group (column).

## Details

For a given treatment group \\g\\ and calendar time \\\tau\\, the number
of patients at risk at analysis time \\t\\ is calculated as \$\$\phi_g
A(\tau - t) S_g(t) G_g(t),\$\$ where \\\phi_g\\ is the probability of
randomization to treatment group \\g\\, \\A(\tau - t)\\ is the number of
patients enrolled by calendar time \\\tau - t\\, \\S_g(t)G_g(t)\\ is the
probability of being at risk at analysis time \\t\\ for a patient in
treatment group \\g\\ after enrollment. Obviously, \\t \< \min(\tau,
T\_{\rm{fmax}})\\.

## Author

Kaifeng Lu, <kaifenglu@gmail.com>

## Examples

``` r
# Piecewise accrual, piecewise exponential survivals, and 5% dropout by
# the end of 1 year.

natrisk(t = c(9, 24), allocationRatioPlanned = 1,
        accrualTime = c(0, 3), accrualIntensity = c(10, 20),
        piecewiseSurvivalTime = c(0, 6),
        lambda1 = c(0.0533, 0.0309), lambda2 = c(0.0533, 0.0533),
        gamma1 = -log(1-0.05)/12, gamma2 = -log(1-0.05)/12,
        accrualDuration = 12, maxFollowupTime = 30, time = 30)
#>          [,1]     [,2]
#> [1,] 66.88605 62.53900
#> [2,] 16.91289 11.30083
```
