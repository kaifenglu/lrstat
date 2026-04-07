# Number of Subjects Having an Event by Calendar Time

Obtains the number of subjects having an event by given calendar times
for each treatment group.

## Usage

``` r
nevent(
  time = NA_real_,
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

- time:

  A vector of calendar times at which to calculate the number of
  patients having an event.

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

A matrix of the number of patients having an event at the specified
calendar times (row) for each treatment group (column).

## Details

For a given treatment group \\g\\ and calendar time \\\tau\\, the number
of patients having an event by calendar time \\\tau\\ is calculated as
\\I_1 + I_2\\, where \$\$I_1 = \phi_g A(\tau - T\_{\rm{fmax}})
P_g(T\_{\rm{fmax}}),\$\$ and \$\$I_2 = \phi_g \int\_{\tau -
T\_{\rm{fmax}}}^{\tau} a(u) P_g(\tau - u) du,\$\$ where \\\phi_g\\ is
the probability of randomization to treatment group \\g\\, \\A(\tau -
T\_{\rm{fmax}})\\ is the number of patients enrolled by calendar time
\\\tau - T\_{\rm{fmax}}\\, \\P_g(T\_{\rm{fmax}})\\ is the probability of
having an event by the maximum follow-up time \\T\_{\rm{fmax}}\\ for a
patient in treatment group \\g\\ after enrollment, \\a(u)\\ is the
accrual intensity at calendar time \\u\\, and \\P_g(\tau - u)\\ is the
probability of having an event by calendar time \\\tau\\ for a patient
in treatment group \\g\\ enrolled at calendar time \\u\\.

## Author

Kaifeng Lu, <kaifenglu@gmail.com>

## Examples

``` r
# Piecewise accrual, piecewise exponential survivals, and 5% dropout by
# the end of 1 year.
nevent(time = c(9, 24), allocationRatioPlanned = 1,
       accrualTime = c(0, 3), accrualIntensity = c(10, 20),
       piecewiseSurvivalTime = c(0, 6),
       lambda1 = c(0.0533, 0.0309), lambda2 = c(0.0533, 0.0533),
       gamma1 = -log(1-0.05)/12, gamma2 = -log(1-0.05)/12,
       accrualDuration = 12, maxFollowupTime = 30)
#>          [,1]     [,2]
#> [1,] 13.10990 13.43671
#> [2,] 49.61969 60.81305
```
