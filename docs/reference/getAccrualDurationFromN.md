# Accrual Duration to Enroll Target Number of Subjects

Obtains the accrual duration to enroll the target number of subjects.

## Usage

``` r
getAccrualDurationFromN(nsubjects, accrualTime, accrualIntensity)
```

## Arguments

- nsubjects:

  The vector of target number of subjects.

- accrualTime:

  A vector that specifies the starting time of piecewise Poisson
  enrollment time intervals. Must start with 0, e.g., `c(0, 3)` breaks
  the time axis into 2 accrual intervals: \\\[0, 3)\\ and \\\[3,
  \infty)\\.

- accrualIntensity:

  A vector of accrual intensities. One for each accrual time interval.

## Value

A vector of accrual durations.

## Author

Kaifeng Lu, <kaifenglu@gmail.com>

## Examples

``` r
getAccrualDurationFromN(nsubjects = c(20, 150), accrualTime = c(0, 3),
                        accrualIntensity = c(10, 20))
#> [1] 2 9
```
