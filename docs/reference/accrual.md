# Number of Enrolled Subjects

Obtains the number of subjects enrolled by given calendar times.

## Usage

``` r
accrual(time, accrualTime, accrualIntensity, accrualDuration)
```

## Arguments

- time:

  A vector of calendar times at which to calculate the number of
  enrolled subjects.

- accrualTime:

  A vector that specifies the starting time of piecewise Poisson
  enrollment time intervals. Must start with 0, e.g., `c(0, 3)` breaks
  the time axis into 2 accrual intervals: \\\[0, 3)\\ and \\\[3,
  \infty)\\.

- accrualIntensity:

  A vector of accrual intensities. One for each accrual time interval.

- accrualDuration:

  Duration of the enrollment period.

## Value

A vector of total number of subjects enrolled by the specified calendar
times.

## Author

Kaifeng Lu, <kaifenglu@gmail.com>

## Examples

``` r
# Example 1: Uniform enrollment with 20 patients per month for 12 months.

accrual(time = 3, accrualTime = 0, accrualIntensity = 20,
        accrualDuration = 12)
#> [1] 60


# Example 2: Piecewise accrual, 10 patients per month for the first
# 3 months, and 20 patients per month thereafter. Patient recruitment
# ends at 12 months for the study.

accrual(time = c(2, 9), accrualTime = c(0, 3),
        accrualIntensity = c(10, 20), accrualDuration = 12)
#> [1]  20 150
```
