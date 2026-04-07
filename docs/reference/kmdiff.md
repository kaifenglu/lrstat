# Estimate of Milestone Survival Difference

Obtains the estimate of milestone survival difference between two
treatment groups.

## Usage

``` r
kmdiff(
  data,
  stratum = "",
  treat = "treat",
  time = "time",
  time2 = "",
  event = "event",
  weight = "",
  milestone = 0,
  survDiffH0 = 0,
  conflev = 0.95
)
```

## Arguments

- data:

  The input data frame that contains the following variables:

  - `stratum`: The stratum.

  - `treat`: The treatment.

  - `time`: The follow-up time for right censored data, or the left end
    of each interval for counting process data.

  - `time2`: The right end of each interval for counting process data.
    Intervals are assumed to be open on the left and closed on the
    right, and event indicates whether an event occurred at the right
    end of each interval.

  - `event`: The event indicator, 1=event, 0=no event.

  - `weight`: The weight for each observation.

- stratum:

  The name of the stratum variable in the input data.

- treat:

  The name of the treatment variable in the input data.

- time:

  The name of the time variable or the left end of each interval for
  counting process data in the input data.

- time2:

  The name of the right end of each interval for counting process data
  in the input data.

- event:

  The name of the event variable in the input data.

- weight:

  The name of the weight variable in the input data.

- milestone:

  The milestone time at which to calculate the survival probability.

- survDiffH0:

  The difference in milestone survival probabilities under the null
  hypothesis. Defaults to 0 for superiority test.

- conflev:

  The level of the two-sided confidence interval for the difference in
  milestone survival probabilities. Defaults to 0.95.

## Value

A data frame with the following variables:

- `milestone`: The milestone time relative to randomization.

- `survDiffH0`: The difference in milestone survival probabilities under
  the null hypothesis.

- `surv1`: The estimated milestone survival probability for the
  treatment group.

- `surv2`: The estimated milestone survival probability for the control
  group.

- `survDiff`: The estimated difference in milestone survival
  probabilities.

- `vsurv1`: The variance for surv1.

- `vsurv2`: The variance for surv2.

- `sesurvDiff`: The standard error for survDiff.

- `survDiffZ`: The Z-statistic value.

- `survDiffPValue`: The two-sided p-value.

- `lower`: The lower bound of confidence interval.

- `upper`: The upper bound of confidence interval.

- `conflev`: The level of confidence interval.

## Author

Kaifeng Lu, <kaifenglu@gmail.com>

## Examples

``` r
kmdiff(data = rawdata[rawdata$iterationNumber == 1, ],
       stratum = "stratum", treat = "treatmentGroup",
       time = "timeUnderObservation", event = "event",
       milestone = 12)
#>   milestone survDiffH0     surv1     surv2   survDiff      vsurv1      vsurv2
#> 1        12          0 0.4848425 0.3859297 0.09891287 0.001001194 0.001004362
#>   sesurvDiff survDiffZ survDiffPValue      lower     upper conflev
#> 1 0.04478344  2.208693     0.02719602 0.01113893 0.1866868    0.95
```
