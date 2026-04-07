# Estimate of Restricted Mean Survival Time Difference

Obtains the estimate of restricted mean survival time difference between
two treatment groups.

## Usage

``` r
rmdiff(
  data,
  stratum = "",
  treat = "treat",
  time = "time",
  event = "event",
  milestone = 0,
  rmstDiffH0 = 0,
  conflev = 0.95,
  biascorrection = FALSE
)
```

## Arguments

- data:

  The input data frame that contains the following variables:

  - `stratum`: The stratum.

  - `treat`: The treatment.

  - `time`: The possibly right-censored survival time.

  - `event`: The event indicator.

- stratum:

  The name of the stratum variable in the input data.

- treat:

  The name of the treatment variable in the input data.

- time:

  The name of the time variable in the input data.

- event:

  The name of the event variable in the input data.

- milestone:

  The milestone time at which to calculate the restricted mean survival
  time.

- rmstDiffH0:

  The difference in restricted mean survival times under the null
  hypothesis. Defaults to 0 for superiority test.

- conflev:

  The level of the two-sided confidence interval for the difference in
  restricted mean survival times. Defaults to 0.95.

- biascorrection:

  Whether to apply bias correction for the variance estimate of
  individual restricted mean survival times. Defaults to no bias
  correction.

## Value

A data frame with the following variables:

- `milestone`: The milestone time relative to randomization.

- `rmstDiffH0`: The difference in restricted mean survival times under
  the null hypothesis.

- `rmst1`: The estimated restricted mean survival time for the treatment
  group.

- `rmst2`: The estimated restricted mean survival time for the control
  group.

- `rmstDiff`: The estimated difference in restricted mean survival
  times.

- `vrmst1`: The variance for rmst1.

- `vrmst2`: The variance for rmst2.

- `sermstDiff`: The standard error for rmstDiff.

- `rmstDiffZ`: The Z-statistic value.

- `rmstDiffPValue`: The two-sided p-value.

- `lower`: The lower bound of confidence interval.

- `upper`: The upper bound of confidence interval.

- `conflev`: The level of confidence interval.

- `biascorrection`: Whether to apply bias correction for the variance
  estimate of individual restricted mean survival times.

## Author

Kaifeng Lu, <kaifenglu@gmail.com>

## Examples

``` r
rmdiff(data = rawdata[rawdata$iterationNumber == 1, ],
       stratum = "stratum", treat = "treatmentGroup",
       time = "timeUnderObservation", event = "event",
       milestone = 12)
#>   milestone rmstDiffH0    rmst1    rmst2  rmstDiff     vrmst1     vrmst2
#> 1        12          0 7.974701 7.595799 0.3789028 0.08235337 0.07767696
#>   sermstDiff rmstDiffZ rmstDiffPValue      lower    upper conflev
#> 1  0.4000379 0.9471671      0.3435536 -0.4051572 1.162963    0.95
#>   biascorrection
#> 1          FALSE
```
