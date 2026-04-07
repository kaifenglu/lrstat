# Estimate of Restricted Mean Survival Time

Obtains the estimate of restricted means survival time for each stratum.

## Usage

``` r
rmest(
  data,
  stratum = "",
  time = "time",
  event = "event",
  milestone = 0,
  conflev = 0.95,
  biascorrection = FALSE
)
```

## Arguments

- data:

  The input data frame that contains the following variables:

  - `stratum`: The stratum.

  - `time`: The possibly right-censored survival time.

  - `event`: The event indicator.

- stratum:

  The name of the stratum variable in the input data.

- time:

  The name of the time variable in the input data.

- event:

  The name of the event variable in the input data.

- milestone:

  The milestone time at which to calculate the restricted mean survival
  time.

- conflev:

  The level of the two-sided confidence interval for the survival
  probabilities. Defaults to 0.95.

- biascorrection:

  Whether to apply bias correction for the variance estimate. Defaults
  to no bias correction.

## Value

A data frame with the following variables:

- `stratum`: The stratum variable.

- `size`: The number of subjects in the stratum.

- `milestone`: The milestone time relative to randomization.

- `rmst`: The estimate of restricted mean survival time.

- `stderr`: The standard error of the estimated rmst.

- `lower`: The lower bound of confidence interval if requested.

- `upper`: The upper bound of confidence interval if requested.

- `conflev`: The level of confidence interval if requested.

- `biascorrection`: Whether to apply bias correction for the variance
  estimate.

## Author

Kaifeng Lu, <kaifenglu@gmail.com>

## Examples

``` r
rmest(data = aml, stratum = "x",
      time = "time", event = "status", milestone = 24)
#>   size milestone     rmst   stderr    lower    upper conflev biascorrection x
#> 1   11        24 20.92045 1.541833 17.89852 23.94239    0.95          FALSE 1
#> 2   12        24 17.06944 2.361346 12.44129 21.69760    0.95          FALSE 2
```
