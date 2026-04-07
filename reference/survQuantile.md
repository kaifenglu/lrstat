# Brookmeyer-Crowley Confidence Interval for Quantiles of Right-Censored Time-to-Event Data

Obtains the Brookmeyer-Crowley confidence interval for quantiles of
right-censored time-to-event data.

## Usage

``` r
survQuantile(
  time,
  event,
  cilevel = 0.95,
  transform = "loglog",
  probs = as.numeric(c(0.25, 0.5, 0.75))
)
```

## Arguments

- time:

  The vector of possibly right-censored survival times.

- event:

  The vector of event indicators.

- cilevel:

  The confidence interval level. Defaults to 0.95.

- transform:

  The transformation of the survival function to use to construct the
  confidence interval. Options include "linear" (alternatively "plain"),
  "log", "loglog" (alternatively "log-log" or "cloglog"), "asinsqrt"
  (alternatively "asin" or "arcsin"), and "logit". Defaults to "loglog".

- probs:

  The vector of probabilities to calculate the quantiles. Defaults to
  c(0.25, 0.5, 0.75).

## Value

A data frame containing the estimated quantile and confidence interval
corresponding to each specified probability. It includes the following
variables:

- `prob`: The probability to calculate the quantile.

- `quantile`: The estimated quantile.

- `lower`: The lower limit of the confidence interval.

- `upper`: The upper limit of the confidence interval.

- `cilevel`: The confidence interval level.

- `transform`: The transformation of the survival function to use to
  construct the confidence interval.

## Author

Kaifeng Lu, <kaifenglu@gmail.com>

## Examples

``` r
survQuantile(
  time = c(33.7, 3.9, 10.5, 5.4, 19.5, 23.8, 7.9, 16.9, 16.6,
           33.7, 17.1, 7.9, 10.5, 38),
  event = c(0, 1, 1, 1, 1, 0, 1, 0, 0, 0, 0, 0, 1, 1),
  probs = c(0.25, 0.5, 0.75))
#>   prob quantile lower upper cilevel transform
#> 1 0.25     10.5   3.9    38    0.95    loglog
#> 2 0.50     38.0   7.9    38    0.95    loglog
#> 3 0.75     38.0  19.5    38    0.95    loglog
```
