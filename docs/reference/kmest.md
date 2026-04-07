# Kaplan-Meier Estimates of Survival Curve

Obtains the Kaplan-Meier estimates of the survival curve.

## Usage

``` r
kmest(
  data,
  stratum = "",
  time = "time",
  time2 = "",
  event = "event",
  weight = "",
  conftype = "log-log",
  conflev = 0.95,
  keep_censor = FALSE
)
```

## Arguments

- data:

  The input data frame that contains the following variables:

  - `stratum`: The stratum.

  - `time`: The follow-up time for right censored data, or the left end
    of each interval for counting process data.

  - `time2`: The right end of each interval for counting process data.
    Intervals are assumed to be open on the left and closed on the
    right, and event indicates whether an event occurred at the right
    end of each interval.

  - `event`: The event indicator, 1=event, 0=no event.

  - `weight`: The weight for each observation.

- stratum:

  The name(s) of the stratum variable(s) in the input data.

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

- conftype:

  The type of the confidence interval. One of "none", "plain", "log",
  "log-log" (the default), or "arcsin". The arcsin option bases the
  intervals on asin(sqrt(survival)).

- conflev:

  The level of the two-sided confidence interval for the survival
  probabilities. Defaults to 0.95.

- keep_censor:

  Whether to retain the censoring time in the output data frame.

## Value

A data frame with the following variables:

- `size`: The number of subjects in the stratum.

- `time`: The event time.

- `nrisk`: The number of subjects at risk.

- `nevent`: The number of subjects having the event.

- `ncensor`: The number of censored subjects.

- `surv`: The Kaplan-Meier estimate of the survival probability.

- `sesurv`: The standard error of the estimated survival probability
  based on the Greendwood formula.

- `lower`: The lower bound of confidence interval if requested.

- `upper`: The upper bound of confidence interval if requested.

- `conflev`: The level of confidence interval if requested.

- `conftype`: The type of confidence interval if requested.

- `stratum`: The stratum.

## Author

Kaifeng Lu, <kaifenglu@gmail.com>

## Examples

``` r
kmest(data = aml, stratum = "x", time = "time", event = "status")
#>    size time nrisk nevent ncensor       surv     sesurv       lower     upper
#> 1    11    9    11      1       0 0.90909091 0.08667842 0.508080206 0.9866738
#> 2    11   13    10      1       1 0.81818182 0.11629130 0.447428615 0.9511622
#> 3    11   18     8      1       0 0.71590909 0.13966497 0.350190386 0.8990240
#> 4    11   23     7      1       0 0.61363636 0.15263233 0.265752040 0.8352992
#> 5    11   31     5      1       0 0.49090909 0.16419327 0.167330910 0.7533998
#> 6    11   34     4      1       0 0.36818182 0.16266889 0.092829575 0.6570408
#> 7    11   48     2      1       0 0.18409091 0.15349275 0.011738480 0.5250148
#> 8    12    5    12      2       0 0.83333333 0.10758287 0.481714942 0.9555094
#> 9    12    8    10      2       0 0.66666667 0.13608276 0.337018933 0.8597118
#> 10   12   12     8      1       0 0.58333333 0.14231876 0.270138924 0.8009402
#> 11   12   23     6      1       0 0.48611111 0.14813006 0.191876620 0.7296716
#> 12   12   27     5      1       0 0.38888889 0.14698618 0.126272012 0.6498174
#> 13   12   30     4      1       0 0.29166667 0.13871517 0.072401609 0.5608861
#> 14   12   33     3      1       0 0.19444444 0.12187451 0.031198643 0.4614295
#> 15   12   43     2      1       0 0.09722222 0.09186636 0.005746306 0.3489039
#> 16   12   45     1      1       0 0.00000000        NaN         NaN       NaN
#>    conflev conftype x
#> 1     0.95  log-log 1
#> 2     0.95  log-log 1
#> 3     0.95  log-log 1
#> 4     0.95  log-log 1
#> 5     0.95  log-log 1
#> 6     0.95  log-log 1
#> 7     0.95  log-log 1
#> 8     0.95  log-log 2
#> 9     0.95  log-log 2
#> 10    0.95  log-log 2
#> 11    0.95  log-log 2
#> 12    0.95  log-log 2
#> 13    0.95  log-log 2
#> 14    0.95  log-log 2
#> 15    0.95  log-log 2
#> 16    0.95  log-log 2
```
