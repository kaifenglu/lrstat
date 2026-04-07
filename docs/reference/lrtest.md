# Log-Rank Test of Survival Curve Difference

Obtains the log-rank test using the Fleming-Harrington family of
weights.

## Usage

``` r
lrtest(
  data,
  stratum = "",
  treat = "treat",
  time = "time",
  time2 = "",
  event = "event",
  weight = "",
  weight_readj = FALSE,
  rho1 = 0,
  rho2 = 0
)
```

## Arguments

- data:

  The input data frame or list of data frames that contains the
  following variables:

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

  The name(s) of the stratum variable(s) in the input data.

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

- weight_readj:

  Whether the weight variable at each event time will be readjusted to
  be proportional to the number at risk by treatment group. Defaults to
  `FALSE`.

- rho1:

  The first parameter of the Fleming-Harrington family of weighted
  log-rank test. Defaults to 0 for conventional log-rank test.

- rho2:

  The second parameter of the Fleming-Harrington family of weighted
  log-rank test. Defaults to 0 for conventional log-rank test.

## Value

A data frame with the following variables:

- `uscore`: The numerator of the log-rank test statistic.

- `vscore`: The variance of the log-rank score test statistic.

- `logRankZ`: The Z-statistic value.

- `logRankPValue`: The two-sided p-value.

- `weight_readj`: Whether the weight variable will be readjusted.

- `rho1`: The first parameter of the Fleming-Harrington weights.

- `rho2`: The second parameter of the Fleming-Harrington weights.

## Author

Kaifeng Lu, <kaifenglu@gmail.com>

## Examples

``` r
lrtest(rawdata[rawdata$iterationNumber == 1, ],
       stratum = "stratum", treat = "treatmentGroup",
       time = "timeUnderObservation", event = "event",
       rho1 = 0.5, rho2 = 0)
#>      uscore   vscore  logRankZ logRankPValue weight_readj rho1 rho2
#> 1 -12.10312 56.85576 -1.605129     0.1084654        FALSE  0.5    0
```
