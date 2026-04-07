# Survival Curve for Proportional Hazards Regression Models

Obtains the predicted survivor function for a proportional hazards
regression model.

## Usage

``` r
survfit_phregr(
  object,
  newdata,
  sefit = TRUE,
  conftype = "log-log",
  conflev = 0.95
)
```

## Arguments

- object:

  The output from the `phregr` call.

- newdata:

  A data frame with the same variable names as those that appear in the
  `phregr` call. For right-censored data, one curve is produced per row
  to represent a cohort whose covariates correspond to the values in
  `newdata`. For counting-process data, one curve is produced per `id`
  in `newdata` to present the survival curve along the path of
  time-dependent covariates at the observed event times in the data used
  to fit `phregr`.

- sefit:

  Whether to compute the standard error of the survival estimates.

- conftype:

  The type of the confidence interval. One of `"none"`, `"plain"`,
  `"log"`, `"log-log"` (the default), or `"arcsin"`. The `arcsin` option
  bases the intervals on `asin(sqrt(surv))`.

- conflev:

  The level of the two-sided confidence interval for the survival
  probabilities. Defaults to 0.95.

## Value

A data frame with the following variables:

- `id`: The id of the subject for counting-process data with
  time-dependent covariates.

- `time`: The observed times in the data used to fit `phregr`.

- `nrisk`: The number of patients at risk at the time point in the data
  used to fit `phregr`.

- `nevent`: The number of patients having event at the time point in the
  data used to fit `phregr`.

- `cumhaz`: The cumulative hazard at the time point.

- `surv`: The estimated survival probability at the time point.

- `sesurv`: The standard error of the estimated survival probability.

- `lower`: The lower confidence limit for survival probability.

- `upper`: The upper confidence limit for survival probability.

- `conflev`: The level of the two-sided confidence interval.

- `conftype`: The type of the confidence interval.

- `covariates`: The values of covariates based on `newdata`.

- `stratum`: The stratum of the subject.

## Details

If `newdata` is not provided and there is no covariate, survival curves
based on the `basehaz` data frame will be produced.

## References

Terry M. Therneau and Patricia M. Grambsch. Modeling Survival Data:
Extending the Cox Model. Springer-Verlag, 2000.

## Author

Kaifeng Lu, <kaifenglu@gmail.com>

## Examples

``` r
library(dplyr)

# Example 1 with right-censored data
fit1 <- phregr(data = rawdata %>% filter(iterationNumber == 1) %>%
                 mutate(treat = 1*(treatmentGroup == 1)),
               stratum = "stratum",
               time = "timeUnderObservation", event = "event",
               covariates = "treat")

surv1 <- survfit_phregr(fit1,
                        newdata = data.frame(
                          stratum = as.integer(c(1,1,2,2)),
                          treat = c(1,0,1,0)))
head(surv1)
#>         time nrisk nevent ncensor      cumhaz      surv      sesurv     lower
#> 1 0.08705128   103      1       0 0.008787705 0.9912508 0.008724575 0.9393399
#> 2 0.20111510   102      1       0 0.017653319 0.9825016 0.012303519 0.9315562
#> 3 0.31332543   101      1       0 0.026598234 0.9737524 0.015025815 0.9204223
#> 4 0.36024791   100      1       0 0.035623884 0.9650032 0.017300688 0.9088955
#> 5 0.48874656    99      1       0 0.044731738 0.9562540 0.019287111 0.8974048
#> 6 0.72067837    98      1       0 0.053941225 0.9474878 0.021072360 0.8860117
#>       upper conflev conftype treat stratum
#> 1 0.9987667    0.95  log-log     1       1
#> 2 0.9956141    0.95  log-log     1       1
#> 3 0.9915047    0.95  log-log     1       1
#> 4 0.9868028    0.95  log-log     1       1
#> 5 0.9816852    0.95  log-log     1       1
#> 6 0.9762450    0.95  log-log     1       1

# Example 2 with counting process data and robust variance estimate
fit2 <- phregr(data = heart %>% mutate(rx = as.numeric(transplant) - 1),
               time = "start", time2 = "stop", event = "event",
               covariates = c("rx", "age"), id = "id", robust = TRUE)

surv2 <- survfit_phregr(fit2,
                        newdata = data.frame(
                          id = c(4,4,11,11),
                          age = c(-7.737,-7.737,-0.019,-0.019),
                          start = c(0,36,0,26),
                          stop = c(36,39,26,153),
                          rx = c(0,1,0,1)))
head(surv2)
#>   time nrisk nevent ncensor      cumhaz      surv      sesurv     lower
#> 1  1.0   103      1       2 0.008019895 0.9920122 0.008005765 0.9439912
#> 2  2.0   102      3       3 0.032646621 0.9678805 0.016181660 0.9147786
#> 3  3.0    99      3       3 0.058091568 0.9435635 0.021576395 0.8819199
#> 4  4.0    96      0       2 0.058091568 0.9435635 0.021576395 0.8819199
#> 5  4.5    96      0       1 0.058091568 0.9435635 0.021576395 0.8819199
#> 6  5.0    96      2       2 0.075508031 0.9272723 0.024572348 0.8605201
#>       upper conflev conftype rx    age id
#> 1 0.9988847    0.95  log-log  0 -7.737  4
#> 2 0.9881058    0.95  log-log  0 -7.737  4
#> 3 0.9735009    0.95  log-log  0 -7.737  4
#> 4 0.9735009    0.95  log-log  0 -7.737  4
#> 5 0.9735009    0.95  log-log  0 -7.737  4
#> 6 0.9627567    0.95  log-log  0 -7.737  4
```
