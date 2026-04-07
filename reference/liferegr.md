# Parametric Regression Models for Failure Time Data

Obtains the parameter estimates from parametric regression models with
uncensored, right censored, left censored, or interval censored data.

## Usage

``` r
liferegr(
  data,
  stratum = "",
  time = "time",
  time2 = "",
  event = "event",
  covariates = "",
  weight = "",
  offset = "",
  id = "",
  dist = "weibull",
  init = NA_real_,
  robust = FALSE,
  plci = FALSE,
  alpha = 0.05,
  maxiter = 50,
  eps = 1e-09
)
```

## Arguments

- data:

  The input data frame that contains the following variables:

  - `stratum`: The stratum.

  - `time`: The follow-up time for right censored data, or the left end
    of each interval for interval censored data.

  - `time2`: The right end of each interval for interval censored data.

  - `event`: The event indicator, 1=event, 0=no event.

  - `covariates`: The values of baseline covariates.

  - `weight`: The weight for each observation.

  - `offset`: The offset for each observation.

  - `id`: The optional subject ID to group the score residuals in
    computing the robust sandwich variance.

- stratum:

  The name(s) of the stratum variable(s) in the input data.

- time:

  The name of the time variable or the left end of each interval for
  interval censored data in the input data.

- time2:

  The name of the right end of each interval for interval censored data
  in the input data.

- event:

  The name of the event variable in the input data for right censored
  data.

- covariates:

  The vector of names of baseline covariates in the input data.

- weight:

  The name of the weight variable in the input data.

- offset:

  The name of the offset variable in the input data.

- id:

  The name of the id variable in the input data.

- dist:

  The assumed distribution for time to event. Options include
  "exponential", "weibull", "lognormal", and "loglogistic" to be modeled
  on the log-scale, and "normal" and "logistic" to be modeled on the
  original scale.

- init:

  A vector of initial values for the model parameters, including
  regression coefficients and the log scale parameter. By default,
  initial values are derived from an intercept-only model. If this
  approach fails, ordinary least squares (OLS) estimates, ignoring
  censoring, are used instead.

- robust:

  Whether a robust sandwich variance estimate should be computed. In the
  presence of the id variable, the score residuals will be aggregated
  for each id when computing the robust sandwich variance estimate.

- plci:

  Whether to obtain profile likelihood confidence interval.

- alpha:

  The two-sided significance level.

- maxiter:

  The maximum number of iterations.

- eps:

  The tolerance to declare convergence.

## Value

A list with the following components:

- `sumstat`: The data frame of summary statistics of model fit with the
  following variables:

  - `n`: The number of observations.

  - `nevents`: The number of events.

  - `loglik0`: The log-likelihood under null.

  - `loglik1`: The maximum log-likelihood.

  - `niter`: The number of Newton-Raphson iterations.

  - `dist`: The assumed distribution.

  - `p`: The number of parameters, including the intercept, regression
    coefficients associated with the covariates, and the log scale
    parameters for the strata.

  - `nvar`: The number of regression coefficients associated with the
    covariates (excluding the intercept).

  - `robust`: Whether the robust sandwich variance estimate is
    requested.

  - `fail`: Whether the model fails to converge.

- `parest`: The data frame of parameter estimates with the following
  variables:

  - `param`: The name of the covariate for the parameter estimate.

  - `beta`: The parameter estimate.

  - `sebeta`: The standard error of parameter estimate.

  - `z`: The Wald test statistic for the parameter.

  - `expbeta`: The exponentiated parameter estimate.

  - `lower`: The lower limit of confidence interval.

  - `upper`: The upper limit of confidence interval.

  - `p`: The p-value from the chi-square test.

  - `method`: The method to compute the confidence interval and p-value.

  - `sebeta_naive`: The naive standard error of parameter estimate if
    robust variance is requested.

- `linear_predictors`: The vector of linear predictors.

- `p`: The number of parameters.

- `nvar`: The number of columns of the design matrix excluding the
  intercept.

- `param`: The parameter names.

- `beta`: The parameter estimate.

- `vbeta`: The covariance matrix for parameter estimates.

- `vbeta_naive`: The naive covariance matrix for parameter estimates.

- `terms`: The terms object.

- `xlevels`: A record of the levels of the factors used in fitting.

- `settings`: A list containing the input parameter values.

## Details

There are two ways to specify the model, one for right censored data
through the time and event variables, and the other for interval
censored data through the time (lower) and time2 (upper) variables. For
the second form, we follow the convention used in SAS PROC LIFEREG:

- If lower is not missing, upper is not missing, and lower is equal to
  upper, then there is no censoring and the event occurred at time
  lower.

- If lower is not missing, upper is not missing, and lower \< upper,
  then the event time is censored within the interval (lower, upper).

- If lower is missing, but upper is not missing, then upper will be used
  as the left censoring value.

- If lower is not missing, but upper is missing, then lower will be used
  as the right censoring value.

- If lower is not missing, upper is not missing, but lower \> upper, or
  if both lower and upper are missing, then the observation will not be
  used.

## References

John D. Kalbfleisch and Ross L. Prentice. The Statistical Analysis of
Failure Time Data. Wiley: New York, 1980.

## Author

Kaifeng Lu, <kaifenglu@gmail.com>

## Examples

``` r
library(dplyr)
#> 
#> Attaching package: ‘dplyr’
#> The following objects are masked from ‘package:stats’:
#> 
#>     filter, lag
#> The following objects are masked from ‘package:base’:
#> 
#>     intersect, setdiff, setequal, union

# right censored data
(fit1 <- liferegr(
  data = rawdata %>% filter(iterationNumber == 1) %>%
         mutate(treat = (treatmentGroup == 1)),
  stratum = "stratum",
  time = "timeUnderObservation", event = "event",
  covariates = "treat", dist = "weibull"))
#>     n nevents   loglik0   loglik1  lrchisq df     pvalue niter    dist
#> 1 491     376 -797.8817 -795.9032 3.957024  1 0.04667617     3 weibull
#> 
#>          param      coef exp(coef)   se(coef)         z             p
#> 1  (Intercept) 2.5744632 13.124270 0.08300898 31.014274 3.461104e-211
#> 2    treatTRUE 0.2376368  1.268248 0.11966101  1.985917  4.704257e-02
#> 3 Log(scale 1) 0.2860613  1.331174 0.10407996  2.748476  5.987305e-03
#> 4 Log(scale 2) 0.1150327  1.121910 0.04811086  2.390992  1.680290e-02

# tobit regression for left censored data
(fit2 <- liferegr(
  data = tobin %>% mutate(time = ifelse(durable>0, durable, NA)),
  time = "time", time2 = "durable",
  covariates = c("age", "quant"), dist = "normal"))
#>    n nevents  loglik0   loglik1  lrchisq df    pvalue niter   dist
#> 1 20       7 -29.4922 -28.94013 1.104133  2 0.5757589     4 normal
#> 
#>         param        coef    exp(coef)    se(coef)          z            p
#> 1 (Intercept) 15.14486644 3.778609e+06 16.07945315  0.9418770 3.462556e-01
#> 2         age -0.12905928 8.789219e-01  0.21858360 -0.5904344 5.548994e-01
#> 3       quant -0.04554166 9.554798e-01  0.05825412 -0.7817759 4.343463e-01
#> 4  Log(scale)  1.71785092 5.572540e+00  0.31032272  5.5356918 3.100023e-08
```
