# Proportional Hazards Regression Models

Obtains the hazard ratio estimates from the proportional hazards
regression model with right censored or counting process data.

## Usage

``` r
phregr(
  data,
  stratum = "",
  time = "time",
  time2 = "",
  event = "event",
  covariates = "",
  weight = "",
  offset = "",
  id = "",
  ties = "efron",
  init = NA_real_,
  robust = FALSE,
  est_basehaz = TRUE,
  est_resid = TRUE,
  firth = FALSE,
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
    of each interval for counting process data.

  - `time2`: The right end of each interval for counting process data.
    Intervals are assumed to be open on the left and closed on the
    right, and event indicates whether an event occurred at the right
    end of each interval.

  - `event`: The event indicator, 1=event, 0=no event.

  - `covariates`: The values of baseline covariates (and time-dependent
    covariates in each interval for counting process data).

  - `weight`: The weight for each observation.

  - `offset`: The offset for each observation.

  - `id`: The optional subject ID for counting process data with
    time-dependent covariates.

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

- covariates:

  The vector of names of baseline and time-dependent covariates in the
  input data.

- weight:

  The name of the weight variable in the input data.

- offset:

  The name of the offset variable in the input data.

- id:

  The name of the id variable in the input data.

- ties:

  The method for handling ties, either "breslow" or "efron" (default).

- init:

  The vector of initial values. Defaults to zero for all variables.

- robust:

  Whether a robust sandwich variance estimate should be computed. In the
  presence of the id variable, the score residuals will be aggregated
  for each id when computing the robust sandwich variance estimate.

- est_basehaz:

  Whether to estimate the baseline hazards. Defaults to `TRUE`.

- est_resid:

  Whether to estimate the martingale residuals. Defaults to `TRUE`.

- firth:

  Whether to use Firth’s penalized likelihood method. Defaults to
  `FALSE`.

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

  - `loglik0`: The (penalized) log-likelihood under null.

  - `loglik1`: The maximum (penalized) log-likelihood.

  - `scoretest`: The score test statistic.

  - `niter`: The number of Newton-Raphson iterations.

  - `ties`: The method for handling ties, either "breslow" or "efron".

  - `p`: The number of columns of the Cox model design matrix.

  - `robust`: Whether to use the robust variance estimate.

  - `firth`: Whether to use Firth's penalized likelihood method.

  - `fail`: Whether the model fails to converge.

  - `loglik0_unpenalized`: The unpenalized log-likelihood under null.

  - `loglik1_unpenalized`: The maximum unpenalized log-likelihood.

- `parest`: The data frame of parameter estimates with the following
  variables:

  - `param`: The name of the covariate for the parameter estimate.

  - `beta`: The log hazard ratio estimate.

  - `sebeta`: The standard error of log hazard ratio estimate.

  - `z`: The Wald test statistic for log hazard ratio.

  - `expbeta`: The hazard ratio estimate.

  - `lower`: The lower limit of confidence interval.

  - `upper`: The upper limit of confidence interval.

  - `p`: The p-value from the chi-square test.

  - `method`: The method to compute the confidence interval and p-value.

  - `sebeta_naive`: The naive standard error of log hazard ratio
    estimate if robust variance is requested.

- `basehaz`: The data frame of baseline hazards with the following
  variables (if est_basehaz is TRUE):

  - `time`: The observed event time.

  - `nrisk`: The number of patients at risk at the time point.

  - `nevent`: The number of events at the time point.

  - `haz`: The baseline hazard at the time point.

  - `varhaz`: The variance of the baseline hazard at the time point
    assuming the parameter beta is known.

  - `gradhaz`: The gradient of the baseline hazard with respect to beta
    at the time point (in the presence of covariates).

  - `stratum`: The stratum.

- `residuals`: The martingale residuals.

- `linear_predictors`: The vector of linear predictors.

- `p`: The number of parameters.

- `param`: The parameter names.

- `beta`: The parameter estimate.

- `vbeta`: The covariance matrix for parameter estimates.

- `vbeta_naive`: The naive covariance matrix for parameter estimates.

- `terms`: The terms object.

- `xlevels`: A record of the levels of the factors used in fitting.

- `settings`: A list containing the input parameter values.

## References

Per K. Anderson and Richard D. Gill. Cox's regression model for counting
processes, a large sample study. Annals of Statistics 1982;
10:1100-1120.

Terry M. Therneau and Patricia M. Grambsch. Modeling Survival Data:
Extending the Cox Model. Springer-Verlag, 2000.

## Author

Kaifeng Lu, <kaifenglu@gmail.com>

## Examples

``` r
library(dplyr)

# Example 1 with right-censored data
(fit1 <- phregr(
  data = rawdata %>% filter(iterationNumber == 1) %>%
    mutate(treat = 1*(treatmentGroup == 1)),
  stratum = "stratum",
  time = "timeUnderObservation", event = "event",
  covariates = "treat", est_basehaz = FALSE, est_resid = FALSE))
#>     n nevents   loglik0   loglik1  lrchisq df    pvalue scoretest niter  ties
#> 1 491     376 -1893.975 -1892.253 3.444581  1 0.0634595  3.449958     2 efron
#> 
#>   param       coef exp(coef)  se(coef)         z          p
#> 1 treat -0.1920025 0.8253048 0.1035279 -1.854596 0.06365391

# Example 2 with counting process data and robust variance estimate
(fit2 <- phregr(
  data = heart %>% mutate(rx = as.numeric(transplant) - 1),
  time = "start", time2 = "stop", event = "event",
  covariates = c("rx", "age"), id = "id",
  robust = TRUE, est_basehaz = TRUE, est_resid = TRUE))
#>     n nevents   loglik0   loglik1  lrchisq df     pvalue scoretest niter  ties
#> 1 172      75 -298.1214 -295.5367 5.169366  2 0.07541998  4.641254     3 efron
#> 
#>   param         coef exp(coef)   se(coef)  robust se           z          p
#> 1    rx -0.004178247 0.9958305 0.31208501 0.31415106 -0.01330012 0.98938835
#> 2   age  0.030742256 1.0312197 0.01449603 0.01518176  2.02494619 0.04287289
```
