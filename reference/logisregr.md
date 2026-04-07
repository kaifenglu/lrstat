# Logistic Regression Models for Binary Data

Obtains the parameter estimates from logistic regression models with
binary data.

## Usage

``` r
logisregr(
  data,
  event = "event",
  covariates = "",
  freq = "",
  weight = "",
  offset = "",
  id = "",
  link = "logit",
  init = NA_real_,
  robust = FALSE,
  firth = FALSE,
  flic = FALSE,
  plci = FALSE,
  alpha = 0.05,
  maxiter = 50,
  eps = 1e-09
)
```

## Arguments

- data:

  The input data frame that contains the following variables:

  - `event`: The event indicator, 1=event, 0=no event.

  - `covariates`: The values of baseline covariates.

  - `freq`: The frequency for each observation.

  - `weight`: The weight for each observation.

  - `offset`: The offset for each observation.

  - `id`: The optional subject ID to group the score residuals in
    computing the robust sandwich variance.

- event:

  The name of the event variable in the input data.

- covariates:

  The vector of names of baseline covariates in the input data.

- freq:

  The name of the frequency variable in the input data. The frequencies
  must be the same for all observations within each cluster as indicated
  by the id. Thus freq is the cluster frequency.

- weight:

  The name of the weight variable in the input data.

- offset:

  The name of the offset variable in the input data.

- id:

  The name of the id variable in the input data.

- link:

  The link function linking the response probabilities to the linear
  predictors. Options include "logit" (default), "probit", and "cloglog"
  (complementary log-log).

- init:

  A vector of initial values for the model parameters. By default,
  initial values are derived from an intercept-only model.

- robust:

  Whether a robust sandwich variance estimate should be computed. In the
  presence of the id variable, the score residuals will be aggregated
  for each id when computing the robust sandwich variance estimate.

- firth:

  Whether the firth's bias reducing penalized likelihood should be used.
  The default is `FALSE`.

- flic:

  Whether to apply intercept correction to obtain more accurate
  predicted probabilities. The default is `FALSE`.

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

  - `n`: The number of subjects.

  - `nevents`: The number of events.

  - `loglik0`: The (penalized) log-likelihood under null.

  - `loglik1`: The maximum (penalized) log-likelihood.

  - `niter`: The number of Newton-Raphson iterations.

  - `p`: The number of parameters, including the intercept, and
    regression coefficients associated with the covariates.

  - `link`: The link function.

  - `robust`: Whether a robust sandwich variance estimate should be
    computed.

  - `firth`: Whether the firth's penalized likelihood is used.

  - `flic`: Whether to apply intercept correction.

  - `fail`: Whether the model fails to converge.

  - `loglik0_unpenalized`: The unpenalized log-likelihood under null.

  - `loglik1_unpenalized`: The maximum unpenalized log-likelihood.

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

  - `sebeta_naive`: The naive standard error of parameter estimate.

- `fitted`: The data frame with the following variables:

  - `linear_predictors`: The linear fit on the link function scale.

  - `fitted_values`: The fitted probabilities of having an event,
    obtained by transforming the linear predictors by the inverse of the
    link function.

- `p`: The number of parameters.

- `link`: The link function.

- `param`: The parameter names.

- `beta`: The parameter estimate.

- `vbeta`: The covariance matrix for parameter estimates.

- `vbeta_naive`: The naive covariance matrix for parameter estimates.

- `linear_predictors`: The linear fit on the link function scale.

- `fitted_values`: The fitted probabilities of having an event.

- `terms`: The terms object.

- `xlevels`: A record of the levels of the factors used in fitting.

- `settings`: A list containing the input parameter values.

## Details

Fitting a logistic regression model using Firth's bias reduction method
is equivalent to penalization of the log-likelihood by the Jeffreys
prior. Firth's penalized log-likelihood is given by \$\$l(\beta) +
\frac{1}{2} \log(\mbox{det}(I(\beta)))\$\$ and the components of the
gradient \\g(\beta)\\ are computed as \$\$g(\beta_j) + \frac{1}{2}
\mbox{trace}\left(I(\beta)^{-1} \frac{\partial I(\beta)}{\partial
\beta_j}\right)\$\$ The Hessian matrix is not modified by this penalty.

Firth's method reduces bias in maximum likelihood estimates of
coefficients, but it introduces a bias toward one-half in the predicted
probabilities.

A straightforward modification to Firth’s logistic regression to achieve
unbiased average predicted probabilities involves a post hoc adjustment
of the intercept. This approach, known as Firth’s logistic regression
with intercept correction (FLIC), preserves the bias-corrected effect
estimates. By excluding the intercept from penalization, it ensures that
we don't sacrifice the accuracy of effect estimates to improve the
predictions.

## References

David Firth. Bias Reduction of Maximum Likelihood Estimates. Biometrika
1993; 80:27–38.

Georg Heinze and Michael Schemper. A solution to the problem of
separation in logistic regression. Statistics in Medicine
2002;21:2409–2419.

Rainer Puhr, Georg Heinze, Mariana Nold, Lara Lusa, and Angelika
Geroldinger. Firth's logistic regression with rare events: accurate
effect estimates and predictions? Statistics in Medicine 2017;
36:2302-2317.

## Author

Kaifeng Lu, <kaifenglu@gmail.com>

## Examples

``` r
(fit1 <- logisregr(
  ingots, event = "NotReady", covariates = "Heat*Soak", freq = "Freq"))
#>     n nevents   loglik0   loglik1  lrchisq df      pvalue niter  link firth
#> 1 387      12 -53.49422 -47.61109 11.76625  3 0.008228328     5 logit FALSE
#>    flic
#> 1 FALSE
#> 
#>         param         coef   exp(coef)   se(coef)          z            p
#> 1 (Intercept) -5.990190966 0.002503186 1.66662245 -3.5942099 0.0003253774
#> 2        Heat  0.096338891 1.101132164 0.04706696  2.0468475 0.0406730611
#> 3        Soak  0.299574022 1.349283920 0.75506786  0.3967511 0.6915509968
#> 4   Heat:Soak -0.008839772 0.991199184 0.02531929 -0.3491319 0.7269902725
```
