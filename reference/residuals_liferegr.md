# Residuals for Parametric Regression Models for Failure Time Data

Obtains the response, martingale, deviance, dfbeta, and likelihood
displacement residuals for a parametric regression model for failure
time data.

## Usage

``` r
residuals_liferegr(
  object,
  type = c("response", "martingale", "deviance", "dfbeta", "dfbetas", "working",
    "ldcase", "ldresp", "ldshape", "matrix"),
  collapse = FALSE,
  weighted = (type %in% c("dfbeta", "dfbetas"))
)
```

## Arguments

- object:

  The output from the `phregr` call.

- type:

  The type of residuals desired, with options including `"response"`,
  `"martingale"`, `"deviance"`, `"dfbeta"`, `"dfbetas"`, `"working"`,
  `"ldcase"`, `"ldresp"`, `"ldshape"`, and `"matrix"`.

- collapse:

  Whether to collapse the residuals by `id`.

- weighted:

  Whether to compute weighted residuals.

## Value

Either a vector or a matrix of residuals, depending on the specified
type:

- `response` residuals are on the scale of the original data.

- `martingale` residuals are event indicators minus the cumulative
  hazards for event or right-censored data.

- `working` residuals are on the scale of the linear predictor.

- `deviance` residuals are on the log-likelihood scale.

- `dfbeta` residuals are returned as a matrix, where the \\i\\-th row
  represents the approximate change in the model coefficients resulting
  from the inclusion of subject \\i\\.

- `dfbetas` residuals are similar to `dfbeta` residuals, but each column
  is scaled by the standard deviation of the corresponding coefficient.

- `matrix` residuals are a matrix of derivatives of the log-likelihood
  function. Let \\L\\ be the log-likelihood, \\p\\ be the linear
  predictor (\\X\beta\\), and \\s\\ be \\log(\sigma)\\. Then the
  resulting matrix contains six columns: \\L\\, \\\partial L/\partial
  p\\, \\\partial^2 L/\partial p^2\\, \\\partial L/\partial s\\,
  \\\partial^2 L/\partial s^2\\, and \\\partial L^2/\partial p\partial
  s\\.

- `ldcase` residulas are likelihood displacement for case weight
  perturbation.

- `ldresp` residuals are likelihood displacement for response value
  perturbation.

- `ldshape` residuals are likelihood displacement related to the shape
  parameter.

## Details

The algorithms follow the `residuals.survreg` function in the `survival`
package, except for martingale residuals, which are defined only for
event or right-censored data for exponential, weibull, lognormal, and
loglogistic distributions.

## References

Escobar, L. A. and Meeker, W. Q. Assessing influence in regression
analysis with censored data. Biometrics 1992; 48:507-528.

## Author

Kaifeng Lu, <kaifenglu@gmail.com>

## Examples

``` r
library(dplyr)

fit1 <- liferegr(
  data = tobin %>% mutate(time = ifelse(durable>0, durable, NA)),
  time = "time", time2 = "durable",
  covariates = c("age", "quant"), dist = "normal")

resid <- residuals_liferegr(fit1, type = "response")
head(resid)
#> [1] 3.0496868 5.0125418 0.5416331 0.2560716 1.8501772 2.4098780
```
