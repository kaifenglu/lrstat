# Confidence Interval After Trial Termination

Obtains the p-value, median unbiased point estimate, and confidence
interval after the end of a group sequential trial.

## Usage

``` r
getCI_mams(
  M = NA_integer_,
  r = 1,
  corr_known = TRUE,
  L = NA_integer_,
  zL = NA_real_,
  IMax = NA_real_,
  informationRates = NA_real_,
  efficacyStopping = NA_integer_,
  criticalValues = NULL,
  alpha = 0.025,
  typeAlphaSpending = "sfOF",
  parameterAlphaSpending = NA_real_,
  spendingTime = NA_real_
)
```

## Arguments

- M:

  Number of active treatment arms.

- r:

  Randomization ratio of each active arm to the common control.

- corr_known:

  Logical. If `TRUE`, the correlation between Wald statistics is derived
  from the randomization ratio `r` as \\r / (r + 1)\\. If `FALSE`, a
  conservative correlation of 0 is used.

- L:

  The termination look.

- zL:

  The vector of z-test statistics at the termination look.

- IMax:

  Maximum information for any active arm versus the common control.

- informationRates:

  The information rates up to look `L`.

- efficacyStopping:

  Indicators of whether efficacy stopping is allowed at each stage up to
  look `L`. Defaults to true if left unspecified.

- criticalValues:

  The matrix of by-level upper boundaries on the z-test statistic scale
  for efficacy stopping up to look `L`. The first column is for level
  `M`, the second column is for level `M - 1`, and so on, with the last
  column for level 1. If left unspecified, the critical values will be
  computed based on the specified alpha spending function.

- alpha:

  The significance level. Defaults to 0.025.

- typeAlphaSpending:

  The type of alpha spending. One of the following: `"OF"` for
  O'Brien-Fleming boundaries, `"P"` for Pocock boundaries, `"WT"` for
  Wang & Tsiatis boundaries, `"sfOF"` for O'Brien-Fleming type spending
  function, `"sfP"` for Pocock type spending function, `"sfKD"` for Kim
  & DeMets spending function, `"sfHSD"` for Hwang, Shi & DeCani spending
  function, `"user"` for user defined spending, and `"none"` for no
  early efficacy stopping. Defaults to `"sfOF"`.

- parameterAlphaSpending:

  The parameter value for the alpha spending. Corresponds to \\\Delta\\
  for `"WT"`, \\\rho\\ for `"sfKD"`, and \\\gamma\\ for `"sfHSD"`.

- spendingTime:

  The error spending time up to look `L`. Defaults to missing, in which
  case, it is the same as `informationRates`.

## Value

A data frame with the following components:

- `level`: Number of individual hypotheses considered for multiplicity.

- `index`: The index of the treatment arm among the M active arms.

- `pvalue`: p-value for rejecting the null hypothesis.

- `thetahat`: Median unbiased point estimate of the parameter.

- `cilevel`: Confidence interval level.

- `lower`: Lower bound of confidence interval.

- `upper`: Upper bound of confidence interval.

## References

Ping Gao, Yingqiu Li. Adaptive multiple comparison sequential design
(AMCSD) for clinical trials. Journal of Biopharmaceutical Statistics,
2024, 34(3), 424-440.

## Author

Kaifeng Lu, <kaifenglu@gmail.com>

## Examples

``` r
getCI_mams(
  L = 2, zL = c(2.075, 2.264),
  M = 2, r = 1, corr_known = FALSE,
  IMax = 300 / 4, informationRates = c(1/2, 1),
  alpha = 0.025, typeAlphaSpending = "sfOF")
#>   level index     pvalue  thetahat cilevel       lower     upper
#> 1     2     2 0.02407215 0.1979009    0.95 0.001703644 0.4874185
#> 2     1     1 0.01957636 0.2389416    0.95 0.011919361 0.4657280
```
