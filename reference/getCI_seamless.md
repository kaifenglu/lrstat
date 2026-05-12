# Confidence Interval After Trial Termination for Phase 2/3 Seamless Design

Obtains the p-value, point estimate, and confidence interval after the
end of a phase 2/3 seamless trial.

## Usage

``` r
getCI_seamless(
  M = NA_integer_,
  r = 1,
  corr_known = TRUE,
  L = NA_integer_,
  zL = NA_real_,
  IMax = NA_real_,
  informationRates = NA_real_,
  efficacyStopping = NA_integer_,
  criticalValues = NA_real_,
  alpha = 0.025,
  typeAlphaSpending = "sfOF",
  parameterAlphaSpending = NA_real_,
  spendingTime = NA_real_
)
```

## Arguments

- M:

  Number of active treatment arms in Phase 2.

- r:

  Randomization ratio of each active arm to the common control in Phase
  2.

- corr_known:

  Logical. If `TRUE`, the correlation between Wald statistics in Phase 2
  is derived from the randomization ratio \\r\\ as \\r / (r + 1)\\. If
  `FALSE`, a conservative correlation of 0 is assumed.

- L:

  The termination look in Phase 3.

- zL:

  The z-test statistic at the termination look.

- IMax:

  Maximum information for any active arm versus the common control.

- informationRates:

  The information rates up to look `L`.

- efficacyStopping:

  Indicators of whether efficacy stopping is allowed at each stage up to
  look `L`. Defaults to `TRUE` if left unspecified.

- criticalValues:

  The upper boundaries on the max z-test statistic scale for Phase 2 and
  the z-test statistics for the selected arm in Phase 3 up to look `L`.
  If missing, boundaries will be computed based on the specified alpha
  spending function.

- alpha:

  The significance level. Defaults to 0.025.

- typeAlphaSpending:

  The type of alpha spending for the trial. One of the following: `"OF"`
  for O'Brien-Fleming boundaries, `"P"` for Pocock boundaries, `"WT"`
  for Wang & Tsiatis boundaries, `"sfOF"` for O'Brien-Fleming type
  spending function, `"sfP"` for Pocock type spending function, `"sfKD"`
  for Kim & DeMets spending function, `"sfHSD"` for Hwang, Shi & DeCani
  spending function, and `"none"` for no early efficacy stopping.
  Defaults to `"sfOF"`.

- parameterAlphaSpending:

  The parameter value for the alpha spending. Corresponds to \\\Delta\\
  for `"WT"`, \\\rho\\ for `"sfKD"`, and \\\gamma\\ for `"sfHSD"`.

- spendingTime:

  The error spending time up to look `L`. Defaults to missing, in which
  case, it is the same as `informationRates`.

## Value

A data frame with the following components:

- `pvalue`: p-value for rejecting the null hypothesis.

- `thetahat`: Point estimate of the parameter.

- `cilevel`: Confidence interval level.

- `lower`: Lower bound of confidence interval.

- `upper`: Upper bound of confidence interval.

## Details

If `typeAlphaSpending` is `"OF"`, `"P"`, `"WT"`, or `"none"`, then
`informationRates`, `efficacyStopping`, and `spendingTime` must be of
full length \\K + 1\\, and `informationRates` and `spendingTime` must
end with 1.

## References

Ping Gao, Yingqiu Li. Adaptive two-stage seamless sequential design for
clinical trials. Journal of Biopharmaceutical Statistics, 2025, 35(4),
565-587.

## Author

Kaifeng Lu, <kaifenglu@gmail.com>

## Examples

``` r
getCI_seamless(
  L = 2, zL = 2.075,
  M = 2, r = 1, corr_known = FALSE,
  IMax = 300 / 4, informationRates = c(1/3, 2/3, 1),
  alpha = 0.025, typeAlphaSpending = "sfOF")
#>       pvalue  thetahat cilevel       lower     upper
#> 1 0.03423737 0.2005519    0.95 -0.01534503 0.4651547
```
