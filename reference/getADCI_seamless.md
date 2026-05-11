# Confidence Interval After Adaptation for a Phase 2/3 Seamless Design

Obtains the p-value, conservative point estimate, and confidence
interval after the end of an adaptive phase 2/3 seamless design.

## Usage

``` r
getADCI_seamless(
  M = NA_integer_,
  r = 1,
  corr_known = TRUE,
  L = NA_integer_,
  zL = NA_real_,
  IMax = NA_real_,
  K = NA_integer_,
  informationRates = NA_real_,
  efficacyStopping = NA_integer_,
  criticalValues = NA_real_,
  alpha = 0.25,
  typeAlphaSpending = "sfOF",
  parameterAlphaSpending = NA_real_,
  spendingTime = NA_real_,
  MullerSchafer = FALSE,
  Lc = NA_integer_,
  zLc = NA_real_,
  INew = NA_real_,
  informationRatesNew = NA_real_,
  efficacyStoppingNew = NA_integer_,
  typeAlphaSpendingNew = "sfOF",
  parameterAlphaSpendingNew = NA_real_,
  spendingTimeNew = NA_real_
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
  `FALSE`, a conservative correlation of 0 is used.

- L:

  The interim adaptation look in Phase 3.

- zL:

  The z-test statistic at the interim adaptation look of Phase 3.

- IMax:

  Maximum information for the active arm versus the common control for
  the original trial. Must be provided.

- K:

  Number of sequential looks in Phase 3.

- informationRates:

  A numeric vector of information rates fixed before the trial. If
  unspecified, defaults to \\(1:(K+1)) / (K+1)\\.

- efficacyStopping:

  Indicators of whether efficacy stopping is allowed at each stage of
  the primary trial. Defaults to `TRUE` if left unspecified.

- criticalValues:

  The upper boundaries on the max z-test statistic scale for Phase 2 and
  the z-test statistics for the selected arm in Phase 3 for the primary
  trial. If missing, boundaries will be computed based on the specified
  alpha spending function.

- alpha:

  The significance level of the primary trial. Defaults to 0.025.

- typeAlphaSpending:

  The type of alpha spending for the primary trial. One of the
  following: `"OF"` for O'Brien-Fleming boundaries, `"P"` for Pocock
  boundaries, `"WT"` for Wang & Tsiatis boundaries, `"sfOF"` for
  O'Brien-Fleming type spending function, `"sfP"` for Pocock type
  spending function, `"sfKD"` for Kim & DeMets spending function,
  `"sfHSD"` for Hwang, Shi & DeCani spending function, and `"none"` for
  no early efficacy stopping. Defaults to `"sfOF"`.

- parameterAlphaSpending:

  The parameter value of alpha spending for the primary trial.
  Corresponds to \\\Delta\\ for `"WT"`, \\\rho\\ for `"sfKD"`, and
  \\\gamma\\ for `"sfHSD"`.

- spendingTime:

  The error spending time of the primary trial. Defaults to missing, in
  which case, it is the same as `informationRates`.

- MullerSchafer:

  Whether to use the Muller and Schafer (2001) method for trial
  adaptation.

- Lc:

  The termination look of the integrated trial.

- zLc:

  The z-test statistic at the termination look of the integrated trial.

- INew:

  The maximum information for the active arm versus the common control
  in the secondary trial.

- informationRatesNew:

  The spacing of looks of the secondary trial.

- efficacyStoppingNew:

  The indicators of whether efficacy stopping is allowed at each look of
  the secondary trial. Defaults to `TRUE` if left unspecified.

- typeAlphaSpendingNew:

  The type of alpha spending for the secondary trial. One of the
  following: `"OF"` for O'Brien-Fleming boundaries, `"P"` for Pocock
  boundaries, `"WT"` for Wang & Tsiatis boundaries, `"sfOF"` for
  O'Brien-Fleming type spending function, `"sfP"` for Pocock type
  spending function, `"sfKD"` for Kim & DeMets spending function,
  `"sfHSD"` for Hwang, Shi & DeCani spending function, and `"none"` for
  no early efficacy stopping. Defaults to `"sfOF"`.

- parameterAlphaSpendingNew:

  The parameter value of alpha spending for the secondary trial.
  Corresponds to \\\Delta\\ for `"WT"`, \\\rho\\ for `"sfKD"`, and
  \\\gamma\\ for `"sfHSD"`.

- spendingTimeNew:

  The error spending time of the secondary trial. Defaults to missing,
  in which case, it is the same as `informationRatesNew`.

## Value

A data frame with the following variables:

- `pvalue`: p-value for rejecting the null hypothesis.

- `thetahat`: Point estimate of the parameter.

- `cilevel`: Confidence interval level.

- `lower`: Lower bound of confidence interval.

- `upper`: Upper bound of confidence interval.

## Details

If typeAlphaSpendingNew is `"OF"`, `"P"`, `"WT"`, or `"none"`, then
`informationRatesNew`, `efficacyStoppingNew`, and `spendingTimeNew` must
be of full length `kNew`, and `informationRatesNew` and
`spendingTimeNew` must end with 1.

## References

Ping Gao, Yingqiu Li. Adaptive multiple comparison sequential design
(AMCSD) for clinical trials. Journal of Biopharmaceutical Statistics,
2024, 34(3), 424-440.

## Author

Kaifeng Lu, <kaifenglu@gmail.com>

## Examples

``` r
getADCI_seamless(
  M = 2, r = 1, corr_known = FALSE,
  L = 1, zL = -log(0.67) * sqrt(80 / 4),
  IMax = 120 / 4, K = 2, informationRates = c(1/3, 2/3, 1),
  alpha = 0.025, typeAlphaSpending = "OF",
  Lc = 2, zLc = -log(0.677) * sqrt(236 / 4), INew = 236 / 4)
#> Error: Element with name 'exitProb' not found.
```
