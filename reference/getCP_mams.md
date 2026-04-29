# Conditional Power for a Multi-Arm Multi-Stage Design

Obtains the conditional power for specified incremental information
given the interim results, parameter values, and data-dependent changes
in the selected treatment(s), the error spending function, as well as
the number and spacing of interim looks.

## Usage

``` r
getCP_mams(
  INew = NA_real_,
  M = NA_integer_,
  r = 1,
  corr_known = TRUE,
  L = NA_integer_,
  zL = NA_real_,
  theta = NA_real_,
  IMax = NA_real_,
  kMax = NA_integer_,
  informationRates = NA_real_,
  efficacyStopping = NA_integer_,
  criticalValues = NA_real_,
  alpha = 0.025,
  typeAlphaSpending = "sfOF",
  parameterAlphaSpending = NA_real_,
  userAlphaSpending = NA_real_,
  spendingTime = NA_real_,
  MullerSchafer = FALSE,
  MNew = NA_integer_,
  selected = NA_integer_,
  rNew = 1,
  kNew = NA_integer_,
  informationRatesNew = NA_real_,
  efficacyStoppingNew = NA_integer_,
  typeAlphaSpendingNew = "sfOF",
  parameterAlphaSpendingNew = NA_real_,
  spendingTimeNew = NA_real_
)
```

## Arguments

- INew:

  The maximum information for any active arm versus the common control
  in the secondary trial.

- M:

  Number of active treatment arms in the primary trial.

- r:

  Randomization ratio of each active arm to the common control in the
  primary trial.

- corr_known:

  Logical. If `TRUE`, the correlation between Wald statistics is derived
  from the randomization ratio \\r\\ as \\r / (r + 1)\\. If `FALSE`, a
  conservative correlation of 0 is assumed.

- L:

  The interim adaptation look of the primary trial.

- zL:

  The z-test statistics at the interim adaptation look of the primary
  trial.

- theta:

  A vector of length \\M\\ representing the assumed treatment effects
  for each active arm versus the common control. The global null is
  \\\theta_i = 0\\ for all \\i\\, and alternatives are one-sided:
  \\\theta_i \> 0\\ for at least one \\i = 1, \ldots, M\\.

- IMax:

  Maximum information for any active arm versus the common control for
  the primary trial. Must be provided.

- kMax:

  The maximum number of stages of the primary trial.

- informationRates:

  The information rates of the primary trial.

- efficacyStopping:

  Indicators of whether efficacy stopping is allowed at each stage of
  the primary trial. Defaults to `TRUE` if left unspecified.

- criticalValues:

  The upper boundaries on the max z-test statistic scale for efficacy
  stopping for the primary trial. If missing, boundaries will be
  computed based on the specified alpha spending function.

- alpha:

  The significance level of the primary trial. Defaults to 0.025.

- typeAlphaSpending:

  The type of alpha spending for the primary trial. One of the
  following: `"OF"` for O'Brien-Fleming boundaries, `"P"` for Pocock
  boundaries, `"WT"` for Wang & Tsiatis boundaries, `"sfOF"` for
  O'Brien-Fleming type spending function, `"sfP"` for Pocock type
  spending function, `"sfKD"` for Kim & DeMets spending function,
  `"sfHSD"` for Hwang, Shi & DeCani spending function, `"user"` for user
  defined spending, and `"none"` for no early efficacy stopping.
  Defaults to `"sfOF"`.

- parameterAlphaSpending:

  The parameter value of alpha spending for the primary trial.
  Corresponds to \\\Delta\\ for `"WT"`, \\\rho\\ for `"sfKD"`, and
  \\\gamma\\ for `"sfHSD"`.

- userAlphaSpending:

  The user-defined alpha spending for the primary trial. Represents the
  cumulative alpha spent up to each stage.

- spendingTime:

  The error spending time of the primary trial. Defaults to missing, in
  which case it is assumed to be the same as `informationRates`.

- MullerSchafer:

  Whether to use the Muller and Schafer (2001) method for trial
  adaptation.

- MNew:

  Number of active treatment arms in the secondary trial.

- selected:

  The indices of the selected active treatment arms for the secondary
  trial.

- rNew:

  Randomization ratio of each active arm to the common control in the
  secondary trial.

- kNew:

  The number of looks of the secondary trial.

- informationRatesNew:

  The spacing of looks of the secondary trial.

- efficacyStoppingNew:

  The indicators of whether efficacy stopping is allowed at each look of
  the secondary trial. Defaults to `TRUE` if left unspecified.

- typeAlphaSpendingNew:

  The type of alpha spending for the secondary trial. One of the
  following: `"OF"` for O'Brien-Fleming boundaries, `"sfOF"` for
  O'Brien-Fleming type spending function, `"sfP"` for Pocock type
  spending function, `"sfKD"` for Kim & DeMets spending function,
  `"sfHSD"` for Hwang, Shi & DeCani spending function, and `"none"` for
  no early efficacy stopping. Defaults to `"sfOF"`.

- parameterAlphaSpendingNew:

  The parameter value of alpha spending for the secondary trial.
  Corresponds to \\\rho\\ for `"sfKD"`, and \\\gamma\\ for `"sfHSD"`.

- spendingTimeNew:

  The error spending time of the secondary trial. Defaults to missing,
  in which case it is assumed to be the same as `informationRatesNew`.

## Value

The conditional power given the interim results, parameter values, and
data-dependent design changes.

## References

Ping Gao, Yingqiu Li. Adaptive multiple comparison sequential design
(AMCSD) for clinical trials. Journal of Biopharmaceutical Statistics,
2024, 34(3), 424-440.

## See also

[`adaptDesign_mams`](https://kaifenglu.github.io/lrstat/reference/adaptDesign_mams.md)

## Author

Kaifeng Lu, <kaifenglu@gmail.com>

## Examples

``` r
getCP_mams(
  INew = 373 / 4, M = 2, r = 1, corr_known = FALSE,
  L = 1, zL = c(-log(0.91), -log(0.78)) * sqrt(324 / 4 / 2),
  theta = c(-log(0.91), -log(0.78)),
  IMax = 324 / 4, kMax = 2, informationRates = c(1/2, 1),
  alpha = 0.025, typeAlphaSpending = "OF",
  MNew = 1, selected = 2, rNew = 1)
#> [1] 0.8000917
```
