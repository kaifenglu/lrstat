# Conditional Power for a Phase 2/3 Seamless Design

Obtains the conditional power for specified incremental information
given the interim results, parameter values, and data-dependent changes
in the error spending function, as well as the number and spacing of
interim looks.

## Usage

``` r
getCP_seamless(
  INew = NA_real_,
  M = NA_integer_,
  r = 1,
  corr_known = TRUE,
  L = NA_integer_,
  zL = NA_real_,
  theta = NA_real_,
  IMax = NA_real_,
  K = NA_integer_,
  informationRates = NA_real_,
  efficacyStopping = NA_integer_,
  futilityStopping = NA_integer_,
  criticalValues = NULL,
  alpha = 0.025,
  typeAlphaSpending = "sfOF",
  parameterAlphaSpending = NA_real_,
  userAlphaSpending = NA_real_,
  futilityBounds = NULL,
  futilityCP = NULL,
  futilityTheta = NULL,
  spendingTime = NA_real_,
  MullerSchafer = FALSE,
  kNew = NA_integer_,
  informationRatesNew = NA_real_,
  efficacyStoppingNew = NA_integer_,
  futilityStoppingNew = NA_integer_,
  typeAlphaSpendingNew = "sfOF",
  parameterAlphaSpendingNew = NA_real_,
  futilityBoundsInt = NULL,
  futilityCPInt = NULL,
  futilityThetaInt = NULL,
  typeBetaSpendingNew = "none",
  parameterBetaSpendingNew = NA_real_,
  spendingTimeNew = NA_real_
)
```

## Arguments

- INew:

  The maximum information for the active arm versus the common control
  in the secondary trial.

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

- theta:

  The assumed treatment effect for the selected arm versus the common
  control.

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

- futilityStopping:

  Indicators of whether futility stopping is allowed at each stage of
  the primary trial. Defaults to true if left unspecified.

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

- futilityBounds:

  The lower boundaries on the max z-test statistic scale for Phase 2 and
  the z-test statistics for the selected arm in Phase 3 for the primary
  trial.

- futilityCP:

  The conditional power-based futility bounds for the primary trial.

- futilityTheta:

  The parameter value-based futility bounds for the primary trial.

- spendingTime:

  The error spending time of the primary trial. Defaults to missing, in
  which case it is assumed to be the same as `informationRates`.

- MullerSchafer:

  Whether to use the Muller and Schafer (2001) method for trial
  adaptation.

- kNew:

  The number of looks of the secondary trial.

- informationRatesNew:

  The spacing of looks of the secondary trial.

- efficacyStoppingNew:

  The indicators of whether efficacy stopping is allowed at each look of
  the secondary trial. Defaults to `TRUE` if left unspecified.

- futilityStoppingNew:

  The indicators of whether futility stopping is allowed at each look of
  the secondary trial. Defaults to true if left unspecified.

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

- futilityBoundsInt:

  The futility boundaries on the z statistic scale for new stages of the
  integrated trial.

- futilityCPInt:

  The conditional power-based futility bounds for new stages of the
  integrated trial.

- futilityThetaInt:

  The parameter value-based futility bounds for the new stages of the
  integrated trial.

- typeBetaSpendingNew:

  The type of beta spending for the secondary trial. One of the
  following: `"sfOF"` for O'Brien-Fleming type spending function,
  `"sfP"` for Pocock type spending function, `"sfKD"` for Kim & DeMets
  spending function, `"sfHSD"` for Hwang, Shi & DeCani spending
  function, and `"none"` for no early futility stopping. Defaults to
  `"none"`.

- parameterBetaSpendingNew:

  The parameter value of beta spending for the secondary trial.
  Corresponds to \\\rho\\ for `"sfKD"`, and \\\gamma\\ for `"sfHSD"`.

- spendingTimeNew:

  The error spending time of the secondary trial. Defaults to missing,
  in which case it is assumed to be the same as `informationRatesNew`.

## Value

A vector of two conditional powers given the interim results and
parameter values, one without design change and the other with
data-dependent design changes.

## References

Ping Gao, Yingqiu Li. Adaptive two-stage seamless sequential design for
clinical trials. Journal of Biopharmaceutical Statistics, 2025, 35(4),
565-587.

## See also

[`adaptDesign_seamless`](https://kaifenglu.github.io/lrstat/reference/adaptDesign_seamless.md)

## Author

Kaifeng Lu, <kaifenglu@gmail.com>

## Examples

``` r
getCP_seamless(
  INew = 198 / 4, M = 2, r = 1, corr_known = FALSE,
  L = 1, zL = -log(0.67) * sqrt(80 / 4), theta = -log(0.691),
  IMax = 120 / 4, K = 2, informationRates = c(1/3, 2/3, 1),
  alpha = 0.025, typeAlphaSpending = "OF", kNew = 1)
#> Error: b must be provided
```
