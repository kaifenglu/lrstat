# Confidence Interval After Adaptation for Multi-Arm Multi-Stage Design

Obtains the p-value, conservative point estimate, and confidence
interval after the end of an adaptive multi-arm multi-stage trial.

## Usage

``` r
getADCI_mams(
  M = NA_integer_,
  r = NA_real_,
  corr_known = TRUE,
  L = NA_integer_,
  zL = NA_real_,
  IMax = NA_real_,
  kMax = NA_integer_,
  informationRates = NA_real_,
  efficacyStopping = NA_integer_,
  criticalValues = NULL,
  alpha = 0.25,
  typeAlphaSpending = "sfOF",
  parameterAlphaSpending = NA_real_,
  spendingTime = NA_real_,
  MullerSchafer = FALSE,
  MNew = NA_integer_,
  selected = NA_integer_,
  rNew = NA_real_,
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

  The matrix of by-level upper boundaries on the z-test statistic scale
  for efficacy stopping up to look `L` for the primary trial. The first
  column is for level `M`, the second column is for level `M - 1`, and
  so on, with the last column for level 1. If left unspecified, the
  critical values will be computed based on the specified alpha spending
  function.

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
  Corresponds to \\\Delta\\ for "WT", \\\rho\\ for "sfKD", and
  \\\gamma\\ for "sfHSD".

- spendingTime:

  The error spending time of the primary trial. Defaults to missing, in
  which case, it is the same as `informationRates`.

- MullerSchafer:

  Whether to use the Muller and Schafer (2001) method for trial
  adaptation.

- MNew:

  The number of active treatment arms in the secondary trial.

- selected:

  The indices of the selected treatment arms for the secondary trial
  among the `M` active arms in the primary trial.

- rNew:

  The randomization ratio of each active arm to the common control in
  the secondary trial.

- Lc:

  The termination look of the integrated trial.

- zLc:

  The z-test statistics at the termination look of the integrated trial.

- INew:

  The maximum information for any active arm versus the common control
  in the secondary trial.

- informationRatesNew:

  The spacing of looks of the secondary trial up to look `L2`.

- efficacyStoppingNew:

  The indicators of whether efficacy stopping is allowed at each look of
  the secondary trial up to look `L2`. Defaults to `TRUE` if left
  unspecified.

- typeAlphaSpendingNew:

  The type of alpha spending for the secondary trial. One of the
  following: `"OF"` for O'Brien-Fleming boundaries, `"sfOF"` for
  O'Brien-Fleming type spending function, `"sfP"` for Pocock type
  spending function, `"sfKD"` for Kim & DeMets spending function,
  `"sfHSD"` for Hwang, Shi & DeCani spending function, and `"none"` for
  no early efficacy stopping. Defaults to `"sfOF"`.

- parameterAlphaSpendingNew:

  The parameter value of alpha spending for the secondary trial.
  Corresponds to \\\Delta\\ for "WT", \\\rho\\ for "sfKD", and
  \\\gamma\\ for "sfHSD".

- spendingTimeNew:

  The error spending time of the secondary trial up to look `L2`.
  Defaults to missing, in which case, it is the same as
  `informationRatesNew`.

## Value

A data frame with the following variables:

- `level`: Number of individual hypotheses considered for multiplicity.

- `index`: The index of the treatment arm among the `M` active arms.

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
getADCI_mams(
  M = 2, r = 1, corr_known = FALSE, L = 1, zL = c(2.075, 2.264),
  IMax = 300 / 4, kMax = 2, informationRates = c(0.5, 1),
  alpha = 0.025, typeAlphaSpending = "sfOF",
  MNew = 1, selected = 2, rNew = 1,
  Lc = 2, zLc = 1.667, INew = 374 / 4)
#>   level index     pvalue  thetahat cilevel         lower     upper
#> 1     1     2 0.02551277 0.1712693    0.95 -0.0007892748 0.3530049
```
