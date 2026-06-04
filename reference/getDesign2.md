# Power and Sample Size for Group Sequential Design With Futility Stopping Under Null Hypothesis

Obtains the maximum information and stopping boundaries for a generic
group sequential design with futility stopping under the null hypothesis
assuming a constant treatment effect, or obtains the power given the
maximum information and stopping boundaries.

## Usage

``` r
getDesign2(
  beta = NA_real_,
  IMax = NA_real_,
  theta = NA_real_,
  kMax = 1L,
  informationRates = NA_real_,
  efficacyStopping = NA_integer_,
  futilityStopping = NA_integer_,
  criticalValues = NULL,
  alpha = 0.025,
  typeAlphaSpending = "sfOF",
  parameterAlphaSpending = NA_real_,
  userAlphaSpending = NA_real_,
  symmetricBounds = TRUE,
  astar = 0.025,
  futilityBounds = NULL,
  typeBetaSpending = "none",
  parameterBetaSpending = NA_real_,
  userBetaSpending = NA_real_,
  spendingTime = NA_real_,
  varianceRatio = 1
)
```

## Arguments

- beta:

  The type II error.

- IMax:

  The maximum information. Either `beta` or `IMax` should be provided
  while the other one should be missing.

- theta:

  The parameter value. Null hypothesis is at `theta = 0`, and the
  alternative hypothesis is one-sided for `theta > 0`.

- kMax:

  The maximum number of stages.

- informationRates:

  The information rates. Fixed prior to the trial. Defaults to
  `(1:kMax) / kMax` if left unspecified.

- efficacyStopping:

  Indicators of whether efficacy stopping is allowed at each stage.
  Defaults to `TRUE` if left unspecified.

- futilityStopping:

  Indicators of whether futility stopping is allowed at each stage.
  Defaults to `TRUE` if left unspecified.

- criticalValues:

  Upper boundaries on the z-test statistic scale for stopping for
  efficacy.

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

- userAlphaSpending:

  The user defined alpha spending. Cumulative alpha spent up to each
  stage.

- symmetricBounds:

  If `TRUE`, futility bounds are set to the negative of efficacy bounds
  at each analysis (subject to `futilityStopping`). If `FALSE`, futility
  bounds are determined by `futilityBounds` or beta spending.

- astar:

  The overall futility stopping probability under the null hypothesis.

- futilityBounds:

  A vector of length `kMax` for the futility stopping boundaries on the
  Z-scale under the null hypothesis. Defaults to `rep(-8, kMax)` if left
  unspecified. The futility bounds are non-binding for the calculation
  of critical values.

- typeBetaSpending:

  The type of beta spending function for determining futility bounds
  under the null hypothesis when `futilityBounds` is not provided. The
  same types as `typeAlphaSpending` are allowed, except that `"none"`
  corresponds to no futility stopping under the null hypothesis.

- parameterBetaSpending:

  The parameter for the beta spending function. Corresponds to
  \\\Delta\\ for `"WT"`, \\\rho\\ for `"sfKD"`, and \\\gamma\\ for
  `"sfHSD"`.

- userBetaSpending:

  A vector of length `kMax` for the user defined beta spending function
  when `typeBetaSpending == "user"`. The last element must be equal to
  `astar` and the vector must be increasing.

- spendingTime:

  A vector of length `kMax` for the error spending time at each
  analysis. Defaults to missing, in which case, it is the same as
  `informationRates`.

- varianceRatio:

  The ratio of the variance under H0 to the variance under H1.

## Value

An S3 class `design` object with three components:

- `overallResults`: A data frame containing the following variables:

  - `overallReject`: The overall rejection probability.

  - `alpha`: The overall significance level.

  - `attainedAlpha`: The attained significance level, which is different
    from the overall significance level in the presence of futility
    stopping.

  - `astar`: The overall futility stopping probability under the null
    hypothesis.

  - `kMax`: The number of stages.

  - `theta`: The parameter value.

  - `information`: The maximum information.

  - `expectedInformationH1`: The expected information under H1.

  - `expectedInformationH0`: The expected information under H0.

  - `drift`: The drift parameter, equal to `theta*sqrt(information)`.

  - `inflationFactor`: The inflation factor (relative to the fixed
    design).

- `byStageResults`: A data frame containing the following variables:

  - `informationRates`: The information rates.

  - `efficacyBounds`: The efficacy boundaries on the Z-scale.

  - `futilityBounds`: The futility boundaries on the Z-scale.

  - `rejectPerStage`: The probability for efficacy stopping.

  - `futilityPerStage`: The probability for futility stopping.

  - `cumulativeRejection`: The cumulative probability for efficacy
    stopping.

  - `cumulativeFutility`: The cumulative probability for futility
    stopping.

  - `cumulativeAlphaSpent`: The cumulative alpha spent.

  - `efficacyTheta`: The efficacy boundaries on the parameter scale.

  - `futilityTheta`: The futility boundaries on the parameter scale.

  - `efficacyP`: The efficacy boundaries on the p-value scale.

  - `futilityP`: The futility boundaries on the p-value scale.

  - `information`: The cumulative information.

  - `efficacyStopping`: Whether to allow efficacy stopping.

  - `futilityStopping`: Whether to allow futility stopping.

  - `rejectPerStageH0`: The probability for efficacy stopping under H0.

  - `futilityPerStageH0`: The probability for futility stopping under
    H0.

  - `cumulativeRejectionH0`: The cumulative probability for efficacy
    stopping under H0.

  - `cumulativeFutilityH0`: The cumulative probability for futility
    stopping under H0.

- `settings`: A list containing the following input parameters:

  - `typeAlphaSpending`: The type of alpha spending.

  - `parameterAlphaSpending`: The parameter value for alpha spending.

  - `userAlphaSpending`: The user defined alpha spending.

  - `typeBetaSpending`: The type of beta spending.

  - `parameterBetaSpending`: The parameter value for beta spending.

  - `userBetaSpending`: The user defined beta spending.

  - `spendingTime`: The error spending time at each analysis.

  - `varianceRatio`: The ratio of the variance under H0 to the variance
    under H1.

## Details

The futility stopping boundaries under the null hypothesis are
non-binding. The function determines efficacy and futility bounds based
on the inputs provided, following a clear priority order.

**Efficacy bounds:** If `criticalValues` are supplied, they take
precedence and all alpha-spending parameters are ignored. Otherwise,
efficacy bounds are derived from the specified alpha-spending function.

**Futility bounds:** Futility inputs are evaluated in the following
order of priority:

1.  If `futilityBounds` are provided, they override beta-spending
    parameters.

2.  If `futilityBounds` are not specified, futility bounds are computed
    using the beta-spending approach with `astar` being the maximum
    futility stopping probability under the null hypothesis. If
    `typeBetaSpending == "none"`, then there is no futility stopping
    under the null hypothesis.

If `symmetricBounds = TRUE`, the futility bounds are set to
`-efficacyBounds` and beta-spending inputs are ignored.

## References

Christopher Jennison, Bruce W. Turnbull. Group Sequential Methods with
Applications to Clinical Trials. Chapman & Hall/CRC: Boca Raton, 2000,
ISBN:0849303168

## Author

Kaifeng Lu, <kaifenglu@gmail.com>

## Examples

``` r
(design1 <- getDesign2(
  beta = 0.149, theta = -log(0.65),
  kMax = 2, informationRates = c(0.87, 1),
  alpha = 0.004, typeAlphaSpending = "sfOF",
  astar = 0.1, typeBetaSpending = "sfHSD",
  parameterBetaSpending = -8))
#>                                                                             
#> Group-sequential design with 2 stages                                       
#> theta: 0.431, maximum information: 74.5                                     
#> Overall power: 0.851, overall alpha (1-sided): 0.004, attained alpha: 0.004 
#> Drift parameter: 3.718, inflation factor: 1.014                             
#> Expected information under H1: 56.39, expected information under H0: 0.56   
#> Alpha spending: Lan-DeMets O'Brien-Fleming, beta spending: HSD(gamma = -8)  
#>                                                                             
#>                               Stage 1 Stage 2
#> Information rate              0.870   1.000  
#> Efficacy boundary (Z)         2.873   2.706  
#> Futility boundary (Z)         -2.873  -2.706 
#> Cumulative rejection          0.7240  0.8510 
#> Cumulative futility           0.0000  0.0000 
#> Cumulative alpha spent        0.0020  0.0040 
#> Efficacy boundary (theta)     0.357   0.314  
#> Futility boundary (theta)     -0.357  -0.314 
#> Efficacy boundary (p)         0.0020  0.0034 
#> Futility boundary (p)         0.9980  0.9966 
#> Information                   64.81   74.50  
#> Cumulative rejection under H0 0.0020  0.0040 
#> Cumulative futility under H0  0.0020  0.0040 
```
