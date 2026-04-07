# Power and Sample Size for a Generic Group Sequential Design

Obtains the maximum information and stopping boundaries for a generic
group sequential design assuming a constant treatment effect, or obtains
the power given the maximum information and stopping boundaries.

## Usage

``` r
getDesign(
  beta = NA_real_,
  IMax = NA_real_,
  theta = NA_real_,
  kMax = 1L,
  informationRates = NA_real_,
  efficacyStopping = NA_integer_,
  futilityStopping = NA_integer_,
  criticalValues = NA_real_,
  alpha = 0.025,
  typeAlphaSpending = "sfOF",
  parameterAlphaSpending = NA_real_,
  userAlphaSpending = NA_real_,
  futilityBounds = NA_real_,
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

- futilityBounds:

  Lower boundaries on the z-test statistic scale for stopping for
  futility at stages `1, ..., kMax-1`. Defaults to `rep(-6, kMax-1)` if
  left unspecified. The futility bounds are non-binding for the
  calculation of critical values.

- typeBetaSpending:

  The type of beta spending. One of the following: `"sfOF"` for
  O'Brien-Fleming type spending function, `"sfP"` for Pocock type
  spending function, `"sfKD"` for Kim & DeMets spending function,
  `"sfHSD"` for Hwang, Shi & DeCani spending function, `"user"` for user
  defined spending, and `"none"` for no early futility stopping.
  Defaults to `"none"`.

- parameterBetaSpending:

  The parameter value for the beta spending. Corresponds to \\\rho\\ for
  `"sfKD"`, and \\\gamma\\ for `"sfHSD"`.

- userBetaSpending:

  The user defined beta spending. Cumulative beta spent up to each
  stage.

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

## References

Christopher Jennison, Bruce W. Turnbull. Group Sequential Methods with
Applications to Clinical Trials. Chapman & Hall/CRC: Boca Raton, 2000,
ISBN:0849303168

## Author

Kaifeng Lu, <kaifenglu@gmail.com>

## Examples

``` r
# Example 1: obtain the maximum information given power
(design1 <- getDesign(
  beta = 0.2, theta = -log(0.7),
  kMax = 2, informationRates = c(0.5,1),
  alpha = 0.025, typeAlphaSpending = "sfOF",
  typeBetaSpending = "sfP"))
#>                                                                              
#> Group-sequential design with 2 stages                                        
#> theta: 0.357, maximum information: 71.97                                     
#> Overall power: 0.8, overall alpha (1-sided): 0.025, attained alpha: 0.0205   
#> Drift parameter: 3.026, inflation factor: 1.167                              
#> Expected information under H1: 60.12, expected information under H0: 41.78   
#> Alpha spending: Lan-DeMets O'Brien-Fleming, beta spending: Lan-DeMets Pocock 
#>                                                                              
#>                               Stage 1 Stage 2
#> Information rate              0.500   1.000  
#> Efficacy boundary (Z)         2.963   1.969  
#> Futility boundary (Z)         0.985   1.969  
#> Cumulative rejection          0.2053  0.8000 
#> Cumulative futility           0.1240  0.2000 
#> Cumulative alpha spent        0.0015  0.0250 
#> Efficacy boundary (theta)     0.494   0.232  
#> Futility boundary (theta)     0.164   0.232  
#> Efficacy boundary (p)         0.0015  0.0245 
#> Futility boundary (p)         0.1624  0.0245 
#> Information                   35.99   71.97  
#> Cumulative rejection under H0 0.0015  0.0205 
#> Cumulative futility under H0  0.8376  0.9795 

# Example 2: obtain power given the maximum information
(design2 <- getDesign(
  IMax = 72.5, theta = -log(0.7),
  kMax = 3, informationRates = c(0.5, 0.75, 1),
  alpha = 0.025, typeAlphaSpending = "sfOF",
  typeBetaSpending = "sfP"))
#>                                                                              
#> Group-sequential design with 3 stages                                        
#> theta: 0.357, maximum information: 72.5                                      
#> Overall power: 0.771, overall alpha (1-sided): 0.025, attained alpha: 0.0187 
#> Drift parameter: 3.037, inflation factor: 1.263                              
#> Expected information under H1: 51.82, expected information under H0: 39.54   
#> Alpha spending: Lan-DeMets O'Brien-Fleming, beta spending: Lan-DeMets Pocock 
#>                                                                              
#>                               Stage 1 Stage 2 Stage 3
#> Information rate              0.500   0.750   1.000  
#> Efficacy boundary (Z)         2.963   2.359   2.014  
#> Futility boundary (Z)         1.076   1.511   2.014  
#> Cumulative rejection          0.2075  0.6018  0.7710 
#> Cumulative futility           0.1420  0.1896  0.2290 
#> Cumulative alpha spent        0.0015  0.0096  0.0250 
#> Efficacy boundary (theta)     0.492   0.320   0.237  
#> Futility boundary (theta)     0.179   0.205   0.237  
#> Efficacy boundary (p)         0.0015  0.0092  0.0220 
#> Futility boundary (p)         0.1410  0.0654  0.0220 
#> Information                   36.25   54.38   72.50  
#> Cumulative rejection under H0 0.0015  0.0093  0.0187 
#> Cumulative futility under H0  0.8590  0.9488  0.9813 
```
