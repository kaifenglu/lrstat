# Group Sequential Design for Paired Mean Ratio

Obtains the power given sample size or obtains the sample size given
power for a group sequential design for paired mean ratio.

## Usage

``` r
getDesignPairedMeanRatio(
  beta = NA_real_,
  n = NA_real_,
  pairedRatioH0 = 1,
  pairedRatio = 1.2,
  CV = 1,
  normalApproximation = TRUE,
  rounding = TRUE,
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
  futilityCP = NA_real_,
  futilityPairedRatio = NA_real_,
  typeBetaSpending = "none",
  parameterBetaSpending = NA_real_,
  userBetaSpending = NA_real_,
  spendingTime = NA_real_
)
```

## Arguments

- beta:

  The type II error.

- n:

  The total sample size.

- pairedRatioH0:

  The paired ratio under the null hypothesis.

- pairedRatio:

  The paired ratio under the alternative hypothesis.

- CV:

  The coefficient of variation for paired ratio.

- normalApproximation:

  The type of computation of the p-values. If `TRUE`, the variance is
  assumed to be known, otherwise the calculations are performed with the
  t distribution. The exact calculation using the t distribution is only
  implemented for the fixed design.

- rounding:

  Whether to round up sample size. Defaults to 1 for sample size
  rounding.

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
  futility at stages `1, ..., kMax-1`. Defaults to `rep(-8, kMax-1)` if
  left unspecified. The futility bounds are non-binding for the
  calculation of critical values.

- futilityCP:

  The futility boundary on the conditional power scale.

- futilityPairedRatio:

  The futility boundary on the paired ratio scale.

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

## Value

An S3 class `designPairedMeanRatio` object with three components:

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

  - `numberOfSubjects`: The maximum number of subjects.

  - `expectedNumberOfSubjectsH1`: The expected number of subjects under
    H1.

  - `expectedNumberOfSubjectsH0`: The expected number of subjects under
    H0.

  - `pairedRatioH0`: The paired ratio under the null hypothesis.

  - `pairedRatio`: The paired ratio under the alternative hypothesis.

  - `CV`: The coefficient of variation for paired ratio.

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

  - `numberOfSubjects`: The number of subjects.

  - `efficacyPairedRatio`: The efficacy boundaries on the paired ratio
    scale.

  - `futilityPairedRatio`: The futility boundaries on the paired ratio
    scale.

- `settings`: A list containing the following input parameters:

  - `typeAlphaSpending`: The type of alpha spending.

  - `parameterAlphaSpending`: The parameter value for alpha spending.

  - `userAlphaSpending`: The user defined alpha spending.

  - `typeBetaSpending`: The type of beta spending.

  - `parameterBetaSpending`: The parameter value for beta spending.

  - `userBetaSpending`: The user defined beta spending.

  - `spendingTime`: The error spending time at each analysis.

  - `normalApproximation`: The type of computation of the p-values. If
    `TRUE`, the variance is assumed to be known, otherwise the
    calculations are performed with the t distribution.

  - `rounding`: Whether to round up sample size.

## Author

Kaifeng Lu, <kaifenglu@gmail.com>

## Examples

``` r
# Example 1: group sequential trial power calculation
(design1 <- getDesignPairedMeanRatio(
  beta = 0.1, n = NA, pairedRatio = 1.2, CV = 0.35,
  kMax = 5, alpha = 0.05, typeAlphaSpending = "sfOF"))
#>                                                                                      
#> Group-sequential design with 5 stages for paired mean ratio                          
#> Paired ratio under H0: 1, paired ratio under H1: 1.2, coefficient of variation: 0.35 
#> Overall power: 0.902, overall alpha: 0.05                                            
#> Drift parameter: 2.986, inflation factor: 1.033                                      
#> Maximum information: 268.26, expected under H1: 191.83, expected under H0: 266.02    
#> Maximum # subjects: 31, expected under H1: 22.2, expected under H0: 30.7             
#> Alpha spending: Lan-DeMets O'Brien-Fleming, beta spending: None                      
#>                                                                                      
#>                                  Stage 1 Stage 2 Stage 3 Stage 4 Stage 5
#> Information rate                 0.194   0.387   0.613   0.806   1.000  
#> Efficacy boundary (Z)            4.304   2.942   2.265   1.955   1.742  
#> Cumulative rejection             0.0014  0.1391  0.5325  0.7758  0.9020 
#> Cumulative alpha spent           0.0000  0.0016  0.0123  0.0291  0.0500 
#> Number of subjects               6.0     12.0    19.0    25.0    31.0   
#> Efficacy boundary (paired ratio) 1.817   1.335   1.193   1.142   1.112  
#> Efficacy boundary (p)            0.0000  0.0016  0.0118  0.0253  0.0407 
#> Information                      51.92   103.84  164.42  216.34  268.26 

# Example 2: sample size calculation for one-sample t-test
(design2 <- getDesignPairedMeanRatio(
  beta = 0.1, n = NA, pairedRatio = 1.2, CV = 0.35,
  normalApproximation = FALSE, alpha = 0.05))
#>                                                                                      
#> Fixed design for paired mean ratio                                                   
#> Paired ratio under H0: 1, paired ratio under H1: 1.2, coefficient of variation: 0.35 
#> Overall power: 0.9069, overall alpha: 0.05                                           
#> Drift parameter: 3.034, inflation factor: 1                                          
#> Information: 276.92                                                                  
#> Number of subjects: 32                                                               
#>                                                                                      
#>                                        
#> Efficacy boundary (t)            1.696 
#> Efficacy boundary (paired ratio) 1.107 
#> Efficacy boundary (p)            0.0500
```
