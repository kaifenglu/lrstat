# Group Sequential Design for One-Sample Proportion

Obtains the power given sample size or obtains the sample size given
power for a group sequential design for one-sample proportion.

## Usage

``` r
getDesignOneProportion(
  beta = NA_real_,
  n = NA_real_,
  piH0 = 0.1,
  pi = 0.2,
  nullVariance = TRUE,
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

- piH0:

  The response probability under the null hypothesis.

- pi:

  The response probability under the alternative hypothesis.

- nullVariance:

  Whether to use the variance under the null or the variance under the
  alternative.

- normalApproximation:

  The type of computation of the p-values. If `TRUE`, the normal
  approximation will be used, otherwise the calculations are performed
  with the binomial distribution. The exact calculation using the
  binomial distribution is only implemented for the fixed design.

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

## Value

An S3 class `designOneProportion` object with three components:

- `overallResults`: A data frame containing the following variables:

  - `overallReject`: The overall rejection probability.

  - `alpha`: The overall significance level.

  - `attainedAlpha`: The attained significance level, which is different
    from the overall significance level in the presence of futility
    stopping as well as for the binomial exact test in a fixed design.

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

  - `piH0`: The response probability under the null hypothesis.

  - `pi`: The response probability under the alternative hypothesis.

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

  - `efficacyResponses`: The efficacy boundaries on the number of
    responses scale.

  - `futilityResponses`: The futility boundaries on the number of
    responses scale.

  - `numberOfSubjects`: The number of subjects.

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

  - `nullVariance`: Whether to use the variance under the null or the
    empirical variance under the alternative.

  - `normalApproximation`: The type of computation of the p-values. If
    `TRUE`, the variance is assumed to be known, otherwise the
    calculations are performed with the binomial distribution.

  - `rounding`: Whether to round up sample size.

## Author

Kaifeng Lu, <kaifenglu@gmail.com>

## Examples

``` r
# Example 1: group sequential trial power calculation
(design1 <- getDesignOneProportion(
  beta = 0.2, n = NA, piH0 = 0.15, pi = 0.25,
  kMax = 3, alpha = 0.05, typeAlphaSpending = "sfOF"))
#>                                                                                   
#> Group-sequential design with 3 stages for one-sample proportion                   
#> Response probability under H0: 0.15, response probability under H1: 0.25          
#> Overall power: 0.8012, overall alpha (1-sided): 0.05                              
#> Drift parameter: 2.203, inflation factor: 1.001                                   
#> Maximum information: 485.33, expected under H1: 388.39, expected under H0: 482.56 
#> Maximum # subjects: 91, expected under H1: 72.8, expected under H0: 90.5          
#> Test statistic: z-test                                                            
#> Alpha spending: Lan-DeMets O'Brien-Fleming, beta spending: None                   
#>                                                                                   
#>                                 Stage 1 Stage 2 Stage 3
#> Information rate                0.330   0.670   1.000  
#> Efficacy boundary (Z)           3.220   2.133   1.696  
#> Cumulative rejection            0.0822  0.5209  0.8012 
#> Cumulative alpha spent          0.0006  0.0167  0.0500 
#> Number of subjects              30.0    61.0    91.0   
#> Efficacy boundary (# responses) 11      16      20     
#> Efficacy boundary (p)           0.0006  0.0165  0.0450 
#> Information                     160.00  325.33  485.33 

# Example 2: sample size calculation for one-sample binomial exact test
(design2 <- getDesignOneProportion(
  beta = 0.2, n = NA, piH0 = 0.15, pi = 0.25,
  normalApproximation = FALSE, alpha = 0.05))
#>                                                                              
#> Fixed design for one-sample proportion                                       
#> Response probability under H0: 0.15, response probability under H1: 0.25     
#> Overall power: 0.8097, overall alpha (1-sided): 0.05, attained alpha: 0.0355 
#> Drift parameter: 2.422, inflation factor: 1                                  
#> Information: 586.67                                                          
#> Number of subjects: 110                                                      
#> Test statistic: exact test                                                   
#>                                                                              
#>                                       
#> Efficacy boundary (Z)           1.651 
#> Efficacy boundary (# responses) 24    
#> Efficacy boundary (p)           0.0355
```
