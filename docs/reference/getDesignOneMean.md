# Group Sequential Design for One-Sample Mean

Obtains the power given sample size or obtains the sample size given
power for a group sequential design for one-sample mean.

## Usage

``` r
getDesignOneMean(
  beta = NA_real_,
  n = NA_real_,
  meanH0 = 0,
  mean = 0.5,
  stDev = 1,
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

- meanH0:

  The mean under the null hypothesis. Defaults to 0.

- mean:

  The mean under the alternative hypothesis.

- stDev:

  The standard deviation.

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

An S3 class `designOneMean` object with three components:

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

  - `meanH0`: The mean under the null hypothesis.

  - `mean`: The mean under the alternative hypothesis.

  - `stDev`: The standard deviation.

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

  - `efficacyMean`: The efficacy boundaries on the mean scale.

  - `futilityMean`: The futility boundaries on the mean scale.

  - `numberOfSubjects`: The number of subjects.

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
(design1 <- getDesignOneMean(
  beta = 0.1, n = NA, meanH0 = 7, mean = 6, stDev = 2.5,
  kMax = 5, alpha = 0.025, typeAlphaSpending = "sfOF",
  typeBetaSpending = "sfP"))
#>                                                                               
#> Group-sequential design with 5 stages for one-sample mean                     
#> Mean under H0: 7, mean under H1: 6, standard deviation: 2.5                   
#> Overall power: 0.9016, overall alpha (1-sided): 0.025, attained alpha: 0.0191 
#> Drift parameter: 3.688, inflation factor: 1.287                               
#> Maximum information: 13.6, expected under H1: 8.74, expected under H0: 5.51   
#> Maximum # subjects: 85, expected under H1: 54.6, expected under H0: 34.4      
#> Alpha spending: Lan-DeMets O'Brien-Fleming, beta spending: Lan-DeMets Pocock  
#>                                                                               
#>                               Stage 1 Stage 2 Stage 3 Stage 4 Stage 5
#> Information rate              0.200   0.400   0.600   0.800   1.000  
#> Efficacy boundary (Z)         4.877   3.357   2.680   2.290   2.031  
#> Futility boundary (Z)         -0.246  0.500   1.075   1.557   2.031  
#> Cumulative rejection          0.0006  0.1528  0.5694  0.8280  0.9016 
#> Cumulative futility           0.0291  0.0515  0.0697  0.0851  0.0984 
#> Cumulative alpha spent        0.0000  0.0004  0.0038  0.0122  0.0250 
#> Number of subjects            17.0    34.0    51.0    68.0    85.0   
#> Efficacy boundary (mean)      4.043   5.561   6.062   6.306   6.449  
#> Futility boundary (mean)      7.149   6.785   6.624   6.528   6.449  
#> Efficacy boundary (p)         0.0000  0.0004  0.0037  0.0110  0.0211 
#> Futility boundary (p)         0.5970  0.3084  0.1411  0.0597  0.0211 
#> Information                   2.72    5.44    8.16    10.88   13.60  
#> Cumulative rejection under H0 0.0000  0.0004  0.0038  0.0114  0.0191 
#> Cumulative futility under H0  0.4030  0.7194  0.8826  0.9541  0.9809 

# Example 2: sample size calculation for one-sample t-test
(design2 <- getDesignOneMean(
  beta = 0.1, n = NA, meanH0 = 7, mean = 6, stDev = 2.5,
  normalApproximation = FALSE, alpha = 0.025))
#>                                                             
#> Fixed design for one-sample mean                            
#> Mean under H0: 7, mean under H1: 6, standard deviation: 2.5 
#> Overall power: 0.9016, overall alpha (1-sided): 0.025       
#> Drift parameter: 3.298, inflation factor: 1                 
#> Information: 10.88                                          
#> Number of subjects: 68                                      
#>                                                             
#>                                
#> Efficacy boundary (t)    1.996 
#> Efficacy boundary (mean) 6.395 
#> Efficacy boundary (p)    0.0250
```
