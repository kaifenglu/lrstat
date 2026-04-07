# Group Sequential Design for Equivalence in Two-Sample Mean Ratio

Obtains the power given sample size or obtains the sample size given
power for a group sequential design for equivalence in two-sample mean
ratio.

## Usage

``` r
getDesignMeanRatioEquiv(
  beta = NA_real_,
  n = NA_real_,
  meanRatioLower = NA_real_,
  meanRatioUpper = NA_real_,
  meanRatio = 1,
  CV = 1,
  allocationRatioPlanned = 1,
  normalApproximation = TRUE,
  rounding = TRUE,
  kMax = 1L,
  informationRates = NA_real_,
  alpha = 0.05,
  typeAlphaSpending = "sfOF",
  parameterAlphaSpending = NA_real_,
  userAlphaSpending = NA_real_,
  spendingTime = NA_real_
)
```

## Arguments

- beta:

  The type II error.

- n:

  The total sample size.

- meanRatioLower:

  The lower equivalence limit of mean ratio.

- meanRatioUpper:

  The upper equivalence limit of mean ratio.

- meanRatio:

  The mean ratio under the alternative hypothesis.

- CV:

  The coefficient of variation.

- allocationRatioPlanned:

  Allocation ratio for the active treatment versus control. Defaults to
  1 for equal randomization.

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

- alpha:

  The significance level for each of the two one-sided tests. Defaults
  to 0.05.

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

- spendingTime:

  A vector of length `kMax` for the error spending time at each
  analysis. Defaults to missing, in which case, it is the same as
  `informationRates`.

## Value

An S3 class `designMeanRatioEquiv` object with three components:

- `overallResults`: A data frame containing the following variables:

  - `overallReject`: The overall rejection probability.

  - `alpha`: The significance level for each of the two one-sided tests.
    Defaults to 0.05.

  - `attainedAlpha`: The attained significance level.

  - `kMax`: The number of stages.

  - `information`: The maximum information.

  - `expectedInformationH1`: The expected information under H1.

  - `expectedInformationH0`: The expected information under H0.

  - `numberOfSubjects`: The maximum number of subjects.

  - `expectedNumberOfSubjectsH1`: The expected number of subjects under
    H1.

  - `expectedNumberOfSubjectsH0`: The expected number of subjects under
    H0.

  - `meanRatioLower`: The lower equivalence limit of mean ratio.

  - `meanRatioUpper`: The upper equivalence limit of mean ratio.

  - `meanRatio`: The mean ratio under the alternative hypothesis.

  - `CV`: The coefficient of variation.

- `byStageResults`: A data frame containing the following variables:

  - `informationRates`: The information rates.

  - `efficacyBounds`: The efficacy boundaries on the Z-scale for each of
    the two one-sided tests.

  - `rejectPerStage`: The probability for efficacy stopping.

  - `cumulativeRejection`: The cumulative probability for efficacy
    stopping.

  - `cumulativeAlphaSpent`: The cumulative alpha for each of the two
    one-sided tests.

  - `cumulativeAttainedAlpha`: The cumulative probability for efficacy
    stopping under H0.

  - `efficacyP`: The efficacy bounds on the p-value scale for each of
    the two one-sided tests.

  - `information`: The cumulative information.

  - `numberOfSubjects`: The number of subjects.

  - `efficacyMeanRatioLower`: The efficacy boundaries on the mean ratio
    scale for the one-sided null hypothesis on the lower equivalence
    limit.

  - `efficacyMeanRatioUpper`: The efficacy boundaries on the mean ratio
    scale for the one-sided null hypothesis on the upper equivalence
    limit.

- `settings`: A list containing the following input parameters:

  - `typeAlphaSpending`: The type of alpha spending.

  - `parameterAlphaSpending`: The parameter value for alpha spending.

  - `userAlphaSpending`: The user defined alpha spending.

  - `spendingTime`: The error spending time at each analysis.

  - `allocationRatioPlanned`: Allocation ratio for the active treatment
    versus control.

  - `normalApproximation`: The type of computation of the p-values. If
    `TRUE`, the variance is assumed to be known, otherwise the
    calculations are performed with the t distribution. The exact
    calculation using the t distribution is only implemented for the
    fixed design.

  - `rounding`: Whether to round up sample size.

## Author

Kaifeng Lu, <kaifenglu@gmail.com>

## Examples

``` r
# Example 1: group sequential trial power calculation
(design1 <- getDesignMeanRatioEquiv(
  beta = 0.1, n = NA, meanRatioLower = 0.8, meanRatioUpper = 1.25,
  meanRatio = 1, CV = 0.35,
  kMax = 4, alpha = 0.05, typeAlphaSpending = "sfOF"))
#>                                                                                   
#> Group-sequential design with 4 stages for equivalence in two-sample mean ratio    
#> Lower limit for mean ratio: 0.8, upper limit for mean ratio: 1.25                 
#> Mean ratio under H1: 1, coefficient of variation: 0.35                            
#> Overall power: 0.9033, overall alpha (1-sided): 0.05, attained alpha: 0.05        
#> Maximum information: 224.99, expected under H1: 189.53, expected under H0: 223.67 
#> Maximum # subjects: 104, expected under H1: 87.6, expected under H0: 103.4        
#> Allocation ratio: 1                                                               
#> Alpha spending: Lan-DeMets O'Brien-Fleming                                        
#>                                                                                   
#>                                        Stage 1 Stage 2 Stage 3 Stage 4
#> Information rate                       0.250   0.500   0.750   1.000  
#> Boundary for each 1-sided test (Z)     3.750   2.540   2.016   1.720  
#> Cumulative rejection                   0.0000  0.0002  0.6303  0.9033 
#> Cumulative alpha for each 1-sided test 0.0001  0.0056  0.0236  0.0500 
#> Cumulative alpha spent                 0.0000  0.0000  0.0236  0.0500 
#> Number of subjects                     26.0    52.0    78.0    104.0  
#> Boundary for lower limit (mean ratio)  1.319   1.016   0.934   0.897  
#> Boundary for upper limit (mean ratio)  0.758   0.984   1.070   1.115  
#> Boundary for each 1-sided test (p)     0.0001  0.0055  0.0219  0.0427 
#> Information                            56.25   112.50  168.75  224.99 

# Example 2: sample size calculation for t-test
(design2 <- getDesignMeanRatioEquiv(
  beta = 0.1, n = NA, meanRatioLower = 0.8, meanRatioUpper = 1.25,
  meanRatio = 1, CV = 0.35,
  normalApproximation = FALSE, alpha = 0.05))
#>                                                                            
#> Fixed design for equivalence in two-sample mean ratio                      
#> Lower limit for mean ratio: 0.8, upper limit for mean ratio: 1.25          
#> Mean ratio under H1: 1, coefficient of variation: 0.35                     
#> Overall power: 0.9005, overall alpha (1-sided): 0.05, attained alpha: 0.05 
#> Information: 220.67                                                        
#> Number of subjects: 102                                                    
#> Allocation ratio: 1                                                        
#>                                                                            
#>                                             
#> Boundary for each 1-sided test (t)    1.660 
#> Boundary for lower limit (mean ratio) 0.895 
#> Boundary for upper limit (mean ratio) 1.118 
#> Boundary for each 1-sided test (p)    0.0500
```
