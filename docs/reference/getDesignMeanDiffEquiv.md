# Group Sequential Design for Equivalence in Two-Sample Mean Difference

Obtains the power given sample size or obtains the sample size given
power for a group sequential design for equivalence in two-sample mean
difference.

## Usage

``` r
getDesignMeanDiffEquiv(
  beta = NA_real_,
  n = NA_real_,
  meanDiffLower = NA_real_,
  meanDiffUpper = NA_real_,
  meanDiff = 0,
  stDev = 1,
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

- meanDiffLower:

  The lower equivalence limit of mean difference.

- meanDiffUpper:

  The upper equivalence limit of mean difference.

- meanDiff:

  The mean difference under the alternative hypothesis.

- stDev:

  The standard deviation.

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
  for `"WT"`, \\rho\\ for `"sfKD"`, and \\\gamma\\ for `"sfHSD"`.

- userAlphaSpending:

  The user defined alpha spending. Cumulative alpha spent up to each
  stage.

- spendingTime:

  A vector of length `kMax` for the error spending time at each
  analysis. Defaults to missing, in which case, it is the same as
  `informationRates`.

## Value

An S3 class `designMeanDiffEquiv` object with three components:

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

  - `meanDiffLower`: The lower equivalence limit of mean difference.

  - `meanDiffUpper`: The upper equivalence limit of mean difference.

  - `meanDiff`: The mean difference under the alternative hypothesis.

  - `stDev`: The standard deviation.

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

  - `efficacyMeanDiffLower`: The efficacy boundaries on the mean
    difference scale for the one-sided null hypothesis on the lower
    equivalence limit.

  - `efficacyMeanDiffUpper`: The efficacy boundaries on the mean
    difference scale for the one-sided null hypothesis on the upper
    equivalence limit.

  - `efficacyP`: The efficacy bounds on the p-value scale for each of
    the two one-sided tests.

  - `information`: The cumulative information.

  - `numberOfSubjects`: The number of subjects.

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
(design1 <- getDesignMeanDiffEquiv(
  beta = 0.1, n = NA, meanDiffLower = -1.3, meanDiffUpper = 1.3,
  meanDiff = 0, stDev = 2.2,
  kMax = 4, alpha = 0.05, typeAlphaSpending = "sfOF"))
#>                                                                                     
#> Group-sequential design with 4 stages for equivalence in two-sample mean difference 
#> Lower limit for mean difference: -1.3, upper limit for mean difference: 1.3         
#> Mean difference under H1: 0, standard deviation: 2.2                                
#> Overall power: 0.9024, overall alpha: 0.05, attained alpha: 0.05                    
#> Maximum information: 6.61, expected under H1: 5.57, expected under H0: 6.57         
#> Maximum # subjects: 128, expected under H1: 107.9, expected under H0: 127.2         
#> Allocation ratio: 1                                                                 
#> Alpha spending: Lan-DeMets O'Brien-Fleming                                          
#>                                                                                     
#>                                        Stage 1 Stage 2 Stage 3 Stage 4
#> Information rate                       0.250   0.500   0.750   1.000  
#> Boundary for each 1-sided test (Z)     3.750   2.540   2.016   1.720  
#> Cumulative rejection                   0.0000  0.0002  0.6283  0.9024 
#> Cumulative alpha for each 1-sided test 0.0001  0.0056  0.0236  0.0500 
#> Cumulative alpha attained under H0     0.0000  0.0000  0.0235  0.0500 
#> Number of subjects                     32.0    64.0    96.0    128.0  
#> Boundary for lower limit (mean diff)   1.616   0.097   -0.395  -0.631 
#> Boundary for upper limit (mean diff)   -1.616  -0.097  0.395   0.631  
#> Boundary for each 1-sided test (p)     0.0001  0.0055  0.0219  0.0427 
#> Information                            1.65    3.31    4.96    6.61   

# Example 2: sample size calculation for t-test
(design2 <- getDesignMeanDiffEquiv(
  beta = 0.1, n = NA, meanDiffLower = -1.3, meanDiffUpper = 1.3,
  meanDiff = 0, stDev = 2.2,
  normalApproximation = FALSE, alpha = 0.05))
#>                                                                             
#> Fixed design for equivalence in two-sample mean difference                  
#> Lower limit for mean difference: -1.3, upper limit for mean difference: 1.3 
#> Mean difference under H1: 0, standard deviation: 2.2                        
#> Overall power: 0.9018, overall alpha: 0.05, attained alpha: 0.05            
#> Information: 6.51                                                           
#> Number of subjects: 126                                                     
#> Allocation ratio: 1                                                         
#>                                                                             
#>                                            
#> Boundary for each 1-sided test (t)   1.657 
#> Boundary for lower limit (mean diff) -0.650
#> Boundary for upper limit (mean diff) 0.650 
#> Boundary for each 1-sided test (p)   0.0500
```
