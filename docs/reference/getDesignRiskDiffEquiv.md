# Group Sequential Design for Equivalence in Two-Sample Risk Difference

Obtains the power given sample size or obtains the sample size given
power for a group sequential design for equivalence in two-sample risk
difference.

## Usage

``` r
getDesignRiskDiffEquiv(
  beta = NA_real_,
  n = NA_real_,
  riskDiffLower = NA_real_,
  riskDiffUpper = NA_real_,
  pi1 = NA_real_,
  pi2 = NA_real_,
  allocationRatioPlanned = 1,
  rounding = TRUE,
  kMax = 1L,
  informationRates = NA_real_,
  criticalValues = NA_real_,
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

- riskDiffLower:

  The lower equivalence limit of risk difference.

- riskDiffUpper:

  The upper equivalence limit of risk difference.

- pi1:

  The assumed probability for the active treatment group.

- pi2:

  The assumed probability for the control group.

- allocationRatioPlanned:

  Allocation ratio for the active treatment versus control. Defaults to
  1 for equal randomization.

- rounding:

  Whether to round up sample size. Defaults to 1 for sample size
  rounding.

- kMax:

  The maximum number of stages.

- informationRates:

  The information rates. Fixed prior to the trial. Defaults to
  `(1:kMax) / kMax` if left unspecified.

- criticalValues:

  Upper boundaries on the z-test statistic scale for stopping for
  efficacy.

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

An S3 class `designRiskDiffEquiv` object with three components:

- `overallResults`: A data frame containing the following variables:

  - `overallReject`: The overall rejection probability.

  - `alpha`: The significance level for each of the two one-sided tests.
    Defaults to 0.05.

  - `attainedAlphaH10`: The attained significance level under H10.

  - `attainedAlphaH20`: The attained significance level under H20.

  - `kMax`: The number of stages.

  - `information`: The maximum information.

  - `expectedInformationH1`: The expected information under H1.

  - `expectedInformationH10`: The expected information under H10.

  - `expectedInformationH20`: The expected information under H20.

  - `numberOfSubjects`: The maximum number of subjects.

  - `expectedNumberOfSubjectsH1`: The expected number of subjects under
    H1.

  - `expectedNumberOfSubjectsH10`: The expected number of subjects under
    H10.

  - `expectedNumberOfSubjectsH20`: The expected number of subjects under
    H20.

  - `riskDiffLower`: The lower equivalence limit of risk difference.

  - `riskDiffUpper`: The upper equivalence limit of risk difference.

  - `pi1`: The assumed probability for the active treatment group.

  - `pi2`: The assumed probability for the control group.

  - `riskDiff`: The risk difference.

- `byStageResults`: A data frame containing the following variables:

  - `informationRates`: The information rates.

  - `efficacyBounds`: The efficacy boundaries on the Z-scale for each of
    the two one-sided tests.

  - `rejectPerStage`: The probability for efficacy stopping.

  - `cumulativeRejection`: The cumulative probability for efficacy
    stopping.

  - `cumulativeAlphaSpent`: The cumulative alpha for each of the two
    one-sided tests.

  - `cumulativeAttainedAlphaH10`: The cumulative alpha attained under
    H10.

  - `cumulativeAttainedAlphaH20`: The cumulative alpha attained under
    H20.

  - `efficacyP`: The efficacy bounds on the p-value scale for each of
    the two one-sided tests.

  - `information`: The cumulative information.

  - `efficacyRiskDiffLower`: The efficacy boundaries on the risk
    difference scale for the one-sided null hypothesis on the lower
    equivalence limit.

  - `efficacyRiskDiffUpper`: The efficacy boundaries on the risk
    difference scale for the one-sided null hypothesis on the upper
    equivalence limit.

  - `numberOfSubjects`: The number of subjects.

- `settings`: A list containing the following input parameters:

  - `typeAlphaSpending`: The type of alpha spending.

  - `parameterAlphaSpending`: The parameter value for alpha spending.

  - `userAlphaSpending`: The user defined alpha spending.

  - `spendingTime`: The error spending time at each analysis.

  - `allocationRatioPlanned`: Allocation ratio for the active treatment
    versus control.

  - `rounding`: Whether to round up sample size.

## Author

Kaifeng Lu, <kaifenglu@gmail.com>

## Examples

``` r
(design1 <- getDesignRiskDiffEquiv(
  beta = 0.2, n = NA, riskDiffLower = -0.1,
  riskDiffUpper = 0.1, pi1 = 0.12, pi2 = 0.12,
  kMax = 3, alpha = 0.05, typeAlphaSpending = "sfOF"))
#>                                                                                         
#> Group-sequential design with 3 stages for equivalence in two-sample risk difference     
#> Lower limit for risk difference: -0.1, upper limit for risk difference: 0.1             
#> Proportion on treatment: 0.12, proportion on control: 0.12, risk difference: 0          
#> Overall power: 0.8007, overall alpha: 0.05, attained under H10: 0.05, under H20: 0.05   
#> Max information: 873.58, expected under H1: 810.5, under H10: 869.86, under H20: 869.86 
#> Max # subjects: 369, expected under H1: 342.4, under H10: 367.4, under H20: 367.4       
#> Allocation ratio: 1, variance of standardized test statistic:                           
#> Alpha spending: Lan-DeMets O'Brien-Fleming                                              
#>                                                                                         
#>                                        Stage 1 Stage 2 Stage 3
#> Information rate                       0.333   0.667   1.000  
#> Boundary for each 1-sided test (Z)     3.200   2.141   1.695  
#> Cumulative rejection                   0.0000  0.2166  0.8007 
#> Cumulative alpha for each 1-sided test 0.0007  0.0164  0.0500 
#> Cumulative alpha attained under H10    0.0000  0.0128  0.0500 
#> Cumulative alpha attained under H20    0.0000  0.0128  0.0500 
#> Number of subjects                     123.0   246.0   369.0  
#> Boundary for lower limit (risk diff)   0.088   -0.011  -0.043 
#> Boundary for upper limit (risk diff)   -0.088  0.011   0.043  
#> Boundary for each 1-sided test (p)     0.0007  0.0161  0.0451 
#> Information                            291.19  582.39  873.58 
```
