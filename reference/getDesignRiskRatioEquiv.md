# Group Sequential Design for Equivalence in Two-Sample Risk Ratio

Obtains the power given sample size or obtains the sample size given
power for a group sequential design for equivalence in two-sample risk
ratio.

## Usage

``` r
getDesignRiskRatioEquiv(
  beta = NA_real_,
  n = NA_real_,
  riskRatioLower = NA_real_,
  riskRatioUpper = NA_real_,
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

- riskRatioLower:

  The lower equivalence limit of risk ratio.

- riskRatioUpper:

  The upper equivalence limit of risk ratio.

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
  for `"WT"`, \\\rho\\ for `"sfKD"`, and \\\gamma\\ for `"sfHSD"`.

- userAlphaSpending:

  The user defined alpha spending. Cumulative alpha spent up to each
  stage.

- spendingTime:

  A vector of length `kMax` for the error spending time at each
  analysis. Defaults to missing, in which case, it is the same as
  `informationRates`.

## Value

An S3 class `designRiskRatioEquiv` object with three components:

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

  - `riskRatioLower`: The lower equivalence limit of risk ratio.

  - `riskRatioUpper`: The upper equivalence limit of risk ratio.

  - `pi1`: The assumed probability for the active treatment group.

  - `pi2`: The assumed probability for the control group.

  - `riskRatio`: The risk ratio.

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

  - `efficacyRiskRatioLower`: The efficacy boundaries on the risk ratio
    scale for the one-sided null hypothesis on the lower equivalence
    limit.

  - `efficacyRiskRatioUpper`: The efficacy boundaries on the risk ratio
    scale for the one-sided null hypothesis on the upper equivalence
    limit.

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

  - `rounding`: Whether to round up sample size.

## Author

Kaifeng Lu, <kaifenglu@gmail.com>

## Examples

``` r
(design1 <- getDesignRiskRatioEquiv(
  beta = 0.2, n = NA, riskRatioLower = 0.8,
  riskRatioUpper = 1.25, pi1 = 0.12, pi2 = 0.12,
  kMax = 3, alpha = 0.05, typeAlphaSpending = "sfOF"))
#>                                                                                          
#> Group-sequential design with 3 stages for equivalence in two-sample risk ratio           
#> Lower limit for risk ratio: 0.8, upper limit for risk ratio: 1.25                        
#> Proportion on treatment: 0.12, proportion on control: 0.12, risk ratio: 1                
#> Overall power: 0.8, overall alpha: 0.05, attained under H10: 0.05, under H20: 0.05       
#> Max information: 175.23, expected under H1: 162.63, under H10: 174.48, under H20: 174.48 
#> Max # subjects: 5140, expected under H1: 4770.6, under H10: 5118.2, under H20: 5118.2    
#> Allocation ratio: 1, variance of standardized test statistic:                            
#> Alpha spending: Lan-DeMets O'Brien-Fleming                                               
#>                                                                                          
#>                                        Stage 1 Stage 2 Stage 3
#> Information rate                       0.333   0.667   1.000  
#> Boundary for each 1-sided test (Z)     3.200   2.141   1.695  
#> Cumulative rejection                   0.0000  0.2157  0.8000 
#> Cumulative alpha for each 1-sided test 0.0007  0.0164  0.0500 
#> Cumulative alpha attained under H10    0.0000  0.0127  0.0500 
#> Cumulative alpha attained under H20    -0.0000 0.0127  0.0500 
#> Number of subjects                     1713.0  3427.0  5140.0 
#> Boundary for lower limit (risk ratio)  1.216   0.975   0.909  
#> Boundary for upper limit (risk ratio)  0.822   1.025   1.100  
#> Boundary for each 1-sided test (p)     0.0007  0.0161  0.0451 
#> Information                            58.40   116.83  175.23 
```
