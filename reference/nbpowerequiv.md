# Power for Equivalence in Negative Binomial Rate Ratio

Obtains the power for equivalence in negative binomial rate ratio.

## Usage

``` r
nbpowerequiv(
  kMax = 1L,
  informationRates = NA_real_,
  criticalValues = NULL,
  alpha = 0.05,
  typeAlphaSpending = "sfOF",
  parameterAlphaSpending = NA_real_,
  userAlphaSpending = NA_real_,
  rateRatioLower = NA_real_,
  rateRatioUpper = NA_real_,
  allocationRatioPlanned = 1,
  accrualTime = 0L,
  accrualIntensity = NA_real_,
  piecewiseSurvivalTime = 0L,
  stratumFraction = 1L,
  kappa1 = NA_real_,
  kappa2 = NA_real_,
  lambda1 = NA_real_,
  lambda2 = NA_real_,
  gamma1 = 0L,
  gamma2 = 0L,
  accrualDuration = NA_real_,
  followupTime = NA_real_,
  fixedFollowup = FALSE,
  spendingTime = NA_real_,
  studyDuration = NA_real_
)
```

## Arguments

- kMax:

  The maximum number of stages.

- informationRates:

  The information rates. Defaults to `(1:kMax) / kMax` if left
  unspecified.

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

- rateRatioLower:

  The lower equivalence limit of rate ratio.

- rateRatioUpper:

  The upper equivalence limit of rate ratio.

- allocationRatioPlanned:

  Allocation ratio for the active treatment versus control. Defaults to
  1 for equal randomization.

- accrualTime:

  A vector that specifies the starting time of piecewise Poisson
  enrollment time intervals. Must start with 0, e.g., `c(0, 3)` breaks
  the time axis into 2 accrual intervals: \\\[0, 3)\\ and \\\[3,
  \infty)\\.

- accrualIntensity:

  A vector of accrual intensities. One for each accrual time interval.

- piecewiseSurvivalTime:

  A vector that specifies the starting time of piecewise exponential
  survival time intervals. Must start with 0, e.g., `c(0, 6)` breaks the
  time axis into 2 event intervals: \\\[0, 6)\\ and \\\[6, \infty)\\.
  Defaults to 0 for exponential distribution.

- stratumFraction:

  A vector of stratum fractions that sum to 1. Defaults to 1 for no
  stratification.

- kappa1:

  The dispersion parameter (reciprocal of the shape parameter of the
  gamma mixing distribution) for the active treatment group by stratum.

- kappa2:

  The dispersion parameter (reciprocal of the shape parameter of the
  gamma mixing distribution) for the control group by stratum.

- lambda1:

  The rate parameter of the negative binomial distribution for the
  active treatment group by stratum.

- lambda2:

  The rate parameter of the negative binomial distribution for the
  control group by stratum.

- gamma1:

  The hazard rate for exponential dropout, a vector of hazard rates for
  piecewise exponential dropout applicable for all strata, or a vector
  of hazard rates for dropout in each analysis time interval by stratum
  for the active treatment group.

- gamma2:

  The hazard rate for exponential dropout, a vector of hazard rates for
  piecewise exponential dropout applicable for all strata, or a vector
  of hazard rates for dropout in each analysis time interval by stratum
  for the control group.

- accrualDuration:

  Duration of the enrollment period.

- followupTime:

  Follow-up time for the last enrolled subject.

- fixedFollowup:

  Whether a fixed follow-up design is used. Defaults to `FALSE` for
  variable follow-up.

- spendingTime:

  A vector of length `kMax` for the error spending time at each
  analysis. Defaults to missing, in which case, it is the same as
  `informationRates`.

- studyDuration:

  Study duration for fixed follow-up design. Defaults to missing, which
  is to be replaced with the sum of `accrualDuration` and
  `followupTime`. If provided, the value is allowed to be less than the
  sum of `accrualDuration` and `followupTime`.

## Value

An S3 class `nbpowerequiv` object with 4 components:

- `overallResults`: A data frame containing the following variables:

  - `overallReject`: The overall rejection probability.

  - `alpha`: The overall significance level.

  - `numberOfEvents`: The total number of events.

  - `numbeOfSubjects`: The total number of subjects.

  - `exposure`: The total exposure.

  - `studyDuration`: The total study duration.

  - `information`: The maximum information.

  - `expectedNumberOfEvents`: The expected number of events.

  - `expectedNumberOfSubjects`: The expected number of subjects.

  - `expectedExposure`: The expected exposure.

  - `expectedStudyDuration`: The expected study duration.

  - `expectedInformation`: The expected information.

  - `kMax`: The number of stages.

  - `rateRatioLower`: The lower equivalence limit of rate ratio.

  - `rateRatioUpper`: The upper equivalence limit of rate ratio.

  - `rateRatio`: The rate ratio.

  - `accrualDuration`: The accrual duration.

  - `followupTime`: The follow-up duration.

  - `fixedFollowup`: Whether a fixed follow-up design is used.

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
    `H10`.

  - `cumulativeAttainedAlphaH20`: The cumulative alpha attained under
    `H20`.

  - `numberOfEvents`: The number of events.

  - `numberOfDropouts`: The number of dropouts.

  - `numberOfSubjects`: The number of subjects.

  - `exposure`: The exposure.

  - `analysisTime`: The average time since trial start.

  - `efficacyRateRatioLower`: The efficacy boundaries on the rate ratio
    scale for the one-sided null hypothesis at the lower equivalence
    limit.

  - `efficacyRateRatioUpper`: The efficacy boundaries on the rate ratio
    scale for the one-sided null hypothesis at the upper equivalence
    limit.

  - `efficacyP`: The efficacy bounds on the p-value scale for each of
    the two one-sided tests.

  - `information`: The cumulative information.

- `settings`: A list containing the following input parameters:
  `typeAlphaSpending`, `parameterAlphaSpending`, `userAlphaSpending`,
  `allocationRatioPlanned`, `accrualTime`, `accuralIntensity`,
  `piecewiseSurvivalTime`, `stratumFraction`, `kappa1`, `kappa2`,
  `lambda1`, `lambda2`, `gamma1`, `gamma2`, `spendingTime`.

- `byTreatmentCounts`: A list containing the following counts by
  treatment group:

  - `numberOfEvents1`: The number of events by stage for the treatment
    group.

  - `numberOfDropouts1`: The number of dropouts by stage for the
    treatment group.

  - `numberOfSubjects1`: The number of subjects by stage for the
    treatment group.

  - `exposure1`: The exposure by stage for the treatment group.

  - `numberOfEvents2`: The number of events by stage for the control
    group.

  - `numberOfDropouts2`: The number of dropouts by stage for the control
    group.

  - `numberOfSubjects2`: The number of subjects by stage for the control
    group.

  - `exposure2`: The exposure by stage for the control group.

  - `expectedNumberOfEvents1`: The expected number of events for the
    treatment group.

  - `expectedNumberOfDropouts1`: The expected number of dropouts for the
    treatment group.

  - `expectedNumberOfSubjects1`: The expected number of subjects for the
    treatment group.

  - `expectedExposure1`: The expected exposure for the treatment group.

  - `expectedNumberOfEvents2`: The expected number of events for control
    group.

  - `expectedNumberOfDropouts2`: The expected number of dropouts for the
    control group.

  - `expectedNumberOfSubjects2`: The expected number of subjects for the
    control group.

  - `expectedExposure2`: The expected exposure for the control group.

## See also

[`nbstat`](https://kaifenglu.github.io/lrstat/reference/nbstat.md)

## Author

Kaifeng Lu, <kaifenglu@gmail.com>

## Examples

``` r
# Example 1: Variable follow-up design
nbpowerequiv(kMax = 2, informationRates = c(0.5, 1),
             alpha = 0.05, typeAlphaSpending = "sfOF",
             rateRatioLower = 2/3, rateRatioUpper = 3/2,
             accrualIntensity = 1956/1.25,
             kappa1 = 5, kappa2 = 5,
             lambda1 = 0.125, lambda2 = 0.125,
             gamma1 = 0, gamma2 = 0,
             accrualDuration = 1.25,
             followupTime = 2.75, fixedFollowup = FALSE)
#>                                                                                       
#> Group-sequential design with 2 stages for equivalence in negative binomial rate ratio 
#> Lower limit for rate ratio: 0.667, upper limit for rate ratio: 1.5, rate ratio: 1     
#> Event rate for treatment: 0.125, event rate for control: 0.125                        
#> Dispersion for treatment: 5, dispersion for control: 5                                
#> Overall power: 0.8995, overall alpha: 0.05                                            
#> Maximum # events: 825.2, expected # events: 825.2                                     
#> Maximum # subjects: 1956, expected # subjects: 1956                                   
#> Maximum exposure: 6601.5, expected exposure: 6601.5                                   
#> Maximum information: 66.18, expected information: 66.18                               
#> Total study duration: 4, expected study duration: 4                                   
#> Accrual duration: 1.2, follow-up duration: 2.8, fixed follow-up: FALSE                
#> Allocation ratio: 1                                                                   
#> Alpha spending: Lan-DeMets O'Brien-Fleming                                            
#>                                                                                       
#>                                        Stage 1 Stage 2
#> Information rate                       0.500   1.000  
#> Boundary for each 1-sided test (Z)     2.538   1.662  
#> Cumulative rejection                   -0.0000 0.8995 
#> Cumulative alpha for each 1-sided test 0.0056  0.0500 
#> Cumulative alpha attained under H10    -0.9833 -0.9500
#> Cumulative alpha attained under H20    0.0000  0.0500 
#> Number of events                       213.1   825.2  
#> Number of dropouts                     0.0     0.0    
#> Number of subjects                     1956.0  1956.0 
#> Exposure                               1705.2  6601.5 
#> Analysis time                          1.5     4.0    
#> Boundary for lower limit (rate ratio)  1.036   0.818  
#> Boundary for upper limit (rate ratio)  0.965   1.223  
#> Boundary for each 1-sided test (p)     0.0056  0.0482 
#> Information                            33.09   66.18  

# Example 2: Fixed follow-up design
nbpowerequiv(kMax = 2, informationRates = c(0.5, 1),
             alpha = 0.05, typeAlphaSpending = "sfOF",
             rateRatioLower = 0.5, rateRatioUpper = 2,
             accrualIntensity = 220/1.5,
             stratumFraction = c(0.2, 0.8),
             kappa1 = 3, kappa2 = 3,
             lambda1 = c(8.4, 10.2),
             lambda2 = c(8.0, 11.5),
             gamma1 = 0, gamma2 = 0,
             accrualDuration = 1.5,
             followupTime = 0.5, fixedFollowup = TRUE)
#>                                                                                       
#> Group-sequential design with 2 stages for equivalence in negative binomial rate ratio 
#> Lower limit for rate ratio: 0.5, upper limit for rate ratio: 2, rate ratio: 0.917     
#> Stratum fraction: 0.2 0.8                                                             
#> Event rate for treatment: 8.4 10.2, event rate for control: 8 11.5                    
#> Dispersion for treatment: 3, dispersion for control: 3                                
#> Overall power: 0.7485, overall alpha: 0.05                                            
#> Maximum # events: 1135.2, expected # events: 1135.2                                   
#> Maximum # subjects: 220, expected # subjects: 220                                     
#> Maximum exposure: 110, expected exposure: 110                                         
#> Maximum information: 17.2, expected information: 17.2                                 
#> Total study duration: 2, expected study duration: 2                                   
#> Accrual duration: 1.5, follow-up duration: 0.5, fixed follow-up: TRUE                 
#> Allocation ratio: 1                                                                   
#> Alpha spending: Lan-DeMets O'Brien-Fleming                                            
#>                                                                                       
#>                                        Stage 1 Stage 2
#> Information rate                       0.500   1.000  
#> Boundary for each 1-sided test (Z)     2.538   1.662  
#> Cumulative rejection                   -0.0000 0.7485 
#> Cumulative alpha for each 1-sided test 0.0056  0.0500 
#> Cumulative alpha attained under H10    -0.9367 -0.9500
#> Cumulative alpha attained under H20    0.0000  0.0500 
#> Number of events                       427.2   1135.2 
#> Number of dropouts                     0.0     0.0    
#> Number of subjects                     119.5   220.0  
#> Exposure                               41.4    110.0  
#> Analysis time                          0.8     2.0    
#> Boundary for lower limit (rate ratio)  1.188   0.746  
#> Boundary for upper limit (rate ratio)  0.842   1.340  
#> Boundary for each 1-sided test (p)     0.0056  0.0482 
#> Information                            8.60    17.20  
```
