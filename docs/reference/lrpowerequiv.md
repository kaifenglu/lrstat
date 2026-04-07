# Power for Equivalence in Hazard Ratio

Obtains the power for equivalence in hazard ratio.

## Usage

``` r
lrpowerequiv(
  kMax = 1L,
  informationRates = NA_real_,
  criticalValues = NA_real_,
  alpha = 0.05,
  typeAlphaSpending = "sfOF",
  parameterAlphaSpending = NA_real_,
  userAlphaSpending = NA_real_,
  hazardRatioLower = NA_real_,
  hazardRatioUpper = NA_real_,
  allocationRatioPlanned = 1,
  accrualTime = 0L,
  accrualIntensity = NA_real_,
  piecewiseSurvivalTime = 0L,
  stratumFraction = 1L,
  lambda1 = NA_real_,
  lambda2 = NA_real_,
  gamma1 = 0L,
  gamma2 = 0L,
  accrualDuration = NA_real_,
  followupTime = NA_real_,
  fixedFollowup = 0L,
  typeOfComputation = "direct",
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

- hazardRatioLower:

  The lower equivalence limit of hazard ratio.

- hazardRatioUpper:

  The upper equivalence limit of hazard ratio.

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

- lambda1:

  A vector of hazard rates for the event in each analysis time interval
  by stratum for the active treatment group.

- lambda2:

  A vector of hazard rates for the event in each analysis time interval
  by stratum for the control group.

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

- typeOfComputation:

  The type of computation, either `"direct"` for the direct
  approximation method, or `"schoenfeld"` for the Schoenfeld method.
  Defaults to empty, which selects the Schoenfeld method under
  proportional hazards and ordinary log-rank test and the direct method
  otherwise.

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

An S3 class `lrpowerequiv` object with 4 components:

- `overallResults`: A data frame containing the following variables:

  - `overallReject`: The overall rejection probability.

  - `alpha`: The overall significance level.

  - `numberOfEvents`: The total number of events.

  - `numberOfDropouts`: The total number of dropouts.

  - `numbeOfSubjects`: The total number of subjects.

  - `studyDuration`: The total study duration.

  - `information`: The maximum information.

  - `expectedNumberOfEvents`: The expected number of events.

  - `expectedNumberOfDropouts`: The expected number of dropouts.

  - `expectedNumberOfSubjects`: The expected number of subjects.

  - `expectedStudyDuration`: The expected study duration.

  - `expectedInformation`: The expected information.

  - `kMax`: The number of stages.

  - `hazardRatioLower`: The lower equivalence limit of hazard ratio.

  - `hazardRatioUpper`: The upper equivalence limit of hazard ratio.

  - `accrualDuration`: The accrual duration.

  - `followupTime`: The follow-up time.

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

  - `analysisTime`: The average time since trial start.

  - `efficacyHRLower`: The efficacy boundaries on the hazard ratio scale
    for the one-sided null hypothesis at the lower equivalence limit.

  - `efficacyHRUpper`: The efficacy boundaries on the hazard ratio scale
    for the one-sided null hypothesis at the upper equivalence limit.

  - `efficacyP`: The efficacy bounds on the p-value scale for each of
    the two one-sided tests.

  - `information`: The cumulative information.

  - `HR`: The average hazard ratio.

- `settings`: A list containing the following input parameters:
  `typeAlphaSpending`, `parameterAlphaSpending`, `userAlphaSpending`,
  `allocationRatioPlanned`, `accrualTime`, `accuralIntensity`,
  `piecewiseSurvivalTime`, `stratumFraction`, `lambda1`, `lambda2`,
  `gamma1`, `gamma2`, `typeOfComputation`, and `spendingTime`.

- `byTreatmentCounts`: A list containing the following counts by
  treatment group:

  - `numberOfEvents1`: The number of events by stage for the treatment
    group.

  - `numberOfDropouts1`: The number of dropouts by stage for the
    treatment group.

  - `numberOfSubjects1`: The number of subjects by stage for the
    treatment group.

  - `numberOfEvents2`: The number of events by stage for the control
    group.

  - `numberOfDropouts2`: The number of dropouts by stage for the control
    group.

  - `numberOfSubjects2`: The number of subjects by stage for the control
    group.

  - `expectedNumberOfEvents1`: The expected number of events for the
    treatment group.

  - `expectedNumberOfDropouts1`: The expected number of dropouts for the
    treatment group.

  - `expectedNumberOfSubjects1`: The expected number of subjects for the
    treatment group.

  - `expectedNumberOfEvents2`: The expected number of events for control
    group.

  - `expectedNumberOfDropouts2`: The expected number of dropouts for the
    control group.

  - `expectedNumberOfSubjects2`: The expected number of subjects for the
    control group.

## See also

[`rmstat`](https://kaifenglu.github.io/lrstat/reference/rmstat.md)

## Author

Kaifeng Lu, <kaifenglu@gmail.com>

## Examples

``` r
lrpowerequiv(kMax = 2, informationRates = c(0.5, 1),
             alpha = 0.05, typeAlphaSpending = "sfOF",
             hazardRatioLower = 0.71, hazardRatioUpper = 1.4,
             allocationRatioPlanned = 1, accrualTime = seq(0, 8),
             accrualIntensity = 100/9*seq(1, 9),
             lambda1 = 0.0533,
             lambda2 = 0.0533,
             gamma1 = -log(1-0.05)/12,
             gamma2 = -log(1-0.05)/12, accrualDuration = 22,
             followupTime = 18, fixedFollowup = FALSE)
#>                                                                       
#> Group-sequential design with 2 stages for equivalence in hazard ratio 
#> Lower limit for hazard ratio: 0.71, upper limit for hazard ratio: 1.4 
#> Overall power: 1, overall alpha: 0.05                                 
#> Maximum # events: 1300.5, expected # events: 698.1                    
#> Maximum # dropouts: 104.3, expected # dropouts: 56                    
#> Maximum # subjects: 1800, expected # subjects: 1800                   
#> Maximum information: 325.14, expected information: 174.53             
#> Total study duration: 40, expected study duration: 23.6               
#> Accrual duration: 22, follow-up duration: 18, fixed follow-up: FALSE  
#> Allocation ratio: 1                                                   
#> Alpha spending: Lan-DeMets O'Brien-Fleming                            
#>                                                                       
#>                                        Stage 1 Stage 2
#> Information rate                       0.500   1.000  
#> Boundary for each 1-sided test (Z)     2.538   1.662  
#> Cumulative rejection                   0.9264  1.0000 
#> Cumulative alpha for each 1-sided test 0.0056  0.0500 
#> Cumulative alpha attained under H10    0.0056  0.0500 
#> Cumulative alpha attained under H20    0.0056  0.0500 
#> Number of events                       650.3   1300.5 
#> Number of dropouts                     52.1    104.3  
#> Number of subjects                     1800.0  1800.0 
#> Analysis time                          22.3    40.0   
#> Boundary for lower limit (HR)          0.866   0.779  
#> Boundary for upper limit (HR)          1.147   1.277  
#> Boundary for each 1-sided test (p)     0.0056  0.0482 
#> Information                            162.57  325.14 
#> HR                                     1.00    1.00   
```
