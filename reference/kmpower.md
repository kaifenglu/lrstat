# Power for Difference in Milestone Survival Probabilities

Estimates the power for testing the difference in milestone survival
probabilities in a two-sample survival design.

## Usage

``` r
kmpower(
  kMax = 1L,
  informationRates = NA_real_,
  efficacyStopping = NA_integer_,
  futilityStopping = NA_integer_,
  criticalValues = NULL,
  alpha = 0.025,
  typeAlphaSpending = "sfOF",
  parameterAlphaSpending = NA_real_,
  userAlphaSpending = NA_real_,
  futilityBounds = NULL,
  futilityCP = NULL,
  futilitySurvDiff = NULL,
  typeBetaSpending = "none",
  parameterBetaSpending = NA_real_,
  milestone = NA_real_,
  survDiffH0 = 0,
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

- futilityCP:

  A vector of length `kMax - 1` for the futility bounds on the
  conditional power scale.

- futilitySurvDiff:

  A vector of length `kMax - 1` for the futility bounds on the milestone
  survival difference scale.

- typeBetaSpending:

  The type of beta spending. One of the following: `"sfOF"` for
  O'Brien-Fleming type spending function, `"sfP"` for Pocock type
  spending function, `"sfKD"` for Kim & DeMets spending function,
  `"sfHSD"` for Hwang, Shi & DeCani spending function, and `"none"` for
  no early futility stopping. Defaults to `"none"`.

- parameterBetaSpending:

  The parameter value for the beta spending. Corresponds to \\\rho\\ for
  `"sfKD"`, and \\\gamma\\ for `"sfHSD"`.

- milestone:

  The milestone time at which to calculate the survival probability.

- survDiffH0:

  The difference in milestone survival probabilities under the null
  hypothesis. Defaults to 0 for superiority test.

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

An S3 class `kmpower` object with 4 components:

- `overallResults`: A data frame containing the following variables:

  - `overallReject`: The overall rejection probability.

  - `alpha`: The overall significance level.

  - `numberOfEvents`: The total number of events.

  - `numberOfDropouts`: The total number of dropouts.

  - `numbeOfSubjects`: The total number of subjects.

  - `numberOfMilestone`: The total number of subjects reaching
    milestone.

  - `studyDuration`: The total study duration.

  - `information`: The maximum information.

  - `expectedNumberOfEvents`: The expected number of events.

  - `expectedNumberOfDropouts`: The expected number of dropouts.

  - `expectedNumberOfSubjects`: The expected number of subjects.

  - `expectedNumberOfMilestone`: The expected number of subjects
    reaching milestone.

  - `expectedStudyDuration`: The expected study duration.

  - `expectedInformation`: The expected information.

  - `accrualDuration`: The accrual duration.

  - `followupTime`: The follow-up duration.

  - `fixedFollowup`: Whether a fixed follow-up design is used.

  - `kMax`: The number of stages.

  - `milestone`: The milestone time relative to randomization.

  - `survDiffH0`: The difference in milestone survival probabilities
    under the null hypothesis.

  - `surv1`: The milestone survival probability for the treatment group.

  - `surv2`: The milestone survival probability for the control group.

  - `survDiff`: The difference in milestone survival probabilities,
    equal to `surv1 - surv2`.

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

  - `numberOfEvents`: The number of events.

  - `numberOfDropouts`: The number of dropouts.

  - `numberOfSubjects`: The number of subjects.

  - `numberOfMilestone`: The number of subjects reaching milestone.

  - `analysisTime`: The average time since trial start.

  - `efficacySurvDiff`: The efficacy boundaries on the survival
    difference scale.

  - `futilitySurvDiff`: The futility boundaries on the survival
    difference scale.

  - `efficacyP`: The efficacy boundaries on the p-value scale.

  - `futilityP`: The futility boundaries on the p-value scale.

  - `information`: The cumulative information.

  - `efficacyStopping`: Whether to allow efficacy stopping.

  - `futilityStopping`: Whether to allow futility stopping.

- `settings`: A list containing the following input parameters:
  `typeAlphaSpending`, `parameterAlphaSpending`, `userAlphaSpending`,
  `typeBetaSpending`, `parameterBetaSpending`, `allocationRatioPlanned`,
  `accrualTime`, `accuralIntensity`, `piecewiseSurvivalTime`,
  `stratumFraction`, `lambda1`, `lambda2`, `gamma1`, `gamma2`, and
  `spendingTime`.

- `byTreatmentCounts`: A list containing the following counts by
  treatment group:

  - `numberOfEvents1`: The number of events by stage for the treatment
    group.

  - `numberOfDropouts1`: The number of dropouts by stage for the
    treatment group.

  - `numberOfSubjects1`: The number of subjects by stage for the
    treatment group.

  - `numberOfMilestone1`: The number of subjects reaching milestone by
    stage for the active treatment group.

  - `numberOfEvents2`: The number of events by stage for the control
    group.

  - `numberOfDropouts2`: The number of dropouts by stage for the control
    group.

  - `numberOfSubjects2`: The number of subjects by stage for the control
    group.

  - `numberOfMilestone2`: The number of subjects reaching milestone by
    stage for the control group.

  - `expectedNumberOfEvents1`: The expected number of events for the
    treatment group.

  - `expectedNumberOfDropouts1`: The expected number of dropouts for the
    active treatment group.

  - `expectedNumberOfSubjects1`: The expected number of subjects for the
    active treatment group.

  - `expectedNumberOfMilestone1`: The expected number of subjects
    reaching milestone for the active treatment group.

  - `expectedNumberOfEvents2`: The expected number of events for control
    group.

  - `expectedNumberOfDropouts2`: The expected number of dropouts for the
    control group.

  - `expectedNumberOfSubjects2`: The expected number of subjects for the
    control group.

  - `expectedNumberOfMilestone2`: The expected number of subjects
    reaching milestone for the control group.

## Author

Kaifeng Lu, <kaifenglu@gmail.com>

## Examples

``` r
# Piecewise accrual, piecewise exponential survival, and 5% dropout by
# the end of 1 year.

kmpower(kMax = 2, informationRates = c(0.8, 1),
        alpha = 0.025, typeAlphaSpending = "sfOF",
        milestone = 18,
        allocationRatioPlanned = 1, accrualTime = seq(0, 8),
        accrualIntensity = 26/9*seq(1, 9),
        piecewiseSurvivalTime = c(0, 6),
        stratumFraction = c(0.2, 0.8),
        lambda1 = c(0.0533, 0.0309, 1.5*0.0533, 1.5*0.0309),
        lambda2 = c(0.0533, 0.0533, 1.5*0.0533, 1.5*0.0533),
        gamma1 = -log(1-0.05)/12,
        gamma2 = -log(1-0.05)/12, accrualDuration = 22,
        followupTime = 18, fixedFollowup = FALSE)
#>                                                                            
#> Group-sequential design with 2 stages for difference in milestone survival 
#> Milestone: 18, survival difference under H0: 0                             
#> Milestone survival on treatment: 0.384, on control: 0.266                  
#> Overall power: 0.7634, overall significance level (1-sided): 0.025         
#> Maximum # events: 356.4, expected # events: 325.3                          
#> Maximum # subjects: 468, expected # subjects: 468                          
#> Maximum information: 527.17, expected information: 467.35                  
#> Total study duration: 40, expected study duration: 34.9                    
#> Accrual duration: 22, follow-up duration: 18, fixed follow-up: FALSE       
#> Allocation ratio: 1                                                        
#> Alpha spending: Lan-DeMets O'Brien-Fleming, beta spending: None            
#>                                                                            
#>                               Stage 1 Stage 2
#> Information rate              0.800   1.000  
#> Efficacy boundary (Z)         2.250   2.025  
#> Cumulative rejection          0.5674  0.7634 
#> Cumulative alpha spent        0.0122  0.0250 
#> Number of events              301.5   356.4  
#> Number of dropouts            19.8    24.2   
#> Number of subjects            468.0   468.0  
#> Number of milestone subjects  71.1    140.9  
#> Analysis time                 31.1    40.0   
#> Efficacy boundary (surv diff) 0.110   0.088  
#> Efficacy boundary (p)         0.0122  0.0214 
#> Information                   421.74  527.17 
```
