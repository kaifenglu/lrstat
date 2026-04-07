# Power for One-Sample Restricted Mean Survival Time

Estimates the power, stopping probabilities, and expected sample size in
a one-group survival design.

## Usage

``` r
rmpower1s(
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
  milestone = NA_real_,
  rmstH0 = NA_real_,
  accrualTime = 0L,
  accrualIntensity = NA_real_,
  piecewiseSurvivalTime = 0L,
  stratumFraction = 1L,
  lambda = NA_real_,
  gamma = 0L,
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

- typeBetaSpending:

  The type of beta spending. One of the following: "sfOF" for
  O'Brien-Fleming type spending function, "sfP" for Pocock type spending
  function, "sfKD" for Kim & DeMets spending function, "sfHSD" for
  Hwang, Shi & DeCani spending function, and "none" for no early
  futility stopping. Defaults to "none".

- parameterBetaSpending:

  The parameter value for the beta spending. Corresponds to \\\rho\\ for
  `"sfKD"`, and \\\gamma\\ for `"sfHSD"`.

- milestone:

  The milestone time at which to calculate the restricted mean survival
  time.

- rmstH0:

  The restricted mean survival time under the null hypothesis.

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

- lambda:

  A vector of hazard rates for the event in each analysis time interval
  by stratum under the alternative hypothesis.

- gamma:

  The hazard rate for exponential dropout or a vector of hazard rates
  for piecewise exponential dropout. Defaults to 0 for no dropout.

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

An S3 class `rmpower1s` object with 3 components:

- `overallResults`: A data frame containing the following variables:

  - `overallReject`: The overall rejection probability.

  - `alpha`: The overall significance level.

  - `numberOfEvents`: The total number of events.

  - `numbeOfSubjects`: The total number of subjects.

  - `numberOfMilestone`: The total number of subjects reaching
    milestone.

  - `studyDuration`: The total study duration.

  - `information`: The maximum information.

  - `expectedNumberOfEvents`: The expected number of events.

  - `expectedNumberOfSubjects`: The expected number of subjects.

  - `expectedNumberOfMilestone`: The expected number of subjects
    reaching milestone.

  - `expectedStudyDuration`: The expected study duration.

  - `expectedInformation`: The expected information.

  - `accrualDuration`: The accrual duration.

  - `followupTime`: The follow-up duration.

  - `fixedFollowup`: Whether a fixed follow-up design is used.

  - `kMax`: The number of stages.

  - `milestone`: The milestone time to calculate the restricted mean
    survival time.

  - `rmstH0`: The restricted mean survival time under the null
    hypothesis.

  - `rmst`: The restricted mean survival time under the alternative
    hypothesis.

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

  - `efficacyRmst`: The efficacy boundaries on the restricted mean
    survival time.

  - `futilityRmst`: The futility boundaries on the restricted mean
    survival time.

  - `efficacyP`: The efficacy boundaries on the p-value scale.

  - `futilityP`: The futility boundaries on the p-value scale.

  - `information`: The cumulative information.

  - `efficacyStopping`: Whether to allow efficacy stopping.

  - `futilityStopping`: Whether to allow futility stopping.

- `settings`: A list containing the following input parameters:
  `typeAlphaSpending`, `parameterAlphaSpending`, `userAlphaSpending`,
  `typeBetaSpending`, `parameterBetaSpending`, `accrualTime`,
  `accuralIntensity`, `piecewiseSurvivalTime`, `stratumFraction`,
  `lambda`, `gamma`, and `spendingTime`.

## See also

[`rmstat`](https://kaifenglu.github.io/lrstat/reference/rmstat.md)

## Author

Kaifeng Lu, <kaifenglu@gmail.com>

## Examples

``` r
rmpower1s(kMax = 2, informationRates = c(0.8, 1),
          alpha = 0.025, typeAlphaSpending = "sfOF",
          milestone = 18, rmstH0 = 10,
          accrualTime = seq(0, 8),
          accrualIntensity = 26/9*seq(1, 9),
          piecewiseSurvivalTime = c(0, 6),
          stratumFraction = c(0.2, 0.8),
          lambda = c(0.0533, 0.0309, 1.5*0.0533, 1.5*0.0309),
          gamma = -log(1-0.05)/12, accrualDuration = 22,
          followupTime = 18, fixedFollowup = FALSE)
#>                                                                                    
#> Group-sequential design with 2 stages for one-sample restricted mean survival time 
#> Milestone: 18, restricted mean survival time under H0: 10, under H1: 10.854        
#> Overall power: 0.7556, overall significance level (1-sided): 0.025                 
#> Maximum # events: 329.7, expected # events: 260.1                                  
#> Maximum # subjects: 468, expected # subjects: 468                                  
#> Maximum # milestone subjects: 166.5, expected # milestone subjects: 83.8           
#> Maximum information: 9.86, expected information: 8.76                              
#> Total study duration: 40, expected study duration: 30.8                            
#> Accrual duration: 22, follow-up duration: 18, fixed follow-up: FALSE               
#> Alpha spending: Lan-DeMets O'Brien-Fleming, beta spending: None                    
#>                                                                                    
#>                              Stage 1 Stage 2
#> Information rate             0.800   1.000  
#> Efficacy boundary (Z)        2.250   2.025  
#> Cumulative rejection         0.5585  0.7556 
#> Cumulative alpha spent       0.0122  0.0250 
#> Number of events             205.2   329.7  
#> Number of dropouts           14.1    26.0   
#> Number of subjects           468.0   468.0  
#> Number of milestone subjects 18.4    166.5  
#> Analysis time                23.5    40.0   
#> Efficacy boundary (rmst)     10.801  10.645 
#> Efficacy boundary (p)        0.0122  0.0214 
#> Information                  7.88    9.86   
```
