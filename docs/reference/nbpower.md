# Power for Negative Binomial Rate Ratio

Estimates the power for negative binomial rate ratio test.

## Usage

``` r
nbpower(
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
  rateRatioH0 = 1,
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
  studyDuration = NA_real_,
  nullVariance = FALSE
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
  for `"WT"`, \\rho\\ for `"sfKD"`, and \\\gamma\\ for `"sfHSD"`.

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

- rateRatioH0:

  Rate ratio under the null hypothesis.

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

- nullVariance:

  Whether to calculate the variance for log rate ratio under the null
  hypothesis.

## Value

An S3 class `nbpower` object with 4 components:

- `overallResults`: A data frame containing the following variables:

  - `overallReject`: The overall rejection probability.

  - `alpha`: The overall significance level.

  - `numberOfEvents`: The total number of events.

  - `numberOfDropouts`: The total number of dropouts.

  - `numbeOfSubjects`: The total number of subjects.

  - `exposure`: The total exposure.

  - `studyDuration`: The total study duration.

  - `information`: The maximum information.

  - `expectedNumberOfEvents`: The expected number of events.

  - `expectedNumberOfDropouts`: The expected number of dropouts.

  - `expectedNumberOfSubjects`: The expected number of subjects.

  - `expectedExposure`: The expected exposure.

  - `expectedStudyDuration`: The expected study duration.

  - `expectedInformation`: The expected information.

  - `accrualDuration`: The accrual duration.

  - `followupTime`: The follow-up duration.

  - `fixedFollowup`: Whether a fixed follow-up design is used.

  - `kMax`: The number of stages.

  - `rateRatioH0`: The rate ratio under the null hypothesis.

  - `rateRatio`: The rate ratio.

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

  - `exposure`: The exposure.

  - `analysisTime`: The average time since trial start.

  - `efficacyRateRatio`: The efficacy boundaries on the rate ratio
    scale.

  - `futilityRateRatio`: The futility boundaries on the rate ratio
    scale.

  - `efficacyP`: The efficacy boundaries on the p-value scale.

  - `futilityP`: The futility boundaries on the p-value scale.

  - `information`: The cumulative information.

  - `efficacyStopping`: Whether to allow efficacy stopping.

  - `futilityStopping`: Whether to allow futility stopping.

- `settings`: A list containing the following input parameters:
  `typeAlphaSpending`, `parameterAlphaSpending`, `userAlphaSpending`,
  `typeBetaSpending`, `parameterBetaSpending`, `allocationRatioPlanned`,
  `accrualTime`, `accuralIntensity`, `piecewiseSurvivalTime`, `kappa1`,
  `kappa2`, `lambda1`, `lambda2`, `gamma1`, `gamma2`, `spendingTime`,
  and `nullVariance`.

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

nbpower(kMax = 2, informationRates = c(0.5, 1),
        alpha = 0.025, typeAlphaSpending = "sfOF",
        accrualIntensity = 1956/1.25,
        stratumFraction = c(0.2, 0.8),
        kappa1 = 5, kappa2 = 5,
        lambda1 = c(0.7*0.125, 0.75*0.25),
        lambda2 = c(0.125, 0.25),
        gamma1 = 0, gamma2 = 0,
        accrualDuration = 1.25,
        followupTime = 2.75, fixedFollowup = FALSE,
        nullVariance = 1)
#>                                                                             
#> Group-sequential design with 2 stages for negative binomial rate ratio      
#> Rate ratio under H0: 1, rate ratio under H1: 0.74                           
#> Stratum fraction: 0.2 0.8                                                   
#> Event rate for treatment: 0.0875 0.1875, event rate for control: 0.125 0.25 
#> Dispersion for treatment: 5, dispersion for control: 5                      
#> Overall power: 0.7316, overall significance level (1-sided): 0.025          
#> Maximum # events: 1295.5, expected # events: 1163.1                         
#> Maximum # dropouts: 0, expected # dropouts: 0                               
#> Maximum # subjects: 1956, expected # subjects: 1956                         
#> Maximum exposure: 6601.5, expected exposure: 5926.4                         
#> Maximum information: 73.03, expected information: 68.23                     
#> Total study duration: 4, expected study duration: 3.7                       
#> Accrual duration: 1.2, follow-up duration: 2.8, fixed follow-up: FALSE      
#> Allocation ratio: 1, variance of standardized test statistic: under H0      
#> Alpha spending: Lan-DeMets O'Brien-Fleming, beta spending: None             
#>                                                                             
#>                                Stage 1 Stage 2
#> Information rate               0.500   1.000  
#> Efficacy boundary (Z)          2.963   1.969  
#> Cumulative rejection           0.1313  0.7316 
#> Cumulative alpha spent         0.0015  0.0250 
#> Number of events               286.5   1295.5 
#> Number of dropouts             0.0     0.0    
#> Number of subjects             1956.0  1956.0 
#> Exposure                       1459.7  6601.5 
#> Analysis time                  1.4     4.0    
#> Efficacy boundary (rate ratio) 0.615   0.795  
#> Efficacy boundary (p)          0.0015  0.0245 
#> Information                    36.51   73.03  

# Example 2: Fixed follow-up design

nbpower(kMax = 2, informationRates = c(0.5, 1),
        alpha = 0.025, typeAlphaSpending = "sfOF",
        accrualIntensity = 220/1.5,
        kappa1 = 3, kappa2 = 3,
        lambda1 = 0.5*8.4, lambda2 = 8.4,
        gamma1 = 0, gamma2 = 0,
        accrualDuration = 1.5,
        followupTime = 0.5, fixedFollowup = TRUE)
#>                                                                        
#> Group-sequential design with 2 stages for negative binomial rate ratio 
#> Rate ratio under H0: 1, rate ratio under H1: 0.5                       
#> Event rate for treatment: 4.2, event rate for control: 8.4             
#> Dispersion for treatment: 3, dispersion for control: 3                 
#> Overall power: 0.7997, overall significance level (1-sided): 0.025     
#> Maximum # events: 693, expected # events: 624                          
#> Maximum # dropouts: 0, expected # dropouts: 0                          
#> Maximum # subjects: 220, expected # subjects: 204.1                    
#> Maximum exposure: 110, expected exposure: 99                           
#> Maximum information: 16.38, expected information: 15.04                
#> Total study duration: 2, expected study duration: 1.8                  
#> Accrual duration: 1.5, follow-up duration: 0.5, fixed follow-up: TRUE  
#> Allocation ratio: 1, variance of standardized test statistic: under H1 
#> Alpha spending: Lan-DeMets O'Brien-Fleming, beta spending: None        
#>                                                                        
#>                                Stage 1 Stage 2
#> Information rate               0.500   1.000  
#> Efficacy boundary (Z)          2.963   1.969  
#> Cumulative rejection           0.1639  0.7997 
#> Cumulative alpha spent         0.0015  0.0250 
#> Number of events               272.0   693.0  
#> Number of dropouts             0.0     0.0    
#> Number of subjects             123.0   220.0  
#> Exposure                       43.2    110.0  
#> Analysis time                  0.8     2.0    
#> Efficacy boundary (rate ratio) 0.355   0.615  
#> Efficacy boundary (p)          0.0015  0.0245 
#> Information                    8.19    16.38  
```
