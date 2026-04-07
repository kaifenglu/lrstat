# Power for One-Sample Negative Binomial Rate

Estimates the power, stopping probabilities, and expected sample size in
a one-group negative binomial design.

## Usage

``` r
nbpower1s(
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
  lambdaH0 = NA_real_,
  accrualTime = 0L,
  accrualIntensity = NA_real_,
  piecewiseSurvivalTime = 0L,
  stratumFraction = 1L,
  kappa = NA_real_,
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

- lambdaH0:

  The rate parameter of the negative binomial distribution under the
  null hypothesis.

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

- kappa:

  The dispersion parameter (reciprocal of the shape parameter of the
  gamma mixing distribution) of the negative binomial distribution by
  stratum.

- lambda:

  The rate parameter of the negative binomial distribution under the
  alternative hypothesis by stratum.

- gamma:

  The hazard rate for exponential dropout or a vector of hazard rates
  for piecewise exponential dropout by stratum. Defaults to 0 for no
  dropout.

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

An S3 class `nbpower1s` object with 3 components:

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

  - `lambdaH0`: The rate parameter of the negative binomial distribution
    under the null hypothesis.

  - `lambda`: The overall rate parameter of the negative binomial
    distribution under the alternative hypothesis.

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

  - `efficacyRate`: The efficacy boundaries on the rate scale.

  - `futilityRate`: The futility boundaries on the rate scale.

  - `efficacyP`: The efficacy boundaries on the p-value scale.

  - `futilityP`: The futility boundaries on the p-value scale.

  - `information`: The cumulative information.

  - `efficacyStopping`: Whether to allow efficacy stopping.

  - `futilityStopping`: Whether to allow futility stopping.

- `settings`: A list containing the following input parameters:
  `typeAlphaSpending`, `parameterAlphaSpending`, `userAlphaSpending`,
  `typeBetaSpending`, `parameterBetaSpending`, `accrualTime`,
  `accuralIntensity`, `piecewiseSurvivalTime`, `stratumFraction`,
  `kappa`, `lambda`, `gamma`, and `spendingTime`.

## See also

[`nbstat`](https://kaifenglu.github.io/lrstat/reference/nbstat.md)

## Author

Kaifeng Lu, <kaifenglu@gmail.com>

## Examples

``` r
# Example 1: Variable follow-up design

nbpower1s(kMax = 2, informationRates = c(0.5, 1),
          alpha = 0.025, typeAlphaSpending = "sfOF",
          lambdaH0 = 0.125, accrualIntensity = 500,
          stratumFraction = c(0.2, 0.8),
          kappa = c(3, 5), lambda = c(0.0875, 0.085),
          gamma = 0, accrualDuration = 1.25,
          followupTime = 2.75, fixedFollowup = FALSE)
#>                                                                             
#> Group-sequential design with 2 stages for one-sample negative binomial rate 
#> Rate under H0: 0.125, rate under H1: 0.0855                                 
#> Stratum fraction: 0.2 0.8, event rate: 0.0875 0.085, dispersion: 3 5        
#> Overall power: 0.9152, overall significance level (1-sided): 0.025          
#> Maximum # events: 180.4, expected # events: 146.3                           
#> Maximum # dropouts: 0, expected # dropouts: 0                               
#> Maximum # subjects: 625, expected # subjects: 625                           
#> Maximum exposure: 2109.4, expected exposure: 1711                           
#> Maximum information: 77.27, expected information: 66.69                     
#> Total study duration: 4, expected study duration: 3.4                       
#> Accrual duration: 1.2, follow-up duration: 2.8, fixed follow-up: FALSE      
#> Alpha spending: Lan-DeMets O'Brien-Fleming, beta spending: None             
#>                                                                             
#>                          Stage 1 Stage 2
#> Information rate         0.500   1.000  
#> Efficacy boundary (Z)    2.963   1.969  
#> Cumulative rejection     0.2738  0.9152 
#> Cumulative alpha spent   0.0015  0.0250 
#> Number of events         55.9    180.4  
#> Number of dropouts       0.0     0.0    
#> Number of subjects       625.0   625.0  
#> Exposure                 654.2   2109.4 
#> Analysis time            1.7     4.0    
#> Efficacy boundary (rate) 0.0776  0.0999 
#> Efficacy boundary (p)    0.0015  0.0245 
#> Information              38.63   77.27  

# Example 2: Fixed follow-up design

nbpower1s(kMax = 2, informationRates = c(0.5, 1),
          alpha = 0.025, typeAlphaSpending = "sfOF",
          lambdaH0 = 8.4, accrualIntensity = 40,
          kappa = 3, lambda = 0.5*8.4,
          gamma = 0, accrualDuration = 1.5,
          followupTime = 0.5, fixedFollowup = TRUE)
#>                                                                             
#> Group-sequential design with 2 stages for one-sample negative binomial rate 
#> Rate under H0: 8.4, rate under H1: 4.2                                      
#> Dispersion: 3                                                               
#> Overall power: 0.8198, overall significance level (1-sided): 0.025          
#> Maximum # events: 126, expected # events: 112.7                             
#> Maximum # dropouts: 0, expected # dropouts: 0                               
#> Maximum # subjects: 60, expected # subjects: 55.4                           
#> Maximum exposure: 30, expected exposure: 26.8                               
#> Maximum information: 17.26, expected information: 15.73                     
#> Total study duration: 2, expected study duration: 1.8                       
#> Accrual duration: 1.5, follow-up duration: 0.5, fixed follow-up: TRUE       
#> Alpha spending: Lan-DeMets O'Brien-Fleming, beta spending: None             
#>                                                                             
#>                          Stage 1 Stage 2
#> Information rate         0.500   1.000  
#> Efficacy boundary (Z)    2.963   1.969  
#> Cumulative rejection     0.1771  0.8198 
#> Cumulative alpha spent   0.0015  0.0250 
#> Number of events         50.7    126.0  
#> Number of dropouts       0.0     0.0    
#> Number of subjects       34.1    60.0   
#> Exposure                 12.1    30.0   
#> Analysis time            0.9     2.0    
#> Efficacy boundary (rate) 3.0641  5.2299 
#> Efficacy boundary (p)    0.0015  0.0245 
#> Information              8.63    17.26  
```
