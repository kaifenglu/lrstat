# Sample Size for Difference in Milestone Survival Probabilities

Obtains the needed accrual duration given power, accrual intensity, and
follow-up time, the needed follow-up time given power, accrual
intensity, and accrual duration, or the needed absolute accrual
intensity given power, relative accrual intensity, accrual duration, and
follow-up time in a two-group survival design.

## Usage

``` r
kmsamplesize(
  beta = 0.2,
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
  userBetaSpending = NA_real_,
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
  rounding = TRUE
)
```

## Arguments

- beta:

  Type II error. Defaults to 0.2.

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

  The type of beta spending. One of the following: `"sfOF"` for
  O'Brien-Fleming type spending function, `"sfP"` for Pocock type
  spending function, `"sfKD"` for Kim & DeMets spending function,
  `"sfHSD"` for Hwang, Shi & DeCani spending function, `"user"` for user
  defined spending, and `"none"` for no early futility stopping.
  Defaults to `"none"`.

- parameterBetaSpending:

  The parameter value for the beta spending. Corresponds to \\\rho\\ for
  `"sfKD"`, and \\\gamma\\ for `"sfHSD"`.

- userBetaSpending:

  The user defined beta spending. Cumulative beta spent up to each
  stage.

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

- rounding:

  Whether to round up sample size. Defaults to 1 for sample size
  rounding.

## Value

A list of two components:

- `resultsUnderH1`: An S3 class `kmpower` object under the alternative
  hypothesis.

- `resultsUnderH0`: An S3 class `kmpower` object under the null
  hypothesis.

## See also

[`kmpower`](https://kaifenglu.github.io/lrstat/reference/kmpower.md)

## Author

Kaifeng Lu, <kaifenglu@gmail.com>

## Examples

``` r
# Example 1: Obtains follow-up time given power, accrual intensity,
# and accrual duration for variable follow-up. Of note, the power
# reaches the maximum when the follow-up time equals milestone.

kmsamplesize(beta = 0.25, kMax = 2, informationRates = c(0.8, 1),
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
             followupTime = NA, fixedFollowup = FALSE)
#> $resultsUnderH1
#>                                                                            
#> Group-sequential design with 2 stages for difference in milestone survival 
#> Milestone: 18, survival difference under H0: 0                             
#> Milestone survival on treatment: 0.384, on control: 0.266                  
#> Overall power: 0.75, overall significance level (1-sided): 0.025           
#> Maximum # events: 338.2, expected # events: 315.4                          
#> Maximum # subjects: 468, expected # subjects: 468                          
#> Maximum information: 510.57, expected information: 454.18                  
#> Total study duration: 36.4, expected study duration: 33.2                  
#> Accrual duration: 22, follow-up duration: 14.4, fixed follow-up: FALSE     
#> Allocation ratio: 1                                                        
#> Alpha spending: Lan-DeMets O'Brien-Fleming, beta spending: None            
#>                                                                            
#>                               Stage 1 Stage 2
#> Information rate              0.800   1.000  
#> Efficacy boundary (Z)         2.250   2.025  
#> Cumulative rejection          0.5522  0.7500 
#> Cumulative alpha spent        0.0122  0.0250 
#> Number of events              297.0   338.2  
#> Number of dropouts            19.4    22.7   
#> Number of subjects            468.0   468.0  
#> Number of milestone subjects  66.8    113.1  
#> Analysis time                 30.5    36.4   
#> Efficacy boundary (surv diff) 0.111   0.090  
#> Efficacy boundary (p)         0.0122  0.0214 
#> Information                   408.46  510.57 
#> 
#> $resultsUnderH0
#>                                                                            
#> Group-sequential design with 2 stages for difference in milestone survival 
#> Milestone: 18, survival difference under H0: 0                             
#> Milestone survival on treatment: 0.266, on control: 0.266                  
#> Overall power: 0.025, overall significance level (1-sided): 0.025          
#> Maximum # events: 345.4, expected # events: 345                            
#> Maximum # subjects: 468, expected # subjects: 468                          
#> Maximum information: 510.57, expected information: 509.33                  
#> Total study duration: 33.5, expected study duration: 33.5                  
#> Accrual duration: 22, follow-up duration: 11.5, fixed follow-up: FALSE     
#> Allocation ratio: 1                                                        
#> Alpha spending: Lan-DeMets O'Brien-Fleming, beta spending: None            
#>                                                                            
#>                               Stage 1 Stage 2
#> Information rate              0.800   1.000  
#> Efficacy boundary (Z)         2.250   2.025  
#> Cumulative rejection          0.0122  0.0250 
#> Cumulative alpha spent        0.0122  0.0250 
#> Number of events              311.9   345.4  
#> Number of dropouts            18.1    20.0   
#> Number of subjects            468.0   468.0  
#> Number of milestone subjects  48.9    73.8   
#> Analysis time                 29.6    33.5   
#> Efficacy boundary (surv diff) 0.111   0.090  
#> Efficacy boundary (p)         0.0122  0.0214 
#> Information                   408.46  510.57 
#> 

# Example 2: Obtains accrual intensity given power, accrual duration, and
# follow-up time for variable follow-up

kmsamplesize(beta = 0.2, kMax = 2, informationRates = c(0.8, 1),
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
#> $resultsUnderH1
#>                                                                            
#> Group-sequential design with 2 stages for difference in milestone survival 
#> Milestone: 18, survival difference under H0: 0                             
#> Milestone survival on treatment: 0.384, on control: 0.266                  
#> Overall power: 0.8, overall significance level (1-sided): 0.025            
#> Maximum # events: 386.9, expected # events: 352.3                          
#> Maximum # subjects: 513, expected # subjects: 513                          
#> Maximum information: 577.11, expected information: 506.6                   
#> Total study duration: 39.3, expected study duration: 34.3                  
#> Accrual duration: 22, follow-up duration: 17.3, fixed follow-up: FALSE     
#> Allocation ratio: 1                                                        
#> Alpha spending: Lan-DeMets O'Brien-Fleming, beta spending: None            
#>                                                                            
#>                               Stage 1 Stage 2
#> Information rate              0.800   1.000  
#> Efficacy boundary (Z)         2.250   2.025  
#> Cumulative rejection          0.6109  0.8000 
#> Cumulative alpha spent        0.0122  0.0250 
#> Number of events              330.3   386.9  
#> Number of dropouts            21.7    26.2   
#> Number of subjects            513.0   513.0  
#> Number of milestone subjects  77.8    148.2  
#> Analysis time                 31.1    39.3   
#> Efficacy boundary (surv diff) 0.105   0.084  
#> Efficacy boundary (p)         0.0122  0.0214 
#> Information                   461.69  577.11 
#> 
#> $resultsUnderH0
#>                                                                            
#> Group-sequential design with 2 stages for difference in milestone survival 
#> Milestone: 18, survival difference under H0: 0                             
#> Milestone survival on treatment: 0.266, on control: 0.266                  
#> Overall power: 0.025, overall significance level (1-sided): 0.025          
#> Maximum # events: 384.9, expected # events: 384.4                          
#> Maximum # subjects: 513, expected # subjects: 513                          
#> Maximum information: 577.11, expected information: 575.7                   
#> Total study duration: 34.3, expected study duration: 34.2                  
#> Accrual duration: 22, follow-up duration: 12.3, fixed follow-up: FALSE     
#> Allocation ratio: 1                                                        
#> Alpha spending: Lan-DeMets O'Brien-Fleming, beta spending: None            
#>                                                                            
#>                               Stage 1 Stage 2
#> Information rate              0.800   1.000  
#> Efficacy boundary (Z)         2.250   2.025  
#> Cumulative rejection          0.0122  0.0250 
#> Cumulative alpha spent        0.0122  0.0250 
#> Number of events              346.4   384.9  
#> Number of dropouts            20.1    22.3   
#> Number of subjects            513.0   513.0  
#> Number of milestone subjects  56.5    86.4   
#> Analysis time                 30.0    34.3   
#> Efficacy boundary (surv diff) 0.105   0.084  
#> Efficacy boundary (p)         0.0122  0.0214 
#> Information                   461.69  577.11 
#> 


# Example 3: Obtains accrual duration given power, accrual intensity, and
# follow-up time for fixed follow-up

kmsamplesize(beta = 0.2, kMax = 2, informationRates = c(0.8, 1),
             alpha = 0.025, typeAlphaSpending = "sfOF",
             milestone = 18,
             allocationRatioPlanned = 1, accrualTime = seq(0, 8),
             accrualIntensity = 26/9*seq(1, 9),
             piecewiseSurvivalTime = c(0, 6),
             stratumFraction = c(0.2, 0.8),
             lambda1 = c(0.0533, 0.0309, 1.5*0.0533, 1.5*0.0309),
             lambda2 = c(0.0533, 0.0533, 1.5*0.0533, 1.5*0.0533),
             gamma1 = -log(1-0.05)/12,
             gamma2 = -log(1-0.05)/12, accrualDuration = NA,
             followupTime = 18, fixedFollowup = TRUE)
#> $resultsUnderH1
#>                                                                            
#> Group-sequential design with 2 stages for difference in milestone survival 
#> Milestone: 18, survival difference under H0: 0                             
#> Milestone survival on treatment: 0.384, on control: 0.266                  
#> Overall power: 0.8, overall significance level (1-sided): 0.025            
#> Maximum # events: 336.3, expected # events: 321.7                          
#> Maximum # subjects: 513, expected # subjects: 513                          
#> Maximum information: 577.11, expected information: 506.6                   
#> Total study duration: 41, expected study duration: 35.6                    
#> Accrual duration: 23.7, follow-up duration: 18, fixed follow-up: TRUE      
#> Allocation ratio: 1                                                        
#> Alpha spending: Lan-DeMets O'Brien-Fleming, beta spending: None            
#>                                                                            
#>                               Stage 1 Stage 2
#> Information rate              0.800   1.000  
#> Efficacy boundary (Z)         2.250   2.025  
#> Cumulative rejection          0.6109  0.8000 
#> Cumulative alpha spent        0.0122  0.0250 
#> Number of events              312.4   336.3  
#> Number of dropouts            20.3    22.1   
#> Number of subjects            513.0   513.0  
#> Number of milestone subjects  80.2    148.5  
#> Analysis time                 32.2    41.0   
#> Efficacy boundary (surv diff) 0.105   0.084  
#> Efficacy boundary (p)         0.0122  0.0214 
#> Information                   461.69  577.11 
#> 
#> $resultsUnderH0
#>                                                                            
#> Group-sequential design with 2 stages for difference in milestone survival 
#> Milestone: 18, survival difference under H0: 0                             
#> Milestone survival on treatment: 0.266, on control: 0.266                  
#> Overall power: 0.025, overall significance level (1-sided): 0.025          
#> Maximum # events: 330.8, expected # events: 330.5                          
#> Maximum # subjects: 464.5, expected # subjects: 464.5                      
#> Maximum information: 577.11, expected information: 575.7                   
#> Total study duration: 39.9, expected study duration: 39.8                  
#> Accrual duration: 21.9, follow-up duration: 18, fixed follow-up: TRUE      
#> Allocation ratio: 1                                                        
#> Alpha spending: Lan-DeMets O'Brien-Fleming, beta spending: None            
#>                                                                            
#>                               Stage 1 Stage 2
#> Information rate              0.800   1.000  
#> Efficacy boundary (Z)         2.250   2.025  
#> Cumulative rejection          0.0122  0.0250 
#> Cumulative alpha spent        0.0122  0.0250 
#> Number of events              310.4   330.8  
#> Number of dropouts            18.0    19.2   
#> Number of subjects            464.5   464.5  
#> Number of milestone subjects  61.0    114.6  
#> Analysis time                 31.5    39.9   
#> Efficacy boundary (surv diff) 0.105   0.084  
#> Efficacy boundary (p)         0.0122  0.0214 
#> Information                   461.69  577.11 
#> 
```
