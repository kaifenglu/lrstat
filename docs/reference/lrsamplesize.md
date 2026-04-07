# Log-Rank Test Sample Size

Obtains the needed accrual duration given power and follow-up time, the
needed follow-up time given power and accrual duration, or the needed
absolute accrual rates given power, accrual duration, follow-up time,
and relative accrual rates in a two-group survival design.

## Usage

``` r
lrsamplesize(
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
  hazardRatioH0 = 1,
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
  rho1 = 0,
  rho2 = 0,
  estimateHazardRatio = TRUE,
  typeOfComputation = "",
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

  The information rates in terms of number of events for the
  conventional log-rank test and in terms of the actual information for
  weighted log-rank tests. Defaults to `(1:kMax) / kMax` if left
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

- hazardRatioH0:

  Hazard ratio under the null hypothesis for the active treatment versus
  control. Defaults to 1 for superiority test.

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

- rho1:

  The first parameter of the Fleming-Harrington family of weighted
  log-rank test. Defaults to 0 for conventional log-rank test.

- rho2:

  The second parameter of the Fleming-Harrington family of weighted
  log-rank test. Defaults to 0 for conventional log-rank test.

- estimateHazardRatio:

  Whether to estimate the hazard ratio from weighted Cox regression
  model and report the stopping boundaries on the hazard ratio scale.

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

- rounding:

  Whether to round up sample size and events. Defaults to 1 for sample
  size rounding.

## Value

A list of two components:

- `resultsUnderH1`: An S3 class `lrpower` object under the alternative
  hypothesis.

- `resultsUnderH0`: An S3 class `lrpower` object under the null
  hypothesis.

## See also

[`lrpower`](https://kaifenglu.github.io/lrstat/reference/lrpower.md)

## Author

Kaifeng Lu, <kaifenglu@gmail.com>

## Examples

``` r
# Piecewise accrual, piecewise exponential survival, and 5% dropout by
# the end of 1 year.

# Example 1: Obtains accrual duration given power and follow-up time

lrsamplesize(beta = 0.2, kMax = 2,
             informationRates = c(0.8, 1),
             alpha = 0.025, typeAlphaSpending = "sfOF",
             accrualTime = seq(0, 8),
             accrualIntensity = 26/9*seq(1, 9),
             piecewiseSurvivalTime = c(0, 6),
             lambda1 = c(0.0533, 0.0309),
             lambda2 = c(0.0533, 0.0533),
             gamma1 = -log(1-0.05)/12,
             gamma2 = -log(1-0.05)/12,
             accrualDuration = NA,
             followupTime = 18, fixedFollowup = FALSE)
#> $resultsUnderH1
#>                                                                          
#> Group-sequential design with 2 stages for log-rank test                  
#> Overall power: 0.8011, overall significance level (1-sided): 0.025       
#> Maximum # events: 322, expected # events: 293.3                          
#> Maximum # dropouts: 30.4, expected # dropouts: 27.4                      
#> Maximum # subjects: 487, expected # subjects: 487                        
#> Maximum information: 79.75, expected information: 72.8                   
#> Total study duration: 40.8, expected study duration: 36.5                
#> Accrual duration: 22.7, follow-up duration: 18.1, fixed follow-up: FALSE 
#> Allocation ratio: 1                                                      
#> Alpha spending: Lan-DeMets O'Brien-Fleming, beta spending: None          
#>                                                                          
#>                        Stage 1 Stage 2
#> Information rate       0.801   1.000  
#> Efficacy boundary (Z)  2.248   2.025  
#> Cumulative rejection   0.4486  0.8011 
#> Cumulative alpha spent 0.0123  0.0250 
#> Number of events       258.0   322.0  
#> Number of dropouts     23.6    30.4   
#> Number of subjects     487.0   487.0  
#> Analysis time          31.3    40.8   
#> Efficacy boundary (HR) 0.755   0.797  
#> Efficacy boundary (p)  0.0123  0.0214 
#> Information            64.26   79.75  
#> HR                     0.767   0.726  
#> 
#> $resultsUnderH0
#>                                                                          
#> Group-sequential design with 2 stages for log-rank test                  
#> Overall power: 0.025, overall significance level (1-sided): 0.025        
#> Maximum # events: 322, expected # events: 321.2                          
#> Maximum # dropouts: 25.8, expected # dropouts: 25.8                      
#> Maximum # subjects: 487, expected # subjects: 487                        
#> Maximum information: 80.5, expected information: 80.3                    
#> Total study duration: 35.9, expected study duration: 35.8                
#> Accrual duration: 22.7, follow-up duration: 13.1, fixed follow-up: FALSE 
#> Allocation ratio: 1                                                      
#> Alpha spending: Lan-DeMets O'Brien-Fleming, beta spending: None          
#>                                                                          
#>                        Stage 1 Stage 2
#> Information rate       0.801   1.000  
#> Efficacy boundary (Z)  2.248   2.025  
#> Cumulative rejection   0.0123  0.0250 
#> Cumulative alpha spent 0.0123  0.0250 
#> Number of events       258.0   322.0  
#> Number of dropouts     20.7    25.8   
#> Number of subjects     487.0   487.0  
#> Analysis time          28.8    35.9   
#> Efficacy boundary (HR) 0.756   0.798  
#> Efficacy boundary (p)  0.0123  0.0214 
#> Information            64.50   80.50  
#> HR                     1.000   1.000  
#> 


# Example 2: Obtains follow-up time given power and accrual duration

lrsamplesize(beta = 0.2, kMax = 2,
             informationRates = c(0.8, 1),
             alpha = 0.025, typeAlphaSpending = "sfOF",
             accrualTime = seq(0, 8),
             accrualIntensity = 26/9*seq(1, 9),
             piecewiseSurvivalTime = c(0, 6),
             lambda1 = c(0.0533, 0.0309),
             lambda2 = c(0.0533, 0.0533),
             gamma1 = -log(1-0.05)/12,
             gamma2 = -log(1-0.05)/12,
             accrualDuration = 22,
             followupTime = NA, fixedFollowup = FALSE)
#> $resultsUnderH1
#>                                                                        
#> Group-sequential design with 2 stages for log-rank test                
#> Overall power: 0.8024, overall significance level (1-sided): 0.025     
#> Maximum # events: 315, expected # events: 286.5                        
#> Maximum # dropouts: 29.8, expected # dropouts: 26.8                    
#> Maximum # subjects: 468, expected # subjects: 468                      
#> Maximum information: 77.96, expected information: 71.09                
#> Total study duration: 41.5, expected study duration: 36.9              
#> Accrual duration: 22, follow-up duration: 19.5, fixed follow-up: FALSE 
#> Allocation ratio: 1                                                    
#> Alpha spending: Lan-DeMets O'Brien-Fleming, beta spending: None        
#>                                                                        
#>                        Stage 1 Stage 2
#> Information rate       0.800   1.000  
#> Efficacy boundary (Z)  2.250   2.025  
#> Cumulative rejection   0.4520  0.8024 
#> Cumulative alpha spent 0.0122  0.0250 
#> Number of events       252.0   315.0  
#> Number of dropouts     23.1    29.8   
#> Number of subjects     468.0   468.0  
#> Analysis time          31.4    41.5   
#> Efficacy boundary (HR) 0.752   0.795  
#> Efficacy boundary (p)  0.0122  0.0214 
#> Information            62.76   77.96  
#> HR                     0.764   0.723  
#> 
#> $resultsUnderH0
#>                                                                        
#> Group-sequential design with 2 stages for log-rank test                
#> Overall power: 0.025, overall significance level (1-sided): 0.025      
#> Maximum # events: 315, expected # events: 314.2                        
#> Maximum # dropouts: 25.3, expected # dropouts: 25.2                    
#> Maximum # subjects: 468, expected # subjects: 468                      
#> Maximum information: 78.75, expected information: 78.56                
#> Total study duration: 36.2, expected study duration: 36.1              
#> Accrual duration: 22, follow-up duration: 14.2, fixed follow-up: FALSE 
#> Allocation ratio: 1                                                    
#> Alpha spending: Lan-DeMets O'Brien-Fleming, beta spending: None        
#>                                                                        
#>                        Stage 1 Stage 2
#> Information rate       0.800   1.000  
#> Efficacy boundary (Z)  2.250   2.025  
#> Cumulative rejection   0.0122  0.0250 
#> Cumulative alpha spent 0.0122  0.0250 
#> Number of events       252.0   315.0  
#> Number of dropouts     20.2    25.3   
#> Number of subjects     468.0   468.0  
#> Analysis time          28.8    36.2   
#> Efficacy boundary (HR) 0.753   0.796  
#> Efficacy boundary (p)  0.0122  0.0214 
#> Information            63.00   78.75  
#> HR                     1.000   1.000  
#> 


# Example 3: Obtains absolute accrual intensity given power,
# accrual duration, follow-up time, and relative accrual intensity

lrsamplesize(beta = 0.2, kMax = 2,
             informationRates = c(0.8, 1),
             alpha = 0.025, typeAlphaSpending = "sfOF",
             accrualTime = seq(0, 8),
             accrualIntensity = 26/9*seq(1, 9),
             piecewiseSurvivalTime = c(0, 6),
             lambda1 = c(0.0533, 0.0309),
             lambda2 = c(0.0533, 0.0533),
             gamma1 = -log(1-0.05)/12,
             gamma2 = -log(1-0.05)/12,
             accrualDuration = 22,
             followupTime = 18, fixedFollowup = FALSE)
#> $resultsUnderH1
#>                                                                      
#> Group-sequential design with 2 stages for log-rank test              
#> Overall power: 0.8003, overall significance level (1-sided): 0.025   
#> Maximum # events: 324, expected # events: 295.3                      
#> Maximum # dropouts: 30.6, expected # dropouts: 27.5                  
#> Maximum # subjects: 493, expected # subjects: 493                    
#> Maximum information: 80.27, expected information: 73.32              
#> Total study duration: 40, expected study duration: 35.9              
#> Accrual duration: 22, follow-up duration: 18, fixed follow-up: FALSE 
#> Allocation ratio: 1                                                  
#> Alpha spending: Lan-DeMets O'Brien-Fleming, beta spending: None      
#>                                                                      
#>                        Stage 1 Stage 2
#> Information rate       0.799   1.000  
#> Efficacy boundary (Z)  2.251   2.025  
#> Cumulative rejection   0.4414  0.8003 
#> Cumulative alpha spent 0.0122  0.0250 
#> Number of events       259.0   324.0  
#> Number of dropouts     23.7    30.6   
#> Number of subjects     493.0   493.0  
#> Analysis time          30.6    40.0   
#> Efficacy boundary (HR) 0.755   0.797  
#> Efficacy boundary (p)  0.0122  0.0214 
#> Information            64.52   80.27  
#> HR                     0.769   0.727  
#> 
#> $resultsUnderH0
#>                                                                        
#> Group-sequential design with 2 stages for log-rank test                
#> Overall power: 0.025, overall significance level (1-sided): 0.025      
#> Maximum # events: 324, expected # events: 323.2                        
#> Maximum # dropouts: 26, expected # dropouts: 25.9                      
#> Maximum # subjects: 493, expected # subjects: 493                      
#> Maximum information: 81, expected information: 80.8                    
#> Total study duration: 35.2, expected study duration: 35.1              
#> Accrual duration: 22, follow-up duration: 13.2, fixed follow-up: FALSE 
#> Allocation ratio: 1                                                    
#> Alpha spending: Lan-DeMets O'Brien-Fleming, beta spending: None        
#>                                                                        
#>                        Stage 1 Stage 2
#> Information rate       0.799   1.000  
#> Efficacy boundary (Z)  2.251   2.025  
#> Cumulative rejection   0.0122  0.0250 
#> Cumulative alpha spent 0.0122  0.0250 
#> Number of events       259.0   324.0  
#> Number of dropouts     20.8    26.0   
#> Number of subjects     493.0   493.0  
#> Analysis time          28.2    35.2   
#> Efficacy boundary (HR) 0.756   0.799  
#> Efficacy boundary (p)  0.0122  0.0214 
#> Information            64.75   81.00  
#> HR                     1.000   1.000  
#> 
```
