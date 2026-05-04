# Schoenfeld Method for Log-Rank Test Sample Size Calculation

Obtains the sample size and study duration by calibrating the number of
events calculated using the Schoenfeld formula under the proportional
hazards assumption.

## Usage

``` r
lrschoenfeld(
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
  futilityCP = NA_real_,
  futilityHR = NA_real_,
  typeBetaSpending = "none",
  parameterBetaSpending = NA_real_,
  userBetaSpending = NA_real_,
  hazardRatioH0 = 1,
  allocationRatioPlanned = 1,
  accrualTime = 0L,
  accrualIntensity = NA_real_,
  piecewiseSurvivalTime = 0L,
  stratumFraction = 1L,
  hazardRatio = NA_real_,
  lambda2 = NA_real_,
  gamma1 = 0L,
  gamma2 = 0L,
  followupTime = NA_real_,
  fixedFollowup = 0L,
  spendingTime = NA_real_,
  rounding = 1L,
  calibrate = 1L,
  maxNumberOfIterations = 10000L,
  maxNumberOfRawDatasetsPerStage = 0L,
  seed = NA_integer_
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

- futilityHR:

  A vector of length `kMax - 1` for the futility bounds on the hazard
  ratio scale.

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

- hazardRatio:

  Hazard ratio under the alternative hypothesis for the active treatment
  versus control.

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

  Whether to round up sample size and events. Defaults to 1 for sample
  size rounding.

- calibrate:

  Whether to use simulations to calibrate the number of events
  calculated using the Schoenfeld formula.

- maxNumberOfIterations:

  The number of simulation iterations. Defaults to 10000.

- maxNumberOfRawDatasetsPerStage:

  The number of raw datasets per stage to extract.

- seed:

  The seed to reproduce the simulation results. The seed from the
  environment will be used if left unspecified.

## Value

A list of two components:

- `analyticalResults`: An S3 class `lrpower` object for the asymptotic
  power.

- `simulationResults`: An S3 class `lrsim` object for the empirical
  power.

## Details

This function calculates the sample size and study duration by
calibrating the number of events estimated using the Schoenfeld formula
under the proportional hazards assumption, particularly when the hazard
ratio is far away from one and/or the allocation between groups is
unequal.

For a fixed design, the Schoenfeld formula for the required number of
events is \$\$D = \frac{(\Phi^{-1}(1-\alpha) + \Phi^{-1}(1-\beta))^2}
{(\theta - \theta_0)^2 r(1-r)}\$\$ where \\D\\ is the total number of
events required, \\\alpha\\ is the type I error rate, \\\beta\\ is the
type II error rate, \\r\\ is the randomization probability for the
active treatment group, \\\theta_0\\ and \\\theta\\ are the log hazard
ratios under the null and alternative hypotheses, respectively.

The function first computes the number of events using the Schoenfeld
formula. If `calibrate` is set to 1, the function uses simulations to
calibrate the number of events, accounting for scenarios where the
Schoenfeld formula may be inaccurate (e.g., when allocation is unequal
or the hazard ratio is extreme).

Let \\D\_{schoenfeld}\\ be the number of events calculated by the
Schoenfeld formula, and \\D\_{calibrated}\\ be the calibrated number of
events. The calibrated number of events is calculated as \#'
\$\$D\_{\text{calibrated}} = \frac{\left\\\Phi^{-1}(1-\alpha) +
\Phi^{-1}(1-\beta)\right\\^2} {\left\\\Phi^{-1}(1-\alpha) +
\Phi^{-1}(1-\beta\_{\text{schoenfeld}})\right\\^2}
D\_{\text{schoenfeld}}\$\$ where \\\beta\_{schoenfeld}\\ is the
empirical type II error estimated via simulation.

A second round of simulation is performed to obtain the empirical power
using the calibrated number of events.

## Author

Kaifeng Lu, <kaifenglu@gmail.com>

## Examples

``` r
(lr1 <- lrschoenfeld(
  beta = 0.1, kMax = 2, alpha = 0.025,
  hazardRatioH0 = 1, allocationRatioPlanned = 1,
  accrualIntensity = 20, hazardRatio = 0.3,
  lambda2 = 1.9/12,
  gamma1 = -log(1-0.1)/24, gamma2 = -log(1-0.1)/24,
  fixedFollowup = 0, rounding = 1,
  calibrate = 0, maxNumberOfIterations = 1000,
  seed = 12345))
#> $analyticalResults
#>                                                                        
#> Group-sequential design with 2 stages for log-rank test                
#> Overall power: 0.9085, overall significance level (1-sided): 0.025     
#> Maximum # events: 30, expected # events: 26                            
#> Maximum # dropouts: 1.4, expected # dropouts: 1.2                      
#> Maximum # subjects: 78, expected # subjects: 78                        
#> Maximum information: 7.5, expected information: 6.51                   
#> Total study duration: 7.2, expected study duration: 6.4                
#> Accrual duration: 3.9, follow-up duration: 3.3, fixed follow-up: FALSE 
#> Allocation ratio: 1                                                    
#> Alpha spending: Lan-DeMets O'Brien-Fleming, beta spending: None        
#>                                                                        
#>                        Stage 1 Stage 2
#> Information rate       0.500   1.000  
#> Efficacy boundary (Z)  2.963   1.969  
#> Cumulative rejection   0.2640  0.9085 
#> Cumulative alpha spent 0.0015  0.0250 
#> Number of events       15.0    30.0   
#> Number of dropouts     0.7     1.4    
#> Number of subjects     78.0    78.0   
#> Analysis time          4.2     7.2    
#> Efficacy boundary (HR) 0.217   0.487  
#> Efficacy boundary (p)  0.0015  0.0245 
#> Information            3.75    7.50   
#> HR                     0.300   0.300  
#> 
#> $simulationResults
#>                                                         
#> Group-sequential design with 2 stages for log-rank test 
#> Empirical power: 0.878                                  
#> Expected # events: 26.8                                 
#> Expected # dropouts: 1.2                                
#> Expected # subjects: 77.5                               
#> Expected study duration: 6.5                            
#> n: 78, fixed follow-up: FALSE                           
#> Number of simulations: 1000                             
#>                                                         
#>                      Stage 1 Stage 2
#> Cumulative rejection 0.2140  0.8780 
#> Cumulative futility  0.0000  0.1220 
#> Number of events     15.0    30.0   
#> Number of dropouts   0.7     1.4    
#> Number of subjects   75.6    78.0   
#> Analysis time        4.2     7.1    
#> 

(lr2 <- lrschoenfeld(
  beta = 0.1, kMax = 2, alpha = 0.025,
  hazardRatioH0 = 1, allocationRatioPlanned = 1,
  accrualIntensity = 20, hazardRatio = 0.3,
  lambda2 = 1.9/12,
  gamma1 = -log(1-0.1)/24, gamma2 = -log(1-0.1)/24,
  fixedFollowup = 0, rounding = 1,
  calibrate = 1, maxNumberOfIterations = 1000,
  seed = 12345))
#> $analyticalResults
#>                                                                        
#> Group-sequential design with 2 stages for log-rank test                
#> Overall power: 0.9321, overall significance level (1-sided): 0.025     
#> Maximum # events: 33, expected # events: 27.7                          
#> Maximum # dropouts: 1.5, expected # dropouts: 1.3                      
#> Maximum # subjects: 78, expected # subjects: 78                        
#> Maximum information: 8.25, expected information: 6.92                  
#> Total study duration: 8, expected study duration: 6.8                  
#> Accrual duration: 3.9, follow-up duration: 4.1, fixed follow-up: FALSE 
#> Allocation ratio: 1                                                    
#> Alpha spending: Lan-DeMets O'Brien-Fleming, beta spending: None        
#>                                                                        
#>                        Stage 1 Stage 2
#> Information rate       0.515   1.000  
#> Efficacy boundary (Z)  2.913   1.970  
#> Cumulative rejection   0.3333  0.9321 
#> Cumulative alpha spent 0.0018  0.0250 
#> Number of events       17.0    33.0   
#> Number of dropouts     0.8     1.5    
#> Number of subjects     78.0    78.0   
#> Analysis time          4.5     8.0    
#> Efficacy boundary (HR) 0.243   0.504  
#> Efficacy boundary (p)  0.0018  0.0244 
#> Information            4.25    8.25   
#> HR                     0.300   0.300  
#> 
#> $simulationResults
#>                                                         
#> Group-sequential design with 2 stages for log-rank test 
#> Empirical power: 0.915                                  
#> Expected # events: 28.3                                 
#> Expected # dropouts: 1.3                                
#> Expected # subjects: 77.7                               
#> Expected study duration: 6.9                            
#> n: 78, fixed follow-up: FALSE                           
#> Number of simulations: 1000                             
#>                                                         
#>                      Stage 1 Stage 2
#> Cumulative rejection 0.2950  0.9150 
#> Cumulative futility  0.0000  0.0850 
#> Number of events     17.0    33.0   
#> Number of dropouts   0.8     1.5    
#> Number of subjects   77.1    78.0   
#> Analysis time        4.5     7.9    
#> 
```
