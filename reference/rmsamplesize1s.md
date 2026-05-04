# Sample Size for One-Sample Restricted Mean Survival Time

Obtains the needed accrual duration given power and follow-up time, the
needed follow-up time given power and accrual duration, or the needed
absolute accrual rates given power, accrual duration, follow-up
duration, and relative accrual rates in a one-group survival design.

## Usage

``` r
rmsamplesize1s(
  beta = 0.2,
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
  futilityRmst = NULL,
  typeBetaSpending = "none",
  parameterBetaSpending = NA_real_,
  userBetaSpending = NA_real_,
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
  for `"WT"`, \\\rho\\ for `"sfKD"`, and \\\gamma\\ for `"sfHSD"`.

- userAlphaSpending:

  The user defined alpha spending. Cumulative alpha spent up to each
  stage.

- futilityBounds:

  Lower boundaries on the z-test statistic scale for stopping for
  futility at stages `1, ..., kMax-1`. Defaults to `rep(-8, kMax-1)` if
  left unspecified. The futility bounds are non-binding for the
  calculation of critical values.

- futilityCP:

  A vector of length `kMax - 1` for the futility bounds on the
  conditional power scale.

- futilityRmst:

  A vector of length `kMax - 1` for the futility bounds on the
  restricted mean survival time scale.

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

  The milestone time at which to calculate the restricted survival time.

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

- rounding:

  Whether to round up sample size. Defaults to 1 for sample size
  rounding.

## Value

A list of two components:

- `resultsUnderH1`: An S3 class `rmpower1s` object under the alternative
  hypothesis.

- `resultsUnderH0`: An S3 class `rmpower1s` object under the null
  hypothesis.

## See also

[`rmpower1s`](https://kaifenglu.github.io/lrstat/reference/rmpower1s.md)

## Author

Kaifeng Lu, <kaifenglu@gmail.com>

## Examples

``` r
# Example 1: Obtains follow-up duration given power, accrual intensity,
# and accrual duration for variable follow-up

rmsamplesize1s(beta = 0.2, kMax = 2,
               informationRates = c(0.8, 1),
               alpha = 0.025, typeAlphaSpending = "sfOF",
               milestone = 18, rmstH0 = 10,
               accrualTime = seq(0, 8),
               accrualIntensity = 26/9*seq(1, 9),
               piecewiseSurvivalTime = c(0, 6),
               stratumFraction = c(0.2, 0.8),
               lambda = c(0.0533, 0.0309, 1.5*0.0533, 1.5*0.0309),
               gamma = -log(1-0.05)/12, accrualDuration = 22,
               followupTime = NA, fixedFollowup = FALSE)
#> $resultsUnderH1
#>                                                                                    
#> Group-sequential design with 2 stages for one-sample restricted mean survival time 
#> Milestone: 18, restricted mean survival time under H0: 10, under H1: 10.854        
#> Overall power: 0.8, overall significance level (1-sided): 0.025                    
#> Maximum # events: 356.8, expected # events: 279.4                                  
#> Maximum # subjects: 522, expected # subjects: 522                                  
#> Maximum # milestone subjects: 159.4, expected # milestone subjects: 79.2           
#> Maximum information: 10.99, expected information: 9.65                             
#> Total study duration: 39.2, expected study duration: 30.5                          
#> Accrual duration: 24.1, follow-up duration: 15.2, fixed follow-up: FALSE           
#> Alpha spending: Lan-DeMets O'Brien-Fleming, beta spending: None                    
#>                                                                                    
#>                              Stage 1 Stage 2
#> Information rate             0.800   1.000  
#> Efficacy boundary (Z)        2.250   2.025  
#> Cumulative rejection         0.6109  0.8000 
#> Cumulative alpha spent       0.0122  0.0250 
#> Number of events             230.1   356.8  
#> Number of dropouts           15.9    27.9   
#> Number of subjects           522.0   522.0  
#> Number of milestone subjects 28.2    159.4  
#> Analysis time                24.9    39.2   
#> Efficacy boundary (rmst)     10.759  10.611 
#> Efficacy boundary (p)        0.0122  0.0214 
#> Information                  8.79    10.99  
#> 
#> $resultsUnderH0
#>                                                                                    
#> Group-sequential design with 2 stages for one-sample restricted mean survival time 
#> Milestone: 18, restricted mean survival time under H0: 10, under H1: 10            
#> Overall power: 0.025, overall significance level (1-sided): 0.025                  
#> Maximum # events: 360.4, expected # events: 359.1                                  
#> Maximum # subjects: 522, expected # subjects: 522                                  
#> Maximum # milestone subjects: 102.9, expected # milestone subjects: 102            
#> Maximum information: 10.99, expected information: 10.97                            
#> Total study duration: 35.2, expected study duration: 35                            
#> Accrual duration: 24.1, follow-up duration: 11.1, fixed follow-up: FALSE           
#> Alpha spending: Lan-DeMets O'Brien-Fleming, beta spending: None                    
#>                                                                                    
#>                              Stage 1 Stage 2
#> Information rate             0.800   1.000  
#> Efficacy boundary (Z)        2.250   2.025  
#> Cumulative rejection         0.0122  0.0250 
#> Cumulative alpha spent       0.0122  0.0250 
#> Number of events             255.2   360.4  
#> Number of dropouts           14.9    23.1   
#> Number of subjects           522.0   522.0  
#> Number of milestone subjects 23.4    102.9  
#> Analysis time                24.9    35.2   
#> Efficacy boundary (rmst)     10.759  10.611 
#> Efficacy boundary (p)        0.0122  0.0214 
#> Information                  8.79    10.99  
#> 

# Example 2: Obtains accrual intensity given power, accrual duration, and
# follow-up duration for variable follow-up

rmsamplesize1s(beta = 0.2, kMax = 2,
               informationRates = c(0.8, 1),
               alpha = 0.025, typeAlphaSpending = "sfOF",
               milestone = 18, rmstH0 = 10,
               accrualTime = seq(0, 8),
               accrualIntensity = 26/9*seq(1, 9),
               piecewiseSurvivalTime = c(0, 6),
               stratumFraction = c(0.2, 0.8),
               lambda = c(0.0533, 0.0309, 1.5*0.0533, 1.5*0.0309),
               gamma = -log(1-0.05)/12, accrualDuration = 22,
               followupTime = 18, fixedFollowup = FALSE)
#> $resultsUnderH1
#>                                                                                    
#> Group-sequential design with 2 stages for one-sample restricted mean survival time 
#> Milestone: 18, restricted mean survival time under H0: 10, under H1: 10.854        
#> Overall power: 0.8, overall significance level (1-sided): 0.025                    
#> Maximum # events: 352.2, expected # events: 276.8                                  
#> Maximum # subjects: 522, expected # subjects: 522                                  
#> Maximum # milestone subjects: 157.2, expected # milestone subjects: 73.7           
#> Maximum information: 10.99, expected information: 9.65                             
#> Total study duration: 37.2, expected study duration: 28.8                          
#> Accrual duration: 22, follow-up duration: 15.2, fixed follow-up: FALSE             
#> Alpha spending: Lan-DeMets O'Brien-Fleming, beta spending: None                    
#>                                                                                    
#>                              Stage 1 Stage 2
#> Information rate             0.800   1.000  
#> Efficacy boundary (Z)        2.250   2.025  
#> Cumulative rejection         0.6109  0.8000 
#> Cumulative alpha spent       0.0122  0.0250 
#> Number of events             228.8   352.2  
#> Number of dropouts           15.7    27.4   
#> Number of subjects           522.0   522.0  
#> Number of milestone subjects 20.5    157.2  
#> Analysis time                23.5    37.2   
#> Efficacy boundary (rmst)     10.759  10.611 
#> Efficacy boundary (p)        0.0122  0.0214 
#> Information                  8.79    10.99  
#> 
#> $resultsUnderH0
#>                                                                                    
#> Group-sequential design with 2 stages for one-sample restricted mean survival time 
#> Milestone: 18, restricted mean survival time under H0: 10, under H1: 10            
#> Overall power: 0.025, overall significance level (1-sided): 0.025                  
#> Maximum # events: 356, expected # events: 354.8                                    
#> Maximum # subjects: 522, expected # subjects: 522                                  
#> Maximum # milestone subjects: 98.6, expected # milestone subjects: 97.6            
#> Maximum information: 10.99, expected information: 10.97                            
#> Total study duration: 33.3, expected study duration: 33.2                          
#> Accrual duration: 22, follow-up duration: 11.3, fixed follow-up: FALSE             
#> Alpha spending: Lan-DeMets O'Brien-Fleming, beta spending: None                    
#>                                                                                    
#>                              Stage 1 Stage 2
#> Information rate             0.800   1.000  
#> Efficacy boundary (Z)        2.250   2.025  
#> Cumulative rejection         0.0122  0.0250 
#> Cumulative alpha spent       0.0122  0.0250 
#> Number of events             254.4   356.0  
#> Number of dropouts           14.7    22.7   
#> Number of subjects           522.0   522.0  
#> Number of milestone subjects 17.0    98.6   
#> Analysis time                23.4    33.3   
#> Efficacy boundary (rmst)     10.759  10.611 
#> Efficacy boundary (p)        0.0122  0.0214 
#> Information                  8.79    10.99  
#> 


# Example 3: Obtains accrual duration given power, accrual intensity, and
# follow-up duration for fixed follow-up

rmsamplesize1s(beta = 0.2, kMax = 2,
               informationRates = c(0.8, 1),
               alpha = 0.025, typeAlphaSpending = "sfOF",
               milestone = 18, rmstH0 = 10,
               accrualTime = seq(0, 8),
               accrualIntensity = 26/9*seq(1, 9),
               piecewiseSurvivalTime = c(0, 6),
               stratumFraction = c(0.2, 0.8),
               lambda = c(0.0533, 0.0309, 1.5*0.0533, 1.5*0.0309),
               gamma = -log(1-0.05)/12, accrualDuration = NA,
               followupTime = 18, fixedFollowup = TRUE)
#> $resultsUnderH1
#>                                                                                    
#> Group-sequential design with 2 stages for one-sample restricted mean survival time 
#> Milestone: 18, restricted mean survival time under H0: 10, under H1: 10.854        
#> Overall power: 0.8, overall significance level (1-sided): 0.025                    
#> Maximum # events: 311.2, expected # events: 260                                    
#> Maximum # subjects: 522, expected # subjects: 522                                  
#> Maximum # milestone subjects: 159.4, expected # milestone subjects: 79.2           
#> Maximum information: 10.99, expected information: 9.65                             
#> Total study duration: 39.2, expected study duration: 30.5                          
#> Accrual duration: 24.1, follow-up duration: 18, fixed follow-up: TRUE              
#> Alpha spending: Lan-DeMets O'Brien-Fleming, beta spending: None                    
#>                                                                                    
#>                              Stage 1 Stage 2
#> Information rate             0.800   1.000  
#> Efficacy boundary (Z)        2.250   2.025  
#> Cumulative rejection         0.6109  0.8000 
#> Cumulative alpha spent       0.0122  0.0250 
#> Number of events             227.4   311.2  
#> Number of dropouts           15.6    23.3   
#> Number of subjects           522.0   522.0  
#> Number of milestone subjects 28.2    159.4  
#> Analysis time                24.9    39.2   
#> Efficacy boundary (rmst)     10.759  10.611 
#> Efficacy boundary (p)        0.0122  0.0214 
#> Information                  8.79    10.99  
#> 
#> $resultsUnderH0
#>                                                                                    
#> Group-sequential design with 2 stages for one-sample restricted mean survival time 
#> Milestone: 18, restricted mean survival time under H0: 10, under H1: 10            
#> Overall power: 0.025, overall significance level (1-sided): 0.025                  
#> Maximum # events: 342.1, expected # events: 341                                    
#> Maximum # subjects: 519.9, expected # subjects: 519.9                              
#> Maximum # milestone subjects: 156.2, expected # milestone subjects: 154.6          
#> Maximum information: 10.99, expected information: 10.97                            
#> Total study duration: 42, expected study duration: 41.8                            
#> Accrual duration: 24, follow-up duration: 18, fixed follow-up: TRUE                
#> Alpha spending: Lan-DeMets O'Brien-Fleming, beta spending: None                    
#>                                                                                    
#>                              Stage 1 Stage 2
#> Information rate             0.800   1.000  
#> Efficacy boundary (Z)        2.250   2.025  
#> Cumulative rejection         0.0122  0.0250 
#> Cumulative alpha spent       0.0122  0.0250 
#> Number of events             252.6   342.1  
#> Number of dropouts           14.6    21.5   
#> Number of subjects           519.9   519.9  
#> Number of milestone subjects 23.5    156.2  
#> Analysis time                24.9    42.0   
#> Efficacy boundary (rmst)     10.759  10.611 
#> Efficacy boundary (p)        0.0122  0.0214 
#> Information                  8.79    10.99  
#> 
```
