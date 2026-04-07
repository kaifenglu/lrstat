# Sample Size for Negative Binomial Rate Ratio

Obtains the needed accrual duration given power and follow-up time, the
needed follow-up time given power and accrual duration, or the needed
absolute accrual rates given power, accrual duration, follow-up
duration, and relative accrual rates in a two-group negative binomial
design.

## Usage

``` r
nbsamplesize(
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
  rounding = TRUE,
  nullVariance = FALSE
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

- rounding:

  Whether to round up sample size. Defaults to 1 for sample size
  rounding.

- nullVariance:

  Whether to calculate the variance for log rate ratio under the null
  hypothesis.

## Value

A list of two components:

- `resultsUnderH1`: An S3 class `nbpower` object under the alternative
  hypothesis.

- `resultsUnderH0`: An S3 class `nbpower` object under the null
  hypothesis.

## See also

[`nbpower`](https://github.com/kaifenglu/lrstat/reference/nbpower.md)

## Author

Kaifeng Lu, <kaifenglu@gmail.com>

## Examples

``` r
# Example 1: Obtains follow-up duration given power, accrual intensity,
# and accrual duration for variable follow-up

nbsamplesize(beta = 0.2, kMax = 2,
             informationRates = c(0.5, 1),
             alpha = 0.025, typeAlphaSpending = "sfOF",
             accrualIntensity = 1956/1.25,
             kappa1 = 5, kappa2 = 5,
             lambda1 = 0.0875, lambda2 = 0.125,
             gamma1 = 0, gamma2 = 0,
             accrualDuration = 1.25,
             followupTime = NA, fixedFollowup = FALSE)
#> $resultsUnderH1
#>                                                                        
#> Group-sequential design with 2 stages for negative binomial rate ratio 
#> Rate ratio under H0: 1, rate ratio under H1: 0.7                       
#> Event rate for treatment: 0.0875, event rate for control: 0.125        
#> Dispersion for treatment: 5, dispersion for control: 5                 
#> Overall power: 0.8, overall significance level (1-sided): 0.025        
#> Maximum # events: 702.1, expected # events: 619.1                      
#> Maximum # dropouts: 0, expected # dropouts: 0                          
#> Maximum # subjects: 1956, expected # subjects: 1956                    
#> Maximum exposure: 6607.9, expected exposure: 5827.1                    
#> Maximum information: 61.93, expected information: 56.85                
#> Total study duration: 4, expected study duration: 3.6                  
#> Accrual duration: 1.2, follow-up duration: 2.8, fixed follow-up: FALSE 
#> Allocation ratio: 1, variance of standardized test statistic: under H1 
#> Alpha spending: Lan-DeMets O'Brien-Fleming, beta spending: None        
#>                                                                        
#>                                Stage 1 Stage 2
#> Information rate               0.500   1.000  
#> Efficacy boundary (Z)          2.963   1.969  
#> Cumulative rejection           0.1641  0.8000 
#> Cumulative alpha spent         0.0015  0.0250 
#> Number of events               196.5   702.1  
#> Number of dropouts             0.0     0.0    
#> Number of subjects             1956.0  1956.0 
#> Exposure                       1849.3  6607.9 
#> Analysis time                  1.6     4.0    
#> Efficacy boundary (rate ratio) 0.587   0.779  
#> Efficacy boundary (p)          0.0015  0.0245 
#> Information                    30.96   61.93  
#> 
#> $resultsUnderH0
#>                                                                        
#> Group-sequential design with 2 stages for negative binomial rate ratio 
#> Rate ratio under H0: 1, rate ratio under H1: 1                         
#> Event rate for treatment: 0.125, event rate for control: 0.125         
#> Dispersion for treatment: 5, dispersion for control: 5                 
#> Overall power: 0.025, overall significance level (1-sided): 0.025      
#> Maximum # events: 682.6, expected # events: 681.9                      
#> Maximum # dropouts: 0, expected # dropouts: 0                          
#> Maximum # subjects: 1956, expected # subjects: 1956                    
#> Maximum exposure: 5460.8, expected exposure: 5454.8                    
#> Maximum information: 61.93, expected information: 61.88                
#> Total study duration: 3.4, expected study duration: 3.4                
#> Accrual duration: 1.2, follow-up duration: 2.2, fixed follow-up: FALSE 
#> Allocation ratio: 1, variance of standardized test statistic: under H1 
#> Alpha spending: Lan-DeMets O'Brien-Fleming, beta spending: None        
#>                                                                        
#>                                Stage 1 Stage 2
#> Information rate               0.500   1.000  
#> Efficacy boundary (Z)          2.963   1.969  
#> Cumulative rejection           0.0015  0.0250 
#> Cumulative alpha spent         0.0015  0.0250 
#> Number of events               194.8   682.6  
#> Number of dropouts             0.0     0.0    
#> Number of subjects             1956.0  1956.0 
#> Exposure                       1558.1  5460.8 
#> Analysis time                  1.4     3.4    
#> Efficacy boundary (rate ratio) 0.587   0.779  
#> Efficacy boundary (p)          0.0015  0.0245 
#> Information                    30.96   61.93  
#> 

# Example 2: Obtains accrual intensity given power, accrual duration, and
# follow-up duration for variable follow-up

nbsamplesize(beta = 0.2, kMax = 2,
             informationRates = c(0.5, 1),
             alpha = 0.025, typeAlphaSpending = "sfOF",
             accrualIntensity = 100,
             kappa1 = 5, kappa2 = 5,
             lambda1 = 0.0875, lambda2 = 0.125,
             gamma1 = 0, gamma2 = 0,
             accrualDuration = 1.25,
             followupTime = 2.25, fixedFollowup = FALSE)
#> $resultsUnderH1
#>                                                                        
#> Group-sequential design with 2 stages for negative binomial rate ratio 
#> Rate ratio under H0: 1, rate ratio under H1: 0.7                       
#> Event rate for treatment: 0.0875, event rate for control: 0.125        
#> Dispersion for treatment: 5, dispersion for control: 5                 
#> Overall power: 0.8, overall significance level (1-sided): 0.025        
#> Maximum # events: 636.2, expected # events: 563.3                      
#> Maximum # dropouts: 0, expected # dropouts: 0                          
#> Maximum # subjects: 2084, expected # subjects: 2084                    
#> Maximum exposure: 5987.5, expected exposure: 5301.9                    
#> Maximum information: 61.93, expected information: 56.85                
#> Total study duration: 3.5, expected study duration: 3.2                
#> Accrual duration: 1.2, follow-up duration: 2.2, fixed follow-up: FALSE 
#> Allocation ratio: 1, variance of standardized test statistic: under H1 
#> Alpha spending: Lan-DeMets O'Brien-Fleming, beta spending: None        
#>                                                                        
#>                                Stage 1 Stage 2
#> Information rate               0.500   1.000  
#> Efficacy boundary (Z)          2.963   1.969  
#> Cumulative rejection           0.1641  0.8000 
#> Cumulative alpha spent         0.0015  0.0250 
#> Number of events               192.2   636.2  
#> Number of dropouts             0.0     0.0    
#> Number of subjects             2084.0  2084.0 
#> Exposure                       1809.1  5987.5 
#> Analysis time                  1.5     3.5    
#> Efficacy boundary (rate ratio) 0.587   0.779  
#> Efficacy boundary (p)          0.0015  0.0245 
#> Information                    30.96   61.93  
#> 
#> $resultsUnderH0
#>                                                                        
#> Group-sequential design with 2 stages for negative binomial rate ratio 
#> Rate ratio under H0: 1, rate ratio under H1: 1                         
#> Event rate for treatment: 0.125, event rate for control: 0.125         
#> Dispersion for treatment: 5, dispersion for control: 5                 
#> Overall power: 0.025, overall significance level (1-sided): 0.025      
#> Maximum # events: 619.2, expected # events: 618.5                      
#> Maximum # dropouts: 0, expected # dropouts: 0                          
#> Maximum # subjects: 2084, expected # subjects: 2084                    
#> Maximum exposure: 4953.3, expected exposure: 4948.1                    
#> Maximum information: 61.93, expected information: 61.88                
#> Total study duration: 3, expected study duration: 3                    
#> Accrual duration: 1.2, follow-up duration: 1.8, fixed follow-up: FALSE 
#> Allocation ratio: 1, variance of standardized test statistic: under H1 
#> Alpha spending: Lan-DeMets O'Brien-Fleming, beta spending: None        
#>                                                                        
#>                                Stage 1 Stage 2
#> Information rate               0.500   1.000  
#> Efficacy boundary (Z)          2.963   1.969  
#> Cumulative rejection           0.0015  0.0250 
#> Cumulative alpha spent         0.0015  0.0250 
#> Number of events               191.0   619.2  
#> Number of dropouts             0.0     0.0    
#> Number of subjects             2084.0  2084.0 
#> Exposure                       1528.3  4953.3 
#> Analysis time                  1.4     3.0    
#> Efficacy boundary (rate ratio) 0.587   0.779  
#> Efficacy boundary (p)          0.0015  0.0245 
#> Information                    30.96   61.93  
#> 


# Example 3: Obtains accrual duration given power, accrual intensity, and
# follow-up duration for fixed follow-up

nbsamplesize(beta = 0.2, kMax = 2,
             informationRates = c(0.5, 1),
             alpha = 0.025, typeAlphaSpending = "sfOF",
             accrualIntensity = 1667,
             stratumFraction = c(0.2, 0.8),
             kappa1 = 5, kappa2 = 5,
             lambda1 = c(0.7*0.125, 0.75*0.25),
             lambda2 = c(0.125, 0.25),
             gamma1 = 0, gamma2 = 0,
             accrualDuration = NA,
             followupTime = 0.5, fixedFollowup = TRUE)
#> $resultsUnderH1
#>                                                                             
#> Group-sequential design with 2 stages for negative binomial rate ratio      
#> Rate ratio under H0: 1, rate ratio under H1: 0.74                           
#> Stratum fraction: 0.2 0.8                                                   
#> Event rate for treatment: 0.0875 0.1875, event rate for control: 0.125 0.25 
#> Dispersion for treatment: 5, dispersion for control: 5                      
#> Overall power: 0.8, overall significance level (1-sided): 0.025             
#> Maximum # events: 556.3, expected # events: 509.8                           
#> Maximum # dropouts: 0, expected # dropouts: 0                               
#> Maximum # subjects: 5670, expected # subjects: 5264.8                       
#> Maximum exposure: 2834.6, expected exposure: 2597.9                         
#> Maximum information: 86.68, expected information: 79.57                     
#> Total study duration: 3.9, expected study duration: 3.6                     
#> Accrual duration: 3.4, follow-up duration: 0.5, fixed follow-up: TRUE       
#> Allocation ratio: 1, variance of standardized test statistic: under H1      
#> Alpha spending: Lan-DeMets O'Brien-Fleming, beta spending: None             
#>                                                                             
#>                                Stage 1 Stage 2
#> Information rate               0.500   1.000  
#> Efficacy boundary (Z)          2.963   1.969  
#> Cumulative rejection           0.1641  0.8000 
#> Cumulative alpha spent         0.0015  0.0250 
#> Number of events               273.1   556.3  
#> Number of dropouts             0.0     0.0    
#> Number of subjects             3200.4  5670.0 
#> Exposure                       1391.8  2834.6 
#> Analysis time                  1.9     3.9    
#> Efficacy boundary (rate ratio) 0.638   0.809  
#> Efficacy boundary (p)          0.0015  0.0245 
#> Information                    43.34   86.68  
#> 
#> $resultsUnderH0
#>                                                                          
#> Group-sequential design with 2 stages for negative binomial rate ratio   
#> Rate ratio under H0: 1, rate ratio under H1: 1                           
#> Stratum fraction: 0.2 0.8                                                
#> Event rate for treatment: 0.125 0.25, event rate for control: 0.125 0.25 
#> Dispersion for treatment: 5, dispersion for control: 5                   
#> Overall power: 0.025, overall significance level (1-sided): 0.025        
#> Maximum # events: 569.5, expected # events: 569                          
#> Maximum # dropouts: 0, expected # dropouts: 0                            
#> Maximum # subjects: 5061.9, expected # subjects: 5058.6                  
#> Maximum exposure: 2531, expected exposure: 2529                          
#> Maximum information: 86.68, expected information: 86.61                  
#> Total study duration: 3.5, expected study duration: 3.5                  
#> Accrual duration: 3, follow-up duration: 0.5, fixed follow-up: TRUE      
#> Allocation ratio: 1, variance of standardized test statistic: under H1   
#> Alpha spending: Lan-DeMets O'Brien-Fleming, beta spending: None          
#>                                                                          
#>                                Stage 1 Stage 2
#> Information rate               0.500   1.000  
#> Efficacy boundary (Z)          2.963   1.969  
#> Cumulative rejection           0.0015  0.0250 
#> Cumulative alpha spent         0.0015  0.0250 
#> Number of events               278.2   569.5  
#> Number of dropouts             0.0     0.0    
#> Number of subjects             2889.3  5061.9 
#> Exposure                       1236.3  2531.0 
#> Analysis time                  1.7     3.5    
#> Efficacy boundary (rate ratio) 0.638   0.809  
#> Efficacy boundary (p)          0.0015  0.0245 
#> Information                    43.34   86.68  
#> 
```
