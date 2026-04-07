# Sample Size for One-Sample Milestone Survival Probability

Obtains the needed accrual duration given power and follow-up time, the
needed follow-up time given power and accrual duration, or the needed
absolute accrual rates given power, accrual duration, follow-up
duration, and relative accrual rates in a one-group survival design.

## Usage

``` r
kmsamplesize1s(
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
  survH0 = NA_real_,
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

- survH0:

  The milestone survival probability under the null hypothesis.

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

- `resultsUnderH1`: An S3 class `kmpower1s` object under the alternative
  hypothesis.

- `resultsUnderH0`: An S3 class `kmpower1s` object under the null
  hypothesis.

## See also

[`kmpower1s`](https://kaifenglu.github.io/lrstat/reference/kmpower1s.md)

## Author

Kaifeng Lu, <kaifenglu@gmail.com>

## Examples

``` r
# Example 1: Obtains follow-up duration given power, accrual intensity,
# and accrual duration for variable follow-up

kmsamplesize1s(beta = 0.2, kMax = 2,
               informationRates = c(0.8, 1),
               alpha = 0.025, typeAlphaSpending = "sfOF",
               milestone = 18, survH0 = 0.30,
               accrualTime = seq(0, 8),
               accrualIntensity = 26/9*seq(1, 9),
               piecewiseSurvivalTime = c(0, 6),
               stratumFraction = c(0.2, 0.8),
               lambda = c(0.0533, 0.0309, 1.5*0.0533, 1.5*0.0309),
               gamma = -log(1-0.05)/12, accrualDuration = 22,
               followupTime = NA, fixedFollowup = FALSE)
#> $resultsUnderH1
#>                                                                                     
#> Group-sequential design with 2 stages for one-sample milestone survival probability 
#> Milestone: 18, survival probability under H0: 0.3, under H1: 0.384                  
#> Overall power: 0.8, overall significance level (1-sided): 0.025                     
#> Maximum # events: 240.2, expected # events: 228.8                                   
#> Maximum # subjects: 468, expected # subjects: 468                                   
#> Maximum # milestone subjects: 41.9, expected # milestone subjects: 33               
#> Maximum information: 1130.95, expected information: 992.77                          
#> Total study duration: 26.5, expected study duration: 25.5                           
#> Accrual duration: 22, follow-up duration: 4.5, fixed follow-up: FALSE               
#> Alpha spending: Lan-DeMets O'Brien-Fleming, beta spending: None                     
#>                                                                                     
#>                              Stage 1 Stage 2
#> Information rate             0.800   1.000  
#> Efficacy boundary (Z)        2.250   2.025  
#> Cumulative rejection         0.6109  0.8000 
#> Cumulative alpha spent       0.0122  0.0250 
#> Number of events             221.6   240.2  
#> Number of dropouts           15.4    17.0   
#> Number of subjects           468.0   468.0  
#> Number of milestone subjects 27.3    41.9   
#> Analysis time                24.8    26.5   
#> Efficacy boundary (surv)     0.375   0.360  
#> Efficacy boundary (p)        0.0122  0.0214 
#> Information                  904.76  1130.95
#> 
#> $resultsUnderH0
#>                                                                                     
#> Group-sequential design with 2 stages for one-sample milestone survival probability 
#> Milestone: 18, survival probability under H0: 0.3, under H1: 0.3                    
#> Overall power: 0.025, overall significance level (1-sided): 0.025                   
#> Maximum # events: 268.9, expected # events: 268.6                                   
#> Maximum # subjects: 468, expected # subjects: 468                                   
#> Maximum # milestone subjects: 27.4, expected # milestone subjects: 27.3             
#> Maximum information: 1130.95, expected information: 1128.19                         
#> Total study duration: 25.8, expected study duration: 25.7                           
#> Accrual duration: 22, follow-up duration: 3.8, fixed follow-up: FALSE               
#> Alpha spending: Lan-DeMets O'Brien-Fleming, beta spending: None                     
#>                                                                                     
#>                              Stage 1 Stage 2
#> Information rate             0.800   1.000  
#> Efficacy boundary (Z)        2.250   2.025  
#> Cumulative rejection         0.0122  0.0250 
#> Cumulative alpha spent       0.0122  0.0250 
#> Number of events             250.3   268.9  
#> Number of dropouts           13.5    14.8   
#> Number of subjects           468.0   468.0  
#> Number of milestone subjects 18.4    27.4   
#> Analysis time                24.3    25.8   
#> Efficacy boundary (surv)     0.375   0.360  
#> Efficacy boundary (p)        0.0122  0.0214 
#> Information                  904.76  1130.95
#> 

# Example 2: Obtains accrual intensity given power, accrual duration, and
# follow-up duration for variable follow-up

kmsamplesize1s(beta = 0.2, kMax = 2,
               informationRates = c(0.8, 1),
               alpha = 0.025, typeAlphaSpending = "sfOF",
               milestone = 18, survH0 = 0.30,
               accrualTime = seq(0, 8),
               accrualIntensity = 26/9*seq(1, 9),
               piecewiseSurvivalTime = c(0, 6),
               stratumFraction = c(0.2, 0.8),
               lambda = c(0.0533, 0.0309, 1.5*0.0533, 1.5*0.0309),
               gamma = -log(1-0.05)/12, accrualDuration = 22,
               followupTime = 18, fixedFollowup = FALSE)
#> $resultsUnderH1
#>                                                                                     
#> Group-sequential design with 2 stages for one-sample milestone survival probability 
#> Milestone: 18, survival probability under H0: 0.3, under H1: 0.384                  
#> Overall power: 0.8, overall significance level (1-sided): 0.025                     
#> Maximum # events: 190.9, expected # events: 172.4                                   
#> Maximum # subjects: 275, expected # subjects: 275                                   
#> Maximum # milestone subjects: 92.4, expected # milestone subjects: 64.3             
#> Maximum information: 1130.95, expected information: 992.77                          
#> Total study duration: 39, expected study duration: 33.8                             
#> Accrual duration: 22, follow-up duration: 17, fixed follow-up: FALSE                
#> Alpha spending: Lan-DeMets O'Brien-Fleming, beta spending: None                     
#>                                                                                     
#>                              Stage 1 Stage 2
#> Information rate             0.800   1.000  
#> Efficacy boundary (Z)        2.250   2.025  
#> Cumulative rejection         0.6109  0.8000 
#> Cumulative alpha spent       0.0122  0.0250 
#> Number of events             160.7   190.9  
#> Number of dropouts           11.9    15.0   
#> Number of subjects           275.0   275.0  
#> Number of milestone subjects 46.4    92.4   
#> Analysis time                30.5    39.0   
#> Efficacy boundary (surv)     0.375   0.360  
#> Efficacy boundary (p)        0.0122  0.0214 
#> Information                  904.76  1130.95
#> 
#> $resultsUnderH0
#>                                                                                     
#> Group-sequential design with 2 stages for one-sample milestone survival probability 
#> Milestone: 18, survival probability under H0: 0.3, under H1: 0.3                    
#> Overall power: 0.025, overall significance level (1-sided): 0.025                   
#> Maximum # events: 193, expected # events: 192.8                                     
#> Maximum # subjects: 275, expected # subjects: 275                                   
#> Maximum # milestone subjects: 46.9, expected # milestone subjects: 46.7             
#> Maximum information: 1130.95, expected information: 1128.19                         
#> Total study duration: 33.1, expected study duration: 33                             
#> Accrual duration: 22, follow-up duration: 11.1, fixed follow-up: FALSE              
#> Alpha spending: Lan-DeMets O'Brien-Fleming, beta spending: None                     
#>                                                                                     
#>                              Stage 1 Stage 2
#> Information rate             0.800   1.000  
#> Efficacy boundary (Z)        2.250   2.025  
#> Cumulative rejection         0.0122  0.0250 
#> Cumulative alpha spent       0.0122  0.0250 
#> Number of events             175.5   193.0  
#> Number of dropouts           10.0    11.4   
#> Number of subjects           275.0   275.0  
#> Number of milestone subjects 29.4    46.9   
#> Analysis time                28.9    33.1   
#> Efficacy boundary (surv)     0.375   0.360  
#> Efficacy boundary (p)        0.0122  0.0214 
#> Information                  904.76  1130.95
#> 


# Example 3: Obtains accrual duration given power, accrual intensity, and
# follow-up duration for fixed follow-up

kmsamplesize1s(beta = 0.2, kMax = 2,
               informationRates = c(0.8, 1),
               alpha = 0.025, typeAlphaSpending = "sfOF",
               milestone = 18, survH0 = 0.30,
               accrualTime = seq(0, 8),
               accrualIntensity = 26/9*seq(1, 9),
               piecewiseSurvivalTime = c(0, 6),
               stratumFraction = c(0.2, 0.8),
               lambda = c(0.0533, 0.0309, 1.5*0.0533, 1.5*0.0309),
               gamma = -log(1-0.05)/12, accrualDuration = NA,
               followupTime = 18, fixedFollowup = TRUE)
#> $resultsUnderH1
#>                                                                                     
#> Group-sequential design with 2 stages for one-sample milestone survival probability 
#> Milestone: 18, survival probability under H0: 0.3, under H1: 0.384                  
#> Overall power: 0.8, overall significance level (1-sided): 0.025                     
#> Maximum # events: 164.7, expected # events: 158.9                                   
#> Maximum # subjects: 275, expected # subjects: 275                                   
#> Maximum # milestone subjects: 90.8, expected # milestone subjects: 57.7             
#> Maximum information: 1130.95, expected information: 992.77                          
#> Total study duration: 31.8, expected study duration: 28.2                           
#> Accrual duration: 14.6, follow-up duration: 18, fixed follow-up: TRUE               
#> Alpha spending: Lan-DeMets O'Brien-Fleming, beta spending: None                     
#>                                                                                     
#>                              Stage 1 Stage 2
#> Information rate             0.800   1.000  
#> Efficacy boundary (Z)        2.250   2.025  
#> Cumulative rejection         0.6109  0.8000 
#> Cumulative alpha spent       0.0122  0.0250 
#> Number of events             155.2   164.7  
#> Number of dropouts           11.4    12.3   
#> Number of subjects           275.0   275.0  
#> Number of milestone subjects 36.5    90.8   
#> Analysis time                25.9    31.8   
#> Efficacy boundary (surv)     0.375   0.360  
#> Efficacy boundary (p)        0.0122  0.0214 
#> Information                  904.76  1130.95
#> 
#> $resultsUnderH0
#>                                                                                     
#> Group-sequential design with 2 stages for one-sample milestone survival probability 
#> Milestone: 18, survival probability under H0: 0.3, under H1: 0.3                    
#> Overall power: 0.025, overall significance level (1-sided): 0.025                   
#> Maximum # events: 166.1, expected # events: 166                                     
#> Maximum # subjects: 243.4, expected # subjects: 243.4                               
#> Maximum # milestone subjects: 67.6, expected # milestone subjects: 67.1             
#> Maximum information: 1130.95, expected information: 1128.19                         
#> Total study duration: 31.4, expected study duration: 31.3                           
#> Accrual duration: 13.4, follow-up duration: 18, fixed follow-up: TRUE               
#> Alpha spending: Lan-DeMets O'Brien-Fleming, beta spending: None                     
#>                                                                                     
#>                              Stage 1 Stage 2
#> Information rate             0.800   1.000  
#> Efficacy boundary (Z)        2.250   2.025  
#> Cumulative rejection         0.0122  0.0250 
#> Cumulative alpha spent       0.0122  0.0250 
#> Number of events             158.5   166.1  
#> Number of dropouts           9.1     9.7    
#> Number of subjects           243.4   243.4  
#> Number of milestone subjects 25.3    67.6   
#> Analysis time                25.4    31.4   
#> Efficacy boundary (surv)     0.375   0.360  
#> Efficacy boundary (p)        0.0122  0.0214 
#> Information                  904.76  1130.95
#> 
```
