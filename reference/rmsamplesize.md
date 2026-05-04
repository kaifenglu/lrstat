# Sample Size for Difference in Restricted Mean Survival Times

Obtains the needed accrual duration given power, accrual intensity, and
follow-up time, the needed follow-up time given power, accrual
intensity, and accrual duration, or the needed absolute accrual
intensity given power, relative accrual intensity, accrual duration, and
follow-up time in a two-group survival design.

## Usage

``` r
rmsamplesize(
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
  futilityRmstDiff = NULL,
  typeBetaSpending = "none",
  parameterBetaSpending = NA_real_,
  userBetaSpending = NA_real_,
  milestone = NA_real_,
  rmstDiffH0 = 0,
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

- futilityRmstDiff:

  A vector of length `kMax - 1` for the futility bounds on the
  restricted mean survival time difference scale.

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

  The milestone time at which to calculate the restricted mean survival
  time.

- rmstDiffH0:

  The difference in restricted mean survival times under the null
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

- `resultsUnderH1`: An S3 class `rmpower` object under the alternative
  hypothesis.

- `resultsUnderH0`: An S3 class `rmpower` object under the null
  hypothesis.

## See also

[`rmpower`](https://kaifenglu.github.io/lrstat/reference/rmpower.md)

## Author

Kaifeng Lu, <kaifenglu@gmail.com>

## Examples

``` r
# Example 1: Obtains follow-up time given power, accrual intensity,
# and accrual duration for variable follow-up. Of note, the power
# reaches the maximum when the follow-up time equals milestone.

rmsamplesize(beta = 0.2, kMax = 2, informationRates = c(0.8, 1),
             alpha = 0.025, typeAlphaSpending = "sfOF",
             milestone = 18,
             allocationRatioPlanned = 1, accrualTime = seq(0, 8),
             accrualIntensity = 100/9*seq(1, 9),
             piecewiseSurvivalTime = c(0, 6),
             stratumFraction = c(0.2, 0.8),
             lambda1 = c(0.0533, 0.0309, 1.5*0.0533, 1.5*0.0309),
             lambda2 = c(0.0533, 0.0533, 1.5*0.0533, 1.5*0.0533),
             gamma1 = -log(1-0.05)/12,
             gamma2 = -log(1-0.05)/12, accrualDuration = 22,
             followupTime = NA, fixedFollowup = FALSE)
#> $resultsUnderH1
#>                                                                                       
#> Group-sequential design with 2 stages for difference in restricted mean survival time 
#> Milestone: 18, restricted mean survival time difference under H0: 0                   
#> Restricted mean survival time on treatment: 10.854, on control: 9.948                 
#> Overall power: 0.8, overall significance level (1-sided): 0.025                       
#> Maximum # events: 1107.6, expected # events: 939.4                                    
#> Maximum # subjects: 1800, expected # subjects: 1800                                   
#> Maximum information: 9.77, expected information: 8.57                                 
#> Total study duration: 29.5, expected study duration: 25.7                             
#> Accrual duration: 22, follow-up duration: 7.5, fixed follow-up: FALSE                 
#> Allocation ratio: 1                                                                   
#> Alpha spending: Lan-DeMets O'Brien-Fleming, beta spending: None                       
#>                                                                                       
#>                               Stage 1 Stage 2
#> Information rate              0.800   1.000  
#> Efficacy boundary (Z)         2.250   2.025  
#> Cumulative rejection          0.6109  0.8000 
#> Cumulative alpha spent        0.0122  0.0250 
#> Number of events              832.3   1107.6 
#> Number of dropouts            52.2    72.1   
#> Number of subjects            1800.0  1800.0 
#> Number of milestone subjects  56.7    225.5  
#> Analysis time                 23.3    29.5   
#> Efficacy boundary (rmst diff) 0.805   0.648  
#> Efficacy boundary (p)         0.0122  0.0214 
#> Information                   7.81    9.77   
#> 
#> $resultsUnderH0
#>                                                                                       
#> Group-sequential design with 2 stages for difference in restricted mean survival time 
#> Milestone: 18, restricted mean survival time difference under H0: 0                   
#> Restricted mean survival time on treatment: 9.948, on control: 9.948                  
#> Overall power: 0.025, overall significance level (1-sided): 0.025                     
#> Maximum # events: 1079, expected # events: 1076.2                                     
#> Maximum # subjects: 1800, expected # subjects: 1800                                   
#> Maximum information: 9.77, expected information: 9.74                                 
#> Total study duration: 26.8, expected study duration: 26.8                             
#> Accrual duration: 22, follow-up duration: 4.8, fixed follow-up: FALSE                 
#> Allocation ratio: 1                                                                   
#> Alpha spending: Lan-DeMets O'Brien-Fleming, beta spending: None                       
#>                                                                                       
#>                               Stage 1 Stage 2
#> Information rate              0.800   1.000  
#> Efficacy boundary (Z)         2.250   2.025  
#> Cumulative rejection          0.0122  0.0250 
#> Cumulative alpha spent        0.0122  0.0250 
#> Number of events              848.5   1079.0 
#> Number of dropouts            49.0    62.4   
#> Number of subjects            1800.0  1800.0 
#> Number of milestone subjects  38.1    119.3  
#> Analysis time                 22.8    26.8   
#> Efficacy boundary (rmst diff) 0.805   0.648  
#> Efficacy boundary (p)         0.0122  0.0214 
#> Information                   7.81    9.77   
#> 

# Example 2: Obtains accrual intensity given power, accrual duration, and
# follow-up time for variable follow-up

rmsamplesize(beta = 0.2, kMax = 2, informationRates = c(0.8, 1),
             alpha = 0.025, typeAlphaSpending = "sfOF",
             milestone = 18,
             allocationRatioPlanned = 1, accrualTime = seq(0, 8),
             accrualIntensity = 100/9*seq(1, 9),
             piecewiseSurvivalTime = c(0, 6),
             stratumFraction = c(0.2, 0.8),
             lambda1 = c(0.0533, 0.0309, 1.5*0.0533, 1.5*0.0309),
             lambda2 = c(0.0533, 0.0533, 1.5*0.0533, 1.5*0.0533),
             gamma1 = -log(1-0.05)/12,
             gamma2 = -log(1-0.05)/12, accrualDuration = 22,
             followupTime = 18, fixedFollowup = FALSE)
#> $resultsUnderH1
#>                                                                                       
#> Group-sequential design with 2 stages for difference in restricted mean survival time 
#> Milestone: 18, restricted mean survival time difference under H0: 0                   
#> Restricted mean survival time on treatment: 10.854, on control: 9.948                 
#> Overall power: 0.8, overall significance level (1-sided): 0.025                       
#> Maximum # events: 1250.7, expected # events: 995                                      
#> Maximum # subjects: 1745, expected # subjects: 1745                                   
#> Maximum information: 9.77, expected information: 8.57                                 
#> Total study duration: 36, expected study duration: 28.5                               
#> Accrual duration: 22, follow-up duration: 14, fixed follow-up: FALSE                  
#> Allocation ratio: 1                                                                   
#> Alpha spending: Lan-DeMets O'Brien-Fleming, beta spending: None                       
#>                                                                                       
#>                               Stage 1 Stage 2
#> Information rate              0.800   1.000  
#> Efficacy boundary (Z)         2.250   2.025  
#> Cumulative rejection          0.6109  0.8000 
#> Cumulative alpha spent        0.0122  0.0250 
#> Number of events              832.1   1250.7 
#> Number of dropouts            52.3    83.8   
#> Number of subjects            1745.0  1745.0 
#> Number of milestone subjects  63.8    408.1  
#> Analysis time                 23.8    36.0   
#> Efficacy boundary (rmst diff) 0.805   0.648  
#> Efficacy boundary (p)         0.0122  0.0214 
#> Information                   7.81    9.77   
#> 
#> $resultsUnderH0
#>                                                                                       
#> Group-sequential design with 2 stages for difference in restricted mean survival time 
#> Milestone: 18, restricted mean survival time difference under H0: 0                   
#> Restricted mean survival time on treatment: 9.948, on control: 9.948                  
#> Overall power: 0.025, overall significance level (1-sided): 0.025                     
#> Maximum # events: 1095.4, expected # events: 1092.4                                   
#> Maximum # subjects: 1745, expected # subjects: 1745                                   
#> Maximum information: 9.77, expected information: 9.74                                 
#> Total study duration: 27.9, expected study duration: 27.9                             
#> Accrual duration: 22, follow-up duration: 5.9, fixed follow-up: FALSE                 
#> Allocation ratio: 1                                                                   
#> Alpha spending: Lan-DeMets O'Brien-Fleming, beta spending: None                       
#>                                                                                       
#>                               Stage 1 Stage 2
#> Information rate              0.800   1.000  
#> Efficacy boundary (Z)         2.250   2.025  
#> Cumulative rejection          0.0122  0.0250 
#> Cumulative alpha spent        0.0122  0.0250 
#> Number of events              847.1   1095.4 
#> Number of dropouts            48.9    63.4   
#> Number of subjects            1745.0  1745.0 
#> Number of milestone subjects  42.5    142.0  
#> Analysis time                 23.2    27.9   
#> Efficacy boundary (rmst diff) 0.805   0.648  
#> Efficacy boundary (p)         0.0122  0.0214 
#> Information                   7.81    9.77   
#> 


# Example 3: Obtains accrual duration given power, accrual intensity, and
# follow-up time for fixed follow-up

rmsamplesize(beta = 0.2, kMax = 2, informationRates = c(0.8, 1),
             alpha = 0.025, typeAlphaSpending = "sfOF",
             milestone = 18,
             allocationRatioPlanned = 1, accrualTime = seq(0, 8),
             accrualIntensity = 100/9*seq(1, 9),
             piecewiseSurvivalTime = c(0, 6),
             stratumFraction = c(0.2, 0.8),
             lambda1 = c(0.0533, 0.0309, 1.5*0.0533, 1.5*0.0309),
             lambda2 = c(0.0533, 0.0533, 1.5*0.0533, 1.5*0.0533),
             gamma1 = -log(1-0.05)/12,
             gamma2 = -log(1-0.05)/12, accrualDuration = NA,
             followupTime = 18, fixedFollowup = TRUE)
#> $resultsUnderH1
#>                                                                                       
#> Group-sequential design with 2 stages for difference in restricted mean survival time 
#> Milestone: 18, restricted mean survival time difference under H0: 0                   
#> Restricted mean survival time on treatment: 10.854, on control: 9.948                 
#> Overall power: 0.8, overall significance level (1-sided): 0.025                       
#> Maximum # events: 1130, expected # events: 944                                        
#> Maximum # subjects: 1745, expected # subjects: 1745                                   
#> Maximum information: 9.77, expected information: 8.57                                 
#> Total study duration: 35.5, expected study duration: 28.1                             
#> Accrual duration: 21.4, follow-up duration: 18, fixed follow-up: TRUE                 
#> Allocation ratio: 1                                                                   
#> Alpha spending: Lan-DeMets O'Brien-Fleming, beta spending: None                       
#>                                                                                       
#>                               Stage 1 Stage 2
#> Information rate              0.800   1.000  
#> Efficacy boundary (Z)         2.250   2.025  
#> Cumulative rejection          0.6109  0.8000 
#> Cumulative alpha spent        0.0122  0.0250 
#> Number of events              825.6   1130.0 
#> Number of dropouts            51.8    74.1   
#> Number of subjects            1745.0  1745.0 
#> Number of milestone subjects  58.6    405.5  
#> Analysis time                 23.4    35.5   
#> Efficacy boundary (rmst diff) 0.805   0.648  
#> Efficacy boundary (p)         0.0122  0.0214 
#> Information                   7.81    9.77   
#> 
#> $resultsUnderH0
#>                                                                                       
#> Group-sequential design with 2 stages for difference in restricted mean survival time 
#> Milestone: 18, restricted mean survival time difference under H0: 0                   
#> Restricted mean survival time on treatment: 9.948, on control: 9.948                  
#> Overall power: 0.025, overall significance level (1-sided): 0.025                     
#> Maximum # events: 1162.6, expected # events: 1158.7                                   
#> Maximum # subjects: 1632.7, expected # subjects: 1632.7                               
#> Maximum information: 9.77, expected information: 9.74                                 
#> Total study duration: 38.3, expected study duration: 38.1                             
#> Accrual duration: 20.3, follow-up duration: 18, fixed follow-up: TRUE                 
#> Allocation ratio: 1                                                                   
#> Alpha spending: Lan-DeMets O'Brien-Fleming, beta spending: None                       
#>                                                                                       
#>                               Stage 1 Stage 2
#> Information rate              0.800   1.000  
#> Efficacy boundary (Z)         2.250   2.025  
#> Cumulative rejection          0.0122  0.0250 
#> Cumulative alpha spent        0.0122  0.0250 
#> Number of events              839.5   1162.6 
#> Number of dropouts            48.5    67.4   
#> Number of subjects            1632.7  1632.7 
#> Number of milestone subjects  41.9    402.6  
#> Analysis time                 23.1    38.3   
#> Efficacy boundary (rmst diff) 0.805   0.648  
#> Efficacy boundary (p)         0.0122  0.0214 
#> Information                   7.81    9.77   
#> 
```
