# Sample Size for Equivalence in Milestone Survival Probability Difference

Obtains the sample size for equivalence in milestone survival
probability difference.

## Usage

``` r
kmsamplesizeequiv(
  beta = 0.2,
  kMax = 1L,
  informationRates = NA_real_,
  criticalValues = NA_real_,
  alpha = 0.05,
  typeAlphaSpending = "sfOF",
  parameterAlphaSpending = NA_real_,
  userAlphaSpending = NA_real_,
  milestone = NA_real_,
  survDiffLower = NA_real_,
  survDiffUpper = NA_real_,
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
  fixedFollowup = 0L,
  spendingTime = NA_real_,
  rounding = 1L
)
```

## Arguments

- beta:

  The type II error.

- kMax:

  The maximum number of stages.

- informationRates:

  The information rates. Defaults to `(1:kMax) / kMax` if left
  unspecified.

- criticalValues:

  Upper boundaries on the z-test statistic scale for stopping for
  efficacy.

- alpha:

  The significance level for each of the two one-sided tests. Defaults
  to 0.05.

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

- milestone:

  The milestone time at which to calculate the survival probability.

- survDiffLower:

  The lower equivalence limit of milestone survival probability
  difference.

- survDiffUpper:

  The upper equivalence limit of milestone survival probability
  difference.

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

An S3 class `kmpowerequiv` object

## See also

[`kmpowerequiv`](https://kaifenglu.github.io/lrstat/reference/kmpowerequiv.md)

## Author

Kaifeng Lu, <kaifenglu@gmail.com>

## Examples

``` r
kmsamplesizeequiv(beta = 0.1, kMax = 2, informationRates = c(0.5, 1),
                  alpha = 0.05, typeAlphaSpending = "sfOF",
                  milestone = 18,
                  survDiffLower = -0.13, survDiffUpper = 0.13,
                  allocationRatioPlanned = 1, accrualTime = seq(0, 8),
                  accrualIntensity = 26/9*seq(1, 9),
                  piecewiseSurvivalTime = c(0, 6),
                  stratumFraction = c(0.2, 0.8),
                  lambda1 = c(0.0533, 0.0533, 1.5*0.0533, 1.5*0.0533),
                  lambda2 = c(0.0533, 0.0533, 1.5*0.0533, 1.5*0.0533),
                  gamma1 = -log(1-0.05)/12,
                  gamma2 = -log(1-0.05)/12, accrualDuration = NA,
                  followupTime = 18, fixedFollowup = FALSE)
#>                                                                                        
#> Group-sequential design with 2 stages for equivalence in milestone survival difference 
#> Milestone: 18, lower limit for survival difference: -0.13, upper limit: 0.13           
#> Milestone survival on treatment: 0.266, on control: 0.266, difference: 0               
#> Overall power: 0.9, overall alpha: 0.05                                                
#> Maximum # events: 426.9, expected # events: 426.9                                      
#> Maximum # subjects: 519, expected # subjects: 519                                      
#> Maximum information: 644.66, expected information: 644.66                              
#> Total study duration: 41.6, expected study duration: 41.6                              
#> Accrual duration: 24, follow-up duration: 17.7, fixed follow-up: FALSE                 
#> Allocation ratio: 1                                                                    
#> Alpha spending: Lan-DeMets O'Brien-Fleming                                             
#>                                                                                        
#>                                        Stage 1 Stage 2
#> Information rate                       0.500   1.000  
#> Boundary for each 1-sided test (Z)     2.538   1.662  
#> Cumulative rejection                   0.0000  0.9000 
#> Cumulative alpha for each 1-sided test 0.0056  0.0500 
#> Cumulative alpha attained under H10    0.0000  0.0500 
#> Cumulative alpha attained under H20    0.0000  0.0500 
#> Number of events                       295.5   426.9  
#> Number of dropouts                     17.1    24.8   
#> Number of subjects                     519.0   519.0  
#> Number of milestone subjects           32.0    125.9  
#> Analysis time                          27.0    41.6   
#> Boundary for lower limit (surv diff)   0.011   -0.065 
#> Boundary for upper limit (surv diff)   -0.011  0.065  
#> Boundary for each 1-sided test (p)     0.0056  0.0482 
#> Information                            322.33  644.66 
```
