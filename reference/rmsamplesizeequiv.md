# Sample Size for Equivalence in Restricted Mean Survival Time Difference

Obtains the sample size for equivalence in restricted mean survival time
difference.

## Usage

``` r
rmsamplesizeequiv(
  beta = 0.2,
  kMax = 1L,
  informationRates = NA_real_,
  criticalValues = NULL,
  alpha = 0.05,
  typeAlphaSpending = "sfOF",
  parameterAlphaSpending = NA_real_,
  userAlphaSpending = NA_real_,
  milestone = NA_real_,
  rmstDiffLower = NA_real_,
  rmstDiffUpper = NA_real_,
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
  for `"WT"`, \\\rho\\ for `"sfKD"`, and \\\gamma\\ for `"sfHSD"`.

- userAlphaSpending:

  The user defined alpha spending. Cumulative alpha spent up to each
  stage.

- milestone:

  The milestone time at which to calculate the restricted mean survival
  time.

- rmstDiffLower:

  The lower equivalence limit of restricted mean survival time
  difference.

- rmstDiffUpper:

  The upper equivalence limit of restricted mean survival time
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

An S3 class `rmpowerequiv` object

## See also

[`rmpowerequiv`](https://kaifenglu.github.io/lrstat/reference/rmpowerequiv.md)

## Author

Kaifeng Lu, <kaifenglu@gmail.com>

## Examples

``` r
rmsamplesizeequiv(beta = 0.1, kMax = 2, informationRates = c(0.5, 1),
                  alpha = 0.05, typeAlphaSpending = "sfOF",
                  milestone = 18,
                  rmstDiffLower = -2, rmstDiffUpper = 2,
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
#> Group-sequential design with 2 stages for equivalence in RMST difference 
#> Milestone: 18, lower limit for RMST difference: -2, upper limit: 2       
#> RMST on treatment: 9.948, on control: 9.948, difference: 0               
#> Overall power: 0.9, overall alpha: 0.05                                  
#> Maximum # events: 345.5, expected # events: 345.5                        
#> Maximum # subjects: 456, expected # subjects: 456                        
#> Maximum information: 2.72, expected information: 2.72                    
#> Total study duration: 34.5, expected study duration: 34.5                
#> Accrual duration: 21.5, follow-up duration: 13, fixed follow-up: FALSE   
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
#> Number of events                       159.2   345.5  
#> Number of dropouts                     9.2     20.1   
#> Number of subjects                     398.3   456.0  
#> Number of milestone subjects           1.2     80.2   
#> Analysis time                          19.3    34.5   
#> Boundary for lower limit (rmst diff)   0.175   -0.993 
#> Boundary for upper limit (rmst diff)   -0.175  0.993  
#> Boundary for each 1-sided test (p)     0.0056  0.0482 
#> Information                            1.36    2.72   
```
