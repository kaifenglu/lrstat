# Sample Size for Equivalence in Hazard Ratio

Obtains the sample size for equivalence in hazard ratio.

## Usage

``` r
lrsamplesizeequiv(
  beta = 0.2,
  kMax = 1L,
  informationRates = NA_real_,
  criticalValues = NA_real_,
  alpha = 0.05,
  typeAlphaSpending = "sfOF",
  parameterAlphaSpending = NA_real_,
  userAlphaSpending = NA_real_,
  hazardRatioLower = NA_real_,
  hazardRatioUpper = NA_real_,
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
  typeOfComputation = "direct",
  spendingTime = NA_real_,
  rounding = TRUE
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

- hazardRatioLower:

  The lower equivalence limit of hazard ratio.

- hazardRatioUpper:

  The upper equivalence limit of hazard ratio.

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

  Whether to round up sample size. Defaults to 1 for sample size
  rounding.

## Value

An S3 class `lrpowerequiv` object

## See also

[`lrpowerequiv`](https://kaifenglu.github.io/lrstat/reference/lrpowerequiv.md)

## Author

Kaifeng Lu, <kaifenglu@gmail.com>

## Examples

``` r
lrsamplesizeequiv(kMax = 2, informationRates = c(0.5, 1),
                  alpha = 0.05, typeAlphaSpending = "sfOF",
                  hazardRatioLower = 0.71, hazardRatioUpper = 1.4,
                  allocationRatioPlanned = 1, accrualTime = seq(0, 8),
                  accrualIntensity = 26/9*seq(1, 9),
                  piecewiseSurvivalTime = c(0, 6),
                  lambda1 = c(0.0533, 0.0533),
                  lambda2 = c(0.0533, 0.0533),
                  gamma1 = -log(1-0.05)/12,
                  gamma2 = -log(1-0.05)/12, accrualDuration = NA,
                  followupTime = 18, fixedFollowup = FALSE)
#>                                                                          
#> Group-sequential design with 2 stages for equivalence in hazard ratio    
#> Lower limit for hazard ratio: 0.71, upper limit for hazard ratio: 1.4    
#> Overall power: 0.8009, overall alpha: 0.05                               
#> Maximum # events: 300, expected # events: 300                            
#> Maximum # dropouts: 24.1, expected # dropouts: 24.1                      
#> Maximum # subjects: 420, expected # subjects: 420                        
#> Maximum information: 75, expected information: 75                        
#> Total study duration: 38.2, expected study duration: 38.2                
#> Accrual duration: 20.2, follow-up duration: 18.1, fixed follow-up: FALSE 
#> Allocation ratio: 1                                                      
#> Alpha spending: Lan-DeMets O'Brien-Fleming                               
#>                                                                          
#>                                        Stage 1 Stage 2
#> Information rate                       0.500   1.000  
#> Boundary for each 1-sided test (Z)     2.538   1.662  
#> Cumulative rejection                   0.0000  0.8009 
#> Cumulative alpha for each 1-sided test 0.0056  0.0500 
#> Cumulative alpha attained under H10    0.0000  0.0500 
#> Cumulative alpha attained under H20    0.0000  0.0500 
#> Number of events                       150.0   300.0  
#> Number of dropouts                     12.0    24.1   
#> Number of subjects                     420.0   420.0  
#> Analysis time                          21.0    38.2   
#> Boundary for lower limit (HR)          1.075   0.860  
#> Boundary for upper limit (HR)          0.925   1.156  
#> Boundary for each 1-sided test (p)     0.0056  0.0482 
#> Information                            37.50   75.00  
#> HR                                     1.00    1.00   
```
