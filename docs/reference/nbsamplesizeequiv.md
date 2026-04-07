# Sample Size for Equivalence in Negative Binomial Rate Ratio

Obtains the sample size for equivalence in negative binomial rate ratio.

## Usage

``` r
nbsamplesizeequiv(
  beta = 0.2,
  kMax = 1L,
  informationRates = NA_real_,
  criticalValues = NA_real_,
  alpha = 0.05,
  typeAlphaSpending = "sfOF",
  parameterAlphaSpending = NA_real_,
  userAlphaSpending = NA_real_,
  rateRatioLower = NA_real_,
  rateRatioUpper = NA_real_,
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

- rateRatioLower:

  The lower equivalence limit of rate ratio.

- rateRatioUpper:

  The upper equivalence limit of rate ratio.

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

## Value

An S3 class `nbpowerequiv` object

## See also

[`nbpowerequiv`](https://kaifenglu.github.io/lrstat/reference/nbpowerequiv.md)

## Author

Kaifeng Lu, <kaifenglu@gmail.com>

## Examples

``` r
# Example 1: Variable follow-up design and solve for follow-up time
nbsamplesizeequiv(beta = 0.1, kMax = 2, informationRates = c(0.5, 1),
                  alpha = 0.05, typeAlphaSpending = "sfOF",
                  rateRatioLower = 2/3, rateRatioUpper = 3/2,
                  accrualIntensity = 1956/1.25,
                  stratumFraction = c(0.2, 0.8),
                  kappa1 = c(3, 5),
                  kappa2 = c(2, 3),
                  lambda1 = c(0.125, 0.165),
                  lambda2 = c(0.135, 0.175),
                  gamma1 = -log(1-0.05),
                  gamma2 = -log(1-0.10),
                  accrualDuration = 1.25,
                  followupTime = NA, fixedFollowup = FALSE)
#>                                                                                       
#> Group-sequential design with 2 stages for equivalence in negative binomial rate ratio 
#> Lower limit for rate ratio: 0.667, upper limit for rate ratio: 1.5, rate ratio: 0.939 
#> Stratum fraction: 0.2 0.8                                                             
#> Event rate for treatment: 0.125 0.165, event rate for control: 0.135 0.175            
#> Dispersion for treatment: 3 5, dispersion for control: 2 3                            
#> Overall power: 0.9, overall alpha: 0.05                                               
#> Maximum # events: 771.7, expected # events: 771.7                                     
#> Maximum # subjects: 1956, expected # subjects: 1956                                   
#> Maximum exposure: 4768.6, expected exposure: 4768.6                                   
#> Maximum information: 75.65, expected information: 75.65                               
#> Total study duration: 3.3, expected study duration: 3.3                               
#> Accrual duration: 1.2, follow-up duration: 2.1, fixed follow-up: FALSE                
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
#> Number of events                       231.7   771.7  
#> Number of dropouts                     111.6   368.9  
#> Number of subjects                     1956.0  1956.0 
#> Exposure                               1430.9  4768.6 
#> Analysis time                          1.4     3.3    
#> Boundary for lower limit (rate ratio)  1.007   0.807  
#> Boundary for upper limit (rate ratio)  0.993   1.239  
#> Boundary for each 1-sided test (p)     0.0056  0.0482 
#> Information                            37.83   75.65  

# Example 2: Fixed follow-up design and solve for accrual duration
nbsamplesizeequiv(beta = 0.2, kMax = 2, informationRates = c(0.5, 1),
                  alpha = 0.05, typeAlphaSpending = "sfOF",
                  rateRatioLower = 0.5, rateRatioUpper = 2,
                  accrualIntensity = 220/1.5,
                  kappa1 = 3, kappa2 = 3,
                  lambda1 = 8.4, lambda2 = 8.4,
                  gamma1 = 0, gamma2 = 0,
                  accrualDuration = NA,
                  followupTime = 0.5, fixedFollowup = TRUE)
#>                                                                                       
#> Group-sequential design with 2 stages for equivalence in negative binomial rate ratio 
#> Lower limit for rate ratio: 0.5, upper limit for rate ratio: 2, rate ratio: 1         
#> Event rate for treatment: 8.4, event rate for control: 8.4                            
#> Dispersion for treatment: 3, dispersion for control: 3                                
#> Overall power: 0.8, overall alpha: 0.05                                               
#> Maximum # events: 958, expected # events: 958                                         
#> Maximum # subjects: 233, expected # subjects: 233                                     
#> Maximum exposure: 114.1, expected exposure: 114.1                                     
#> Maximum information: 17.95, expected information: 17.95                               
#> Total study duration: 1.9, expected study duration: 1.9                               
#> Accrual duration: 1.6, follow-up duration: 0.5, fixed follow-up: TRUE                 
#> Allocation ratio: 1                                                                   
#> Alpha spending: Lan-DeMets O'Brien-Fleming                                            
#>                                                                                       
#>                                        Stage 1 Stage 2
#> Information rate                       0.500   1.000  
#> Boundary for each 1-sided test (Z)     2.538   1.662  
#> Cumulative rejection                   0.0000  0.8000 
#> Cumulative alpha for each 1-sided test 0.0056  0.0500 
#> Cumulative alpha attained under H10    0.0000  0.0500 
#> Cumulative alpha attained under H20    0.0000  0.0500 
#> Number of events                       378.7   958.0  
#> Number of dropouts                     0.0     0.0    
#> Number of subjects                     126.8   233.0  
#> Exposure                               45.1    114.1  
#> Analysis time                          0.9     1.9    
#> Boundary for lower limit (rate ratio)  1.166   0.740  
#> Boundary for upper limit (rate ratio)  0.857   1.351  
#> Boundary for each 1-sided test (p)     0.0056  0.0482 
#> Information                            8.98    17.95  
```
