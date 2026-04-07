# Sample Size for One-Sample Negative Binomial Rate

Obtains the needed accrual duration given power and follow-up time, the
needed follow-up time given power and accrual duration, or the needed
absolute accrual rates given power, accrual duration, follow-up
duration, and relative accrual rates in a one-group negative binomial
design.

## Usage

``` r
nbsamplesize1s(
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
  lambdaH0 = NA_real_,
  accrualTime = 0L,
  accrualIntensity = NA_real_,
  piecewiseSurvivalTime = 0L,
  stratumFraction = 1L,
  kappa = NA_real_,
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

- lambdaH0:

  The rate parameter of the negative binomial distribution under the
  null hypothesis.

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

- kappa:

  The dispersion parameter (reciprocal of the shape parameter of the
  gamma mixing distribution) of the negative binomial distribution by
  stratum.

- lambda:

  The rate parameter of the negative binomial distribution under the
  alternative hypothesis by stratum.

- gamma:

  The hazard rate for exponential dropout or a vector of hazard rates
  for piecewise exponential dropout by stratum. Defaults to 0 for no
  dropout.

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

- `resultsUnderH1`: An S3 class `nbpower1s` object under the alternative
  hypothesis.

- `resultsUnderH0`: An S3 class `nbpower1s` object under the null
  hypothesis.

## See also

[`nbpower1s`](https://github.com/kaifenglu/lrstat/reference/nbpower1s.md)

## Author

Kaifeng Lu, <kaifenglu@gmail.com>

## Examples

``` r
# Example 1: Obtains follow-up duration given power, accrual intensity,
# and accrual duration for variable follow-up

nbsamplesize1s(beta = 0.2, kMax = 2,
               informationRates = c(0.5, 1),
               alpha = 0.025, typeAlphaSpending = "sfOF",
               lambdaH0 = 0.125, accrualIntensity = 500,
               stratumFraction = c(0.2, 0.8),
               kappa = c(3, 5), lambda = c(0.0875, 0.085),
               gamma = 0, accrualDuration = 1.25,
               followupTime = NA, fixedFollowup = FALSE)
#> $resultsUnderH1
#>                                                                             
#> Group-sequential design with 2 stages for one-sample negative binomial rate 
#> Rate under H0: 0.125, rate under H1: 0.0855                                 
#> Stratum fraction: 0.2 0.8, event rate: 0.0875 0.085, dispersion: 3 5        
#> Overall power: 0.8, overall significance level (1-sided): 0.025             
#> Maximum # events: 92.9, expected # events: 83.6                             
#> Maximum # dropouts: 0, expected # dropouts: 0                               
#> Maximum # subjects: 625, expected # subjects: 625                           
#> Maximum exposure: 1086.6, expected exposure: 978.1                          
#> Maximum information: 54.6, expected information: 50.12                      
#> Total study duration: 2.4, expected study duration: 2.2                     
#> Accrual duration: 1.2, follow-up duration: 1.1, fixed follow-up: FALSE      
#> Alpha spending: Lan-DeMets O'Brien-Fleming, beta spending: None             
#>                                                                             
#>                          Stage 1 Stage 2
#> Information rate         0.500   1.000  
#> Efficacy boundary (Z)    2.963   1.969  
#> Cumulative rejection     0.1641  0.8000 
#> Cumulative alpha spent   0.0015  0.0250 
#> Number of events         36.3    92.9   
#> Number of dropouts       0.0     0.0    
#> Number of subjects       625.0   625.0  
#> Exposure                 425.0   1086.6 
#> Analysis time            1.3     2.4    
#> Efficacy boundary (rate) 0.0709  0.0958 
#> Efficacy boundary (p)    0.0015  0.0245 
#> Information              27.30   54.60  
#> 
#> $resultsUnderH0
#>                                                                             
#> Group-sequential design with 2 stages for one-sample negative binomial rate 
#> Rate under H0: 0.125, rate under H1: 0.125                                  
#> Stratum fraction: 0.2 0.8, event rate: 0.125 0.125, dispersion: 3 5         
#> Overall power: 0.025, overall significance level (1-sided): 0.025           
#> Maximum # events: 94.7, expected # events: 94.6                             
#> Maximum # dropouts: 0, expected # dropouts: 0                               
#> Maximum # subjects: 625, expected # subjects: 624.9                         
#> Maximum exposure: 757.9, expected exposure: 757.2                           
#> Maximum information: 54.6, expected information: 54.55                      
#> Total study duration: 1.8, expected study duration: 1.8                     
#> Accrual duration: 1.2, follow-up duration: 0.6, fixed follow-up: FALSE      
#> Alpha spending: Lan-DeMets O'Brien-Fleming, beta spending: None             
#>                                                                             
#>                          Stage 1 Stage 2
#> Information rate         0.500   1.000  
#> Efficacy boundary (Z)    2.963   1.969  
#> Cumulative rejection     0.0015  0.0250 
#> Cumulative alpha spent   0.0015  0.0250 
#> Number of events         38.4    94.7   
#> Number of dropouts       0.0     0.0    
#> Number of subjects       554.4   625.0  
#> Exposure                 307.4   757.9  
#> Analysis time            1.1     1.8    
#> Efficacy boundary (rate) 0.0709  0.0958 
#> Efficacy boundary (p)    0.0015  0.0245 
#> Information              27.30   54.60  
#> 

# Example 2: Obtains accrual intensity given power, accrual duration, and
# follow-up duration for variable follow-up

nbsamplesize1s(beta = 0.2, kMax = 2,
               informationRates = c(0.5, 1),
               alpha = 0.025, typeAlphaSpending = "sfOF",
               lambdaH0 = 0.125, accrualIntensity = 100,
               kappa = 5, lambda = 0.0875,
               gamma = 0, accrualDuration = 1.25,
               followupTime = 2.25, fixedFollowup = FALSE)
#> $resultsUnderH1
#>                                                                             
#> Group-sequential design with 2 stages for one-sample negative binomial rate 
#> Rate under H0: 0.125, rate under H1: 0.0875                                 
#> Dispersion: 5                                                               
#> Overall power: 0.8, overall significance level (1-sided): 0.025             
#> Maximum # events: 140.4, expected # events: 124.7                           
#> Maximum # dropouts: 0, expected # dropouts: 0                               
#> Maximum # subjects: 558, expected # subjects: 558                           
#> Maximum exposure: 1604.2, expected exposure: 1425.1                         
#> Maximum information: 61.93, expected information: 56.85                     
#> Total study duration: 3.5, expected study duration: 3.2                     
#> Accrual duration: 1.2, follow-up duration: 2.2, fixed follow-up: FALSE      
#> Alpha spending: Lan-DeMets O'Brien-Fleming, beta spending: None             
#>                                                                             
#>                          Stage 1 Stage 2
#> Information rate         0.500   1.000  
#> Efficacy boundary (Z)    2.963   1.969  
#> Cumulative rejection     0.1641  0.8000 
#> Cumulative alpha spent   0.0015  0.0250 
#> Number of events         44.9    140.4  
#> Number of dropouts       0.0     0.0    
#> Number of subjects       558.0   558.0  
#> Exposure                 512.7   1604.2 
#> Analysis time            1.5     3.5    
#> Efficacy boundary (rate) 0.0734  0.0973 
#> Efficacy boundary (p)    0.0015  0.0245 
#> Information              30.96   61.93  
#> 
#> $resultsUnderH0
#>                                                                             
#> Group-sequential design with 2 stages for one-sample negative binomial rate 
#> Rate under H0: 0.125, rate under H1: 0.125                                  
#> Dispersion: 5                                                               
#> Overall power: 0.025, overall significance level (1-sided): 0.025           
#> Maximum # events: 141.7, expected # events: 141.5                           
#> Maximum # dropouts: 0, expected # dropouts: 0                               
#> Maximum # subjects: 558, expected # subjects: 558                           
#> Maximum exposure: 1133.2, expected exposure: 1132                           
#> Maximum information: 61.93, expected information: 61.88                     
#> Total study duration: 2.7, expected study duration: 2.7                     
#> Accrual duration: 1.2, follow-up duration: 1.4, fixed follow-up: FALSE      
#> Alpha spending: Lan-DeMets O'Brien-Fleming, beta spending: None             
#>                                                                             
#>                          Stage 1 Stage 2
#> Information rate         0.500   1.000  
#> Efficacy boundary (Z)    2.963   1.969  
#> Cumulative rejection     0.0015  0.0250 
#> Cumulative alpha spent   0.0015  0.0250 
#> Number of events         46.9    141.7  
#> Number of dropouts       0.0     0.0    
#> Number of subjects       558.0   558.0  
#> Exposure                 375.5   1133.2 
#> Analysis time            1.3     2.7    
#> Efficacy boundary (rate) 0.0734  0.0973 
#> Efficacy boundary (p)    0.0015  0.0245 
#> Information              30.96   61.93  
#> 


# Example 3: Obtains accrual duration given power, accrual intensity, and
# follow-up duration for fixed follow-up

nbsamplesize1s(beta = 0.2, kMax = 2,
               informationRates = c(0.5, 1),
               alpha = 0.025, typeAlphaSpending = "sfOF",
               lambdaH0 = 8.4, accrualIntensity = 40,
               kappa = 3, lambda = 4.2,
               gamma = 0, accrualDuration = NA,
               followupTime = 0.5, fixedFollowup = TRUE)
#> $resultsUnderH1
#>                                                                             
#> Group-sequential design with 2 stages for one-sample negative binomial rate 
#> Rate under H0: 8.4, rate under H1: 4.2                                      
#> Dispersion: 3                                                               
#> Overall power: 0.8, overall significance level (1-sided): 0.025             
#> Maximum # events: 112.7, expected # events: 102                             
#> Maximum # dropouts: 0, expected # dropouts: 0                               
#> Maximum # subjects: 58, expected # subjects: 53.8                           
#> Maximum exposure: 26.8, expected exposure: 24.3                             
#> Maximum information: 16.4, expected information: 15.05                      
#> Total study duration: 1.6, expected study duration: 1.5                     
#> Accrual duration: 1.4, follow-up duration: 0.5, fixed follow-up: TRUE       
#> Alpha spending: Lan-DeMets O'Brien-Fleming, beta spending: None             
#>                                                                             
#>                          Stage 1 Stage 2
#> Information rate         0.500   1.000  
#> Efficacy boundary (Z)    2.963   1.969  
#> Cumulative rejection     0.1641  0.8000 
#> Cumulative alpha spent   0.0015  0.0250 
#> Number of events         47.5    112.7  
#> Number of dropouts       0.0     0.0    
#> Number of subjects       32.6    58.0   
#> Exposure                 11.3    26.8   
#> Analysis time            0.8     1.6    
#> Efficacy boundary (rate) 2.9849  5.1659 
#> Efficacy boundary (p)    0.0015  0.0245 
#> Information              8.20    16.40  
#> 
#> $resultsUnderH0
#>                                                                             
#> Group-sequential design with 2 stages for one-sample negative binomial rate 
#> Rate under H0: 8.4, rate under H1: 8.4                                      
#> Dispersion: 3                                                               
#> Overall power: 0.025, overall significance level (1-sided): 0.025           
#> Maximum # events: 223, expected # events: 222.8                             
#> Maximum # dropouts: 0, expected # dropouts: 0                               
#> Maximum # subjects: 53.1, expected # subjects: 53.1                         
#> Maximum exposure: 26.5, expected exposure: 26.5                             
#> Maximum information: 16.4, expected information: 16.38                      
#> Total study duration: 1.8, expected study duration: 1.8                     
#> Accrual duration: 1.3, follow-up duration: 0.5, fixed follow-up: TRUE       
#> Alpha spending: Lan-DeMets O'Brien-Fleming, beta spending: None             
#>                                                                             
#>                          Stage 1 Stage 2
#> Information rate         0.500   1.000  
#> Efficacy boundary (Z)    2.963   1.969  
#> Cumulative rejection     0.0015  0.0250 
#> Cumulative alpha spent   0.0015  0.0250 
#> Number of events         81.6    223.0  
#> Number of dropouts       0.0     0.0    
#> Number of subjects       29.4    53.1   
#> Exposure                 9.7     26.5   
#> Analysis time            0.7     1.8    
#> Efficacy boundary (rate) 2.9849  5.1659 
#> Efficacy boundary (p)    0.0015  0.0245 
#> Information              8.20    16.40  
#> 
```
