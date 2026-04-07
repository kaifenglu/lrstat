# Sample Size Calculation Under Non-Proportional Hazards

``` r
library(lrstat)
```

This R Markdown document illustrates the sample size calculation for a
delayed effect model using the lrsamplesize function from lrstat and
verifies the result using the simulation tool.

Suppose that the survival distribution of the control group is
exponential with a median survival time of 13 months. The survival
distribution of the active treatment group is piecewise exponential with
the same hazard rate as the control group for the first 6 months and
with a hazard ratio of 0.58 afterwards. The accrual has a ramp-up period
of 9 months to reach 26 patients per months thereafter, and the total
duration of the enrollment is 22 months. In addition, the annual drop
rate is 5% for each treatment group.

We would like to know the total number of events needed to achieve 80%
power for a two-stage group sequential trial with O’Brien-Fleming
spending function and with the interim analysis to be conducted after
observing 80% of the target total number of events.

We use the lrsamplesize function to obtain the follow-up time and the
target number of events.

``` r
lrsamplesize(beta = 0.2, kMax = 2, 
             informationRates = c(0.8, 1),
             alpha = 0.025, typeAlphaSpending = "sfOF", 
             accrualTime = seq(0, 8),
             accrualIntensity = 26/9*seq(1, 9),
             piecewiseSurvivalTime = c(0, 6),
             lambda2 = rep(log(2)/13, 2),
             lambda1 = c(log(2)/13, 0.58*log(2)/13),
             gamma1 = -log(1-0.05)/12, 
             gamma2 = -log(1-0.05)/12,
             accrualDuration = 22, followupTime = NA)$resultsUnderH1
#>                                                                        
#> Group-sequential design with 2 stages for log-rank test                
#> Overall power: 0.8016, overall significance level (1-sided): 0.025     
#> Maximum # events: 315, expected # events: 286.6                        
#> Maximum # dropouts: 29.8, expected # dropouts: 26.8                    
#> Maximum # subjects: 468, expected # subjects: 468                      
#> Maximum information: 77.96, expected information: 71.1                 
#> Total study duration: 41.5, expected study duration: 36.9              
#> Accrual duration: 22, follow-up duration: 19.5, fixed follow-up: FALSE 
#> Allocation ratio: 1                                                    
#> Alpha spending: Lan-DeMets O'Brien-Fleming, beta spending: None        
#>                                                                        
#>                        Stage 1 Stage 2
#> Information rate       0.800   1.000  
#> Efficacy boundary (Z)  2.250   2.025  
#> Cumulative rejection   0.4511  0.8016 
#> Cumulative alpha spent 0.0122  0.0250 
#> Number of events       252.0   315.0  
#> Number of dropouts     23.1    29.8   
#> Number of subjects     468.0   468.0  
#> Analysis time          31.4    41.5   
#> Efficacy boundary (HR) 0.752   0.795  
#> Efficacy boundary (p)  0.0122  0.0214 
#> Information            62.76   77.96  
#> HR                     0.764   0.723
```

Thus we need to observe 315 events with 468 subjects, and the maximum
study duration is 41.5 months with an expected study duration of 37
months.

To verify this requirement, we resort to the lrsim function.

``` r
lrsim(kMax = 2, criticalValues = c(2.250, 2.025), 
      accrualTime = seq(0, 9),
      accrualIntensity = c(26/9*seq(1, 9), 26),
      piecewiseSurvivalTime = c(0, 6),
      lambda2 = rep(log(2)/13, 2),
      lambda1 = c(log(2)/13, 0.58*log(2)/13),
      gamma1 = -log(1-0.05)/12, 
      gamma2 = -log(1-0.05)/12,
      n = 468,
      plannedEvents = c(252, 315), 
      maxNumberOfIterations = 10000, seed = 314159)
#>                                                         
#> Group-sequential design with 2 stages for log-rank test 
#> Empirical power: 0.8016                                 
#> Expected # events: 286.3                                
#> Expected # dropouts: 26.7                               
#> Expected # subjects: 468                                
#> Expected study duration: 36.8                           
#> n: 468, fixed follow-up: FALSE                          
#> Number of simulations: 10000                            
#>                                                         
#>                      Stage 1 Stage 2
#> Cumulative rejection 0.4553  0.8016 
#> Cumulative futility  0.0000  0.1984 
#> Number of events     252.0   315.0  
#> Number of dropouts   23.1    29.6   
#> Number of subjects   468.0   468.0  
#> Analysis time        31.4    41.2
```

The simulation results confirm the analytic calculations.
