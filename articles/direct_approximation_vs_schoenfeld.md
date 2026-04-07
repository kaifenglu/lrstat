# Comparing Direct Approximation and Schoenfeld Methods

``` r
library(lrstat)
```

This R Markdown document compares the direct approximation method for
sample size calculation of survival trials as implemented in the lrstat
package with the Schoenfeld method for sample size calculation of
survival trials.

Consider a fixed design with \\\pi_1=0.2\\ for the active treatment
group and \\\pi_2=0.4\\ for the control group at month 12, an enrollment
period lasting 12 months to enroll 200 patients, and a target number of
events of 40.

First, we find the follow-up time at which 40 events will be observed.
This can be accomplished by the caltime function using a large upper
bound of the follow-up time for bracketed search:

``` r
caltime(nevents = 40, accrualDuration = 12, accrualIntensity = 200/12,
        lambda1 = -log(1-0.2)/12, lambda2 = -log(1-0.4)/12, 
        followupTime = 100)
#> [1] 13.63169
```

Therefore the followup time for the last enrolled patient is 13.63 - 12
= 1.63 months. Next we use the lrpower function to estimate the power.

``` r
lrpower(kMax = 1, criticalValues = 1.96, accrualDuration = 12, 
        accrualIntensity = 200/12, lambda1 = -log(1-0.2)/12, 
        lambda2 = -log(1-0.4)/12,  followupTime = 1.63, 
        typeOfComputation = "direct")
#>                                                                       
#> Fixed design for log-rank test                                        
#> Overall power: 0.7138, overall significance level (1-sided): 0.025    
#> Number of events: 40                                                  
#> Number of dropouts: 0                                                 
#> Number of subjects: 200                                               
#> Information: 9.96                                                     
#> Study duration: 13.6                                                  
#> Accrual duration: 12, follow-up duration: 1.6, fixed follow-up: FALSE 
#> Allocation ratio: 1                                                   
#>                                                                       
#>                              
#> Efficacy boundary (Z)  1.960 
#> Efficacy boundary (HR) 0.516 
#> Efficacy boundary (p)  0.0250
#> HR                     0.437
```

Therefore, the power is estimated to be 71.4%. In comparison, the
Schoenfeld formula yields a power of 74.5%.

``` r
hazardRatio = log(1-0.2)/log(1-0.4)
pnorm(abs(log(hazardRatio))*sqrt(40/4) - qnorm(0.975))
#> [1] 0.7450763
```

To see which method yields a more accurate power estimate, we run a
lrsim function, which reports a power of 72.2% with 10000 replications.

``` r
lrsim(kMax = 1, criticalValues = 1.96, 
      accrualIntensity = 200/12, 
      lambda1 = -log(1-0.2)/12, lambda2 = -log(1-0.4)/12,
      n = 200, 
      plannedEvents = 40, 
      maxNumberOfIterations = 10000, seed = 314159)
#>                                
#> Fixed design for log-rank test 
#> Empirical power: 0.7224        
#> Expected # events: 40          
#> Expected # dropouts: 0         
#> Expected # subjects: 199.3     
#> Expected study duration: 13.7  
#> n: 200, fixed follow-up: FALSE 
#> Number of simulations: 10000   
#> 
```

Therefore, the direct approximation method yields a more accurate power
estimate than the Schoenfeld method in this case. The reason for this
discrepancy is due to the hazard ratio being 0.437, which is rather far
away from the null hazard ratio of 1. It is well-known that the
Schoenfeld method yields an optimistic power estimate with equal
randomization in this case.
