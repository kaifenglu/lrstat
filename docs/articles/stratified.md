# Power Calculation With Stratification Variables

``` r
library(lrstat)
```

This R Markdown document illustrates the power calculation in the
presence of stratification variables. This example is taken from EAST
6.4 section 56.7 on lung cancer patients comparing two treatment groups
in a target patient population with some prior therapy. There are three
stratification variables:

- type of cancer cell (small, adeno, large, squamous)

- age in years (\<=50, \>50)

- performance status score (\<=50, \>50-\<=70, \>70)

We consider a three stage Lan-DeMets O’Brien-Fleming group sequential
design. The stratum fractions are

``` r
p1 = c(0.28, 0.13, 0.25, 0.34)
p2 = c(0.28, 0.72)
p3 = c(0.43, 0.37, 0.2)
stratumFraction = p1 %x% p2 %x% p3
stratumFraction = stratumFraction/sum(stratumFraction)
```

Using the small cancer cell, age \<=50, and performance status score
\<=50 as the reference stratum, the hazard ratios are

``` r
theta1 = c(1, 2.127, 0.528, 0.413)
theta2 = c(1, 0.438)
theta3 = c(1, 0.614, 0.159)
```

If the hazard rate of the reference stratum is 0.009211, then the hazard
rate for the control group is

``` r
lambda2 = 0.009211*exp(log(theta1) %x% log(theta2) %x% log(theta3))
```

The hazard ratio of the active treatment group versus the control group
is 0.4466.

In addition, we assume an enrollment period of 24 months with a constant
enrollment rate of 12 patients per month to enroll 288 patients, and the
target number of events of 66.

First we obtain the calendar time at which 66 events will occur.

``` r
caltime(nevents = 66, accrualDuration = 24, accrualIntensity = 12,
        stratumFraction = stratumFraction, 
        lambda1 = 0.4466*lambda2, lambda2 = lambda2, 
        followupTime = 100)
#> [1] 54.92196
```

Therefore, the follow-up time for the last enrolled patient is 30.92
months. Now we can evaluate the power using the lrpower function.

``` r
lrpower(kMax = 3, 
        informationRates = c(0.333, 0.667, 1), 
        alpha = 0.025, typeAlphaSpending = "sfOF", 
        accrualIntensity = 12,
        stratumFraction = stratumFraction,
        lambda1 = 0.4466*lambda2, 
        lambda2 = lambda2, 
        accrualDuration = 24, 
        followupTime = 30.92, 
        typeOfComputation = "schoenfeld")
#>                                                                        
#> Group-sequential design with 3 stages for log-rank test                
#> Overall power: 0.9024, overall significance level (1-sided): 0.025     
#> Maximum # events: 66, expected # events: 52.8                          
#> Maximum # dropouts: 0, expected # dropouts: 0                          
#> Maximum # subjects: 288, expected # subjects: 288                      
#> Maximum information: 16.5, expected information: 13.21                 
#> Total study duration: 54.9, expected study duration: 45.5              
#> Accrual duration: 24, follow-up duration: 30.9, fixed follow-up: FALSE 
#> Allocation ratio: 1                                                    
#> Alpha spending: Lan-DeMets O'Brien-Fleming, beta spending: None        
#>                                                                        
#>                        Stage 1 Stage 2 Stage 3
#> Information rate       0.333   0.667   1.000  
#> Efficacy boundary (Z)  3.712   2.511   1.993  
#> Cumulative rejection   0.0342  0.5652  0.9024 
#> Cumulative alpha spent 0.0001  0.0061  0.0250 
#> Number of events       22.0    44.0    66.0   
#> Number of dropouts     0.0     0.0     0.0    
#> Number of subjects     288.0   288.0   288.0  
#> Analysis time          24.9    39.0    54.9   
#> Efficacy boundary (HR) 0.205   0.469   0.612  
#> Efficacy boundary (p)  0.0001  0.0060  0.0231 
#> Information            5.49    11.01   16.50  
#> HR                     0.447   0.447   0.447
```

Therefore, the overall power is about 90% for the stratified analysis.
This is confirmed by the simulation below.

``` r
lrsim(kMax = 3, 
      informationRates = c(0.333, 0.667, 1), 
      criticalValues = c(3.712, 2.511, 1.993), 
      accrualIntensity = 12,
      stratumFraction = stratumFraction,
      lambda1 = 0.4466*lambda2, 
      lambda2 = lambda2, 
      n = 288, 
      followupTime = 30.92,
      plannedEvents = c(22, 44, 66),
      maxNumberOfIterations = 1000, 
      seed = 314159)
#>                                                         
#> Group-sequential design with 3 stages for log-rank test 
#> Empirical power: 0.889                                  
#> Expected # events: 54.6                                 
#> Expected # dropouts: 0                                  
#> Expected # subjects: 287.8                              
#> Expected study duration: 46.8                           
#> n: 288, fixed follow-up: FALSE                          
#> Number of simulations: 1000                             
#>                                                         
#>                      Stage 1 Stage 2 Stage 3
#> Cumulative rejection 0.0140  0.5060  0.8890 
#> Cumulative futility  0.0000  0.0000  0.1110 
#> Number of events     22.0    44.0    66.0   
#> Number of dropouts   0.0     0.0     0.0    
#> Number of subjects   279.8   288.0   288.0  
#> Analysis time        25.0    39.1    54.8
```
