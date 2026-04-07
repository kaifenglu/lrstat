# Sample Size Calculation With Fixed Follow-up

``` r
library(lrstat)
```

This R Markdown document illustrates the sample size calculation for a
fixed follow-up design, in which the treatment allocation is 3:1 and the
hazard ratio is 0.3. This is a case for which neither the Schoenfeld
method nor the Lakatos method provides an accurate sample size estimate,
and simulation tools are needed to obtain a more accurate result.

Consider a fixed design with the hazard rate of the control group being
0.95 per year, a hazard ratio of the experimental group to the control
group being 0.3, a randomization ratio of 3:1, an enrollment rate of 5
patients per month, a 2-year drop-out rate of 10%, and a planned fixed
follow-up of 26 weeks for each patient. The target power is 90%, and we
are interested in the number of patients to enroll to achieve the target
90% power.

Using the Schoenfeld formula, the required number of events is 39. This
requires 191 patients enrolled over 38.2 months. Denote this design as
design 1.

``` r
lrsamplesize(beta = 0.1, kMax = 1, criticalValues = 1.96, 
             allocationRatioPlanned = 3, accrualIntensity = 5, 
             lambda2 = 0.95/12, lambda1 = 0.3*0.95/12, 
             gamma1 = -log(1-0.1)/24, gamma2 = -log(1-0.1)/24, 
             accrualDuration = NA, followupTime = 26/4, 
             fixedFollowup = TRUE,
             typeOfComputation = "schoenfeld")
#> $resultsUnderH1
#>                                                                        
#> Fixed design for log-rank test                                         
#> Overall power: 0.9025, overall significance level (1-sided): 0.025     
#> Number of events: 39                                                   
#> Number of dropouts: 4.8                                                
#> Number of subjects: 191                                                
#> Information: 7.31                                                      
#> Study duration: 43.1                                                   
#> Accrual duration: 38.2, follow-up duration: 6.5, fixed follow-up: TRUE 
#> Allocation ratio: 3                                                    
#>                                                                        
#>                              
#> Efficacy boundary (Z)  1.960 
#> Efficacy boundary (HR) 0.484 
#> Efficacy boundary (p)  0.0250
#> HR                     0.300 
#> 
#> $resultsUnderH0
#>                                                                        
#> Fixed design for log-rank test                                         
#> Overall power: 0.025, overall significance level (1-sided): 0.025      
#> Number of events: 39                                                   
#> Number of dropouts: 2.2                                                
#> Number of subjects: 98.2                                               
#> Information: 7.31                                                      
#> Study duration: 26.1                                                   
#> Accrual duration: 19.6, follow-up duration: 6.5, fixed follow-up: TRUE 
#> Allocation ratio: 3                                                    
#>                                                                        
#>                              
#> Efficacy boundary (Z)  1.960 
#> Efficacy boundary (HR) 0.484 
#> Efficacy boundary (p)  0.0250
#> HR                     1.000
```

On the other hand, the output from the lrsamplesize call using the
‘direct’ method implies that we only need 26 events with 127 subjects
enrolled over 25.4 months, a dramatic difference from the Schoenfeld
formula. Denote this design as design 2.

``` r
lrsamplesize(beta = 0.1, kMax = 1, criticalValues = 1.96, 
             allocationRatioPlanned = 3, accrualIntensity = 5, 
             lambda2 = 0.95/12, lambda1 = 0.3*0.95/12, 
             gamma1 = -log(1-0.1)/24, gamma2 = -log(1-0.1)/24, 
             accrualDuration = NA, followupTime = 26/4, 
             fixedFollowup = TRUE,
             typeOfComputation = "direct")
#> $resultsUnderH1
#>                                                                        
#> Fixed design for log-rank test                                         
#> Overall power: 0.9017, overall significance level (1-sided): 0.025     
#> Number of events: 26                                                   
#> Number of dropouts: 3.2                                                
#> Number of subjects: 127                                                
#> Information: 4.46                                                      
#> Study duration: 31.1                                                   
#> Accrual duration: 25.4, follow-up duration: 6.5, fixed follow-up: TRUE 
#> Allocation ratio: 3                                                    
#>                                                                        
#>                              
#> Efficacy boundary (Z)  1.960 
#> Efficacy boundary (HR) 0.463 
#> Efficacy boundary (p)  0.0250
#> HR                     0.300 
#> 
#> $resultsUnderH0
#>                                                                        
#> Fixed design for log-rank test                                         
#> Overall power: 0.025, overall significance level (1-sided): 0.025      
#> Number of events: 26                                                   
#> Number of dropouts: 1.4                                                
#> Number of subjects: 65.5                                               
#> Information: 4.88                                                      
#> Study duration: 19.6                                                   
#> Accrual duration: 13.1, follow-up duration: 6.5, fixed follow-up: TRUE 
#> Allocation ratio: 3                                                    
#>                                                                        
#>                              
#> Efficacy boundary (Z)  1.960 
#> Efficacy boundary (HR) 0.412 
#> Efficacy boundary (p)  0.0250
#> HR                     1.000
```

To check the accuracy of either solution, we run simulations using the
lrsim function.

``` r
lrsim(kMax = 1, criticalValues = 1.96,  
      allocation1 = 3, allocation2 = 1,
      accrualIntensity = 5, 
      lambda2 = 0.95/12, lambda1 = 0.3*0.95/12, 
      gamma1 = -log(1-0.1)/24, gamma2 = -log(1-0.1)/24,
      n = 191, followupTime = 6.5, 
      fixedFollowup = TRUE,  
      plannedEvents = 39, 
      maxNumberOfIterations = 10000, seed = 12345)
#>                                
#> Fixed design for log-rank test 
#> Empirical power: 0.9502        
#> Expected # events: 37          
#> Expected # dropouts: 4.4       
#> Expected # subjects: 186.2     
#> Expected study duration: 39.2  
#> n: 191, fixed follow-up: TRUE  
#> Number of simulations: 10000   
#> 

lrsim(kMax = 1, criticalValues = 1.96,  
      allocation1 = 3, allocation2 = 1,
      accrualIntensity = 5, 
      lambda2 = 0.95/12, lambda1 = 0.3*0.95/12, 
      gamma1 = -log(1-0.1)/24, gamma2 = -log(1-0.1)/24,
      n = 127, followupTime = 6.5, 
      fixedFollowup = TRUE,  
      plannedEvents = 26, 
      maxNumberOfIterations = 10000, seed = 12345)
#>                                
#> Fixed design for log-rank test 
#> Empirical power: 0.8361        
#> Expected # events: 24.3        
#> Expected # dropouts: 2.9       
#> Expected # subjects: 124       
#> Expected study duration: 26.9  
#> n: 127, fixed follow-up: TRUE  
#> Number of simulations: 10000   
#> 
```

The simulated power is about 95% for design 1, and 84% for design 2.
Neither is close to the target 90% power.

We use the following formula to adjust the sample size to attain the
target power,
``` math

D = D_0 \left( \frac{\Phi^{-1}(1-\alpha) + \Phi^{-1}(1-\beta)} {\Phi^{-1}(1-\alpha) + \Phi^{-1}(1-\beta_0)} \right)^2
```
where $`D_0`$ and $`\beta_0`$ are the initial event number and the
correponding type II error, and $`D`$ and $`\beta`$ are the required
event number and the target type II error, respectively. For
$`\alpha=0.025`$ and $`\beta=0.1`$, plugging in
$`(D_0=39, \beta_0=0.05)`$ and $`(D_0=26, \beta_0=0.17)`$ would yield
$`D=32`$ and $`D=32`$, respectively. For $`D=32`$, we need about 156
patients for an enrollment period of 31.2 months,  
``` math

N = \frac{D}{ \frac{r}{1+r}\frac{\lambda_1}{\lambda_1+\gamma_1} (1 - \exp(-(\lambda_1+\gamma_1)T_f)) + \frac{1}{1+r}\frac{\lambda_2}{\lambda_2+\gamma_2} (1 - \exp(-(\lambda_2+\gamma_2)T_f)) }
```
Simulation results confirmed the accuracy of this sample size estimate.

``` r
lrsim(kMax = 1, criticalValues = 1.96,  
      allocation1 = 3, allocation2 = 1,
      accrualIntensity = 5, 
      lambda2 = 0.95/12, lambda1 = 0.3*0.95/12, 
      gamma1 = -log(1-0.1)/24, gamma2 = -log(1-0.1)/24,
      n = 156, followupTime = 6.5, 
      fixedFollowup = TRUE,  
      plannedEvents = 32, 
      maxNumberOfIterations = 10000, seed = 12345)
#>                                
#> Fixed design for log-rank test 
#> Empirical power: 0.8984        
#> Expected # events: 30.1        
#> Expected # dropouts: 3.6       
#> Expected # subjects: 152.3     
#> Expected study duration: 32.5  
#> n: 156, fixed follow-up: TRUE  
#> Number of simulations: 10000   
#> 
```
