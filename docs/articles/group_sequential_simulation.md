# Simulation for Group Sequential Trials

``` r
library(lrstat)
```

This R Markdown document illustrates the simulation tool for group
sequential survival trials. This is useful for validating the analytic
calculation, which might be inaccurate when the allocation ratio is not
1:1 or the hazard ratio is far from 1.

Consider a three-stage O’Brien-Fleming group sequential design with two
interim looks conducted at 50% and 75% of the target total number of
events. The first interim is for futility only, and the second interim
is for efficacy only. The hazard rate of the control group is 0.95 per
year. The hazard ratio of the experimental group to the control group is
0.3. The experimental versus control group randomization ratio is 3:1.
The enrollment rate is 5 patients per month. The 2-year drop-out rate is
10%. The study uses a fixed follow-up design and each patient is to be
followed up for 26 weeks. If we use an enrollment duration of 32 months,
then the maximum number of events is expected to be 32.8.

``` r
lrstat(time=c(20, 25, 30, 35, 38.5), allocationRatioPlanned = 3, 
       accrualIntensity = 5, 
       lambda2 = 0.95/12, lambda1 = 0.3*0.95/12, 
       gamma1 = -log(1-0.1)/24, gamma2 = -log(1-0.1)/24, 
       accrualDuration = 32, followupTime = 6.5, fixedFollowup = TRUE)
#>   time subjects  nevents nevents1  nevents2 ndropouts ndropouts1 ndropouts2
#> 1 20.0      100 17.37244  8.91394  8.458504  2.116728   1.647680  0.4690486
#> 2 25.0      125 22.49924 11.55918 10.940058  2.743292   2.136634  0.6066580
#> 3 30.0      150 27.62603 14.20442 13.421612  3.369855   2.625588  0.7442674
#> 4 35.0      160 31.95277 16.45992 15.492843  3.901625   3.042502  0.8591231
#> 5 38.5      160 32.81148 16.92953 15.881948  4.010006   3.129306  0.8807001
#>       nfmax   nfmax1    nfmax2    uscore   vscore  logRankZ hazardRatioH0
#> 1  51.96593 42.16268  9.803258 -4.600670 2.995801 -2.658059             1
#> 2  71.21258 57.77848 13.434094 -5.953764 3.874267 -3.024801             1
#> 3  90.45922 73.39429 17.064930 -7.306858 4.752732 -3.351652             1
#> 4 109.70586 89.01009 20.695767 -8.441729 5.485262 -3.604398             1
#> 5 123.17851 99.94116 23.237352 -8.659802 5.622182 -3.652208             1
```

Suppose we run the trial for a target maximum number of 32 events. The
trial will stop once 32 events have been observed or the trial is
stopped early for futility or efficacy. Due to the fixed follow-up
design, there might be situations where a total of 160 patients with
each followed-up for 6.5 month do not yield 32 events, in which case,
the trial will stop and we allocate all remaining alpha to the final
look. The simulation below demonstrates that the study has approximately
90% under the alternative hypothesis.

``` r
lrsim(kMax = 3, informationRates = c(0.5, 0.75, 1), 
      criticalValues = c(6, 2.34, 2.012), 
      futilityBounds = c(0.282, -6), 
      allocation1 = 3, allocation2 = 1,
      accrualTime = 0, accrualIntensity = 5, 
      piecewiseSurvivalTime = 0, 
      stratumFraction = 1, 
      lambda1 = 0.3*0.95/12, lambda2 = 0.95/12, 
      gamma1 = -log(1-0.1)/24, gamma2 = -log(1-0.1)/24, 
      n = 160, followupTime = 6.5, 
      fixedFollowup = TRUE, 
      rho1 = 0, rho2 = 0, 
      plannedEvents = c(16, 24, 32), 
      maxNumberOfIterations = 1000, 
      maxNumberOfRawDatasetsPerStage = 0, 
      seed = 12345)
#>                                                         
#> Group-sequential design with 3 stages for log-rank test 
#> Empirical power: 0.901                                  
#> Expected # events: 25.4                                 
#> Expected # dropouts: 3                                  
#> Expected # subjects: 136.1                              
#> Expected study duration: 27.7                           
#> n: 160, fixed follow-up: TRUE                           
#> Number of simulations: 1000                             
#>                                                         
#>                      Stage 1 Stage 2 Stage 3
#> Cumulative rejection 0.0020  0.7520  0.9010 
#> Cumulative futility  0.0180  0.0220  0.0990 
#> Number of events     16.0    24.0    31.0   
#> Number of dropouts   1.9     2.8     3.5    
#> Number of subjects   93.3    130.9   154.7  
#> Analysis time        18.7    26.3    32.6
```

The simulation below shows that the probability of futility stopping at
the first look is about 62% under the null hypothesis.

``` r
lrsim(kMax = 3, informationRates = c(0.5, 0.75, 1), 
      criticalValues = c(6, 2.34, 2.012), 
      futilityBounds = c(0.282, -6), 
      allocation1 = 3, allocation2 = 1,
      accrualTime = 0, accrualIntensity = 5, 
      piecewiseSurvivalTime = 0, 
      stratumFraction = 1, 
      lambda1 = 0.95/12, lambda2 = 0.95/12, 
      gamma1 = -log(1-0.1)/24, gamma2 = -log(1-0.1)/24, 
      n = 160, followupTime = 6.5, 
      fixedFollowup = TRUE, 
      rho1 = 0, rho2 = 0, 
      plannedEvents = c(16, 24, 32), 
      maxNumberOfIterations = 1000, 
      maxNumberOfRawDatasetsPerStage = 0, 
      seed = 12345)
#>                                                         
#> Group-sequential design with 3 stages for log-rank test 
#> Empirical power: 0.038                                  
#> Expected # events: 22                                   
#> Expected # dropouts: 1.2                                
#> Expected # subjects: 70.3                               
#> Expected study duration: 14                             
#> n: 160, fixed follow-up: TRUE                           
#> Number of simulations: 1000                             
#>                                                         
#>                      Stage 1 Stage 2 Stage 3
#> Cumulative rejection 0.0000  0.0170  0.0380 
#> Cumulative futility  0.6170  0.6170  0.9620 
#> Number of events     16.0    24.0    32.0   
#> Number of dropouts   0.9     1.3     1.7    
#> Number of subjects   55.2    75.1    95.2   
#> Analysis time        11.0    14.9    19.0
```
