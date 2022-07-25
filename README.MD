
# lrstat

<!-- badges: start -->
<!-- badges: end -->

The goal of lrstat is to calculate power and sample size under non-proportional
hazards using weighted log-rank tests both analytically and using simulations
in a time-to-event group-sequential trial.

It has built-in capabilities to use error-spending functions and can calcualte
power, accrual duration, follow-up duration, and absolute accrual rates for 
the Fleming-Harrington's class of weighting functions. 

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(lrstat)
## basic example code for a two-stage group sequential trial with interim
## analysis at 80% of total number of events using Lan-DeMets O'Brien-Fleming
## error-spending. The accrual has a ramp-up periof of 9 months before 
## reaching 26 patients per month. The survival distribution for the treatment
## group has a delay effect of 6 months and a hazard ratio 0.58 after the delay. 
## The annual dropout rate is 5%. The accrual duration is 22 months. 
## The follow-up duration is 18 months for the last randomized patients. 
## The FH(0,1) weighted log-rank test is used for power calculation. 

lrpower(kMax = 2, informationRates = c(0.8, 1),
        alpha = 0.025, typeAlphaSpending = "sfOF",
        allocationRatioPlanned = 1, accrualTime = seq(0, 9),
        accrualIntensity = c(26/9*seq(1, 9), 26),
        piecewiseSurvivalTime = c(0, 6),
        lambda1 = c(0.0533, 0.0309),
        lambda2 = c(0.0533, 0.0533),
        gamma1 = -log(1-0.05)/12,
        gamma2 = -log(1-0.05)/12, 
        accrualDuration = 22,
        followupTime = 18, fixedFollowup = FALSE, 
        rho1 = 0, rho2 = 1,
        numSubintervals = 2000)

```

