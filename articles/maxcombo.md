# Power Calculation Using Max-Combo Tests

``` r
library(lrstat)
```

This R Markdown document illustrates the power calculation using the
maximum of weighted log-rank statistics in a group sequential analysis
for the delayed effect example from Prior (2020). The hazards in both
arms are 0.25 per month for the first 1.5 months, and the hazard of the
active arm drops to 0.125 per month thereafter, while the hazard of the
placebo arm remains at 0.25 per month. Assume that there is no
censoring. The enrollment rate is 25 patients per month, and the total
number of patients to enroll is 100. Suppose also that an interim
analysis occurs after observing 50 events and the final analysis takes
place after all 100 events have been observed. The conventional log-rank
test will be used for the interim analysis, and the maximum of the
standardized weighted log-rank test statistics of FH(0,0) and FH(0,1)
will be used at the final analysis. For a total one-sided 2.5%
significance level with an interim alpha of 0.0015, Prior (2020) reports
that the resulting power is abut 73%. We now verify this through
analytic calculation and simulation.

First, we derive the critical values at the interim and final analyses.
For the interim analysis, only one log-rank test will be used, hence the
critical value is \\u_1 = \Phi^{-1}(1-0.0015) = 2.968\\. For the final
analysis, the critical value \\u_2\\ satisfies \\ P_0(Z\_{1,1} \< u_1,
\max(Z\_{1,2}, Z\_{2,2}) \< u_2) = 1 - 0.025 \\ which is the same as \\
P_0(Z\_{1,1} \< u_1, Z\_{1,2} \< u_2, Z\_{2,2} \< u_2) = 0.975 \\
(Eq. 1) \\ Let \\U\_{FH(p,q),i}\\ denote the numerator of the FH(p,q)
weighted log-rank test statistic at analysis \\i\\, and
\\V\_{FH(p,q),i}\\ its variance. Then \\W\_{FH(p,q),i} =
U\_{FH(p,q),i}/{V\_{FH(p,q),i}^{1/2}}\\. In addition, similar to
Karrison (2016), we can show that \\ Cov(U\_{FH(p_1,q_1),1},
U\_{FH(p_2,q_2),2}) =
V\_{FH\left(\frac{p_1+p_2}{2},\frac{q_1+q_2}{2}\right),1} \\

First, we find the time when 50 and 99.9 events are expected to have
occurred.

``` r
(time = caltime(nevents=c(50, 99.9), accrualIntensity=25,
                piecewiseSurvivalTime=c(0, 1.5), 
                lambda1=c(0.25, 0.125), lambda2=c(0.25, 0.25), 
                accrualDuration=4, followupTime=60))
#> [1]  5.362939 50.323682
```

Then, we obtain the the means and variances of weighted log-rank test
statistics at the interim and final analyses for relevant FH weights.

``` r
(lr00 = lrstat(time=c(5.363, 50.324), accrualIntensity=25,
               piecewiseSurvivalTime=c(0, 1.5), 
               lambda1=c(0.25, 0.125), lambda2=c(0.25, 0.25), 
               accrualDuration=4, followupTime=60, 
               rho1=0, rho2=0))
#>     time subjects  nevents nevents1 nevents2 ndropouts ndropouts1 ndropouts2
#> 1  5.363      100 50.00056 22.47993 27.52063         0          0          0
#> 2 50.324      100 99.90000 49.90030 49.99970         0          0          0
#>   nfmax nfmax1 nfmax2     uscore   vscore  logRankZ hazardRatioH0
#> 1     0      0      0  -3.178997 12.46380 -0.900461             1
#> 2     0      0      0 -10.544685 22.26989 -2.234470             1

(lr01 = lrstat(time=c(5.363, 50.324), accrualIntensity=25,
               piecewiseSurvivalTime=c(0, 1.5), 
               lambda1=c(0.25, 0.125), lambda2=c(0.25, 0.25), 
               accrualDuration=4, followupTime=60, 
               rho1=0, rho2=1))
#>     time subjects  nevents nevents1 nevents2 ndropouts ndropouts1 ndropouts2
#> 1  5.363      100 50.00056 22.47993 27.52063         0          0          0
#> 2 50.324      100 99.90000 49.90030 49.99970         0          0          0
#>   nfmax nfmax1 nfmax2    uscore   vscore  logRankZ hazardRatioH0
#> 1     0      0      0 -1.385268 1.164020 -1.283966             1
#> 2     0      0      0 -6.608297 6.161009 -2.662341             1

(lr0h = lrstat(time=c(5.363, 50.324), accrualIntensity=25,
               piecewiseSurvivalTime=c(0, 1.5), 
               lambda1=c(0.25, 0.125), lambda2=c(0.25, 0.25), 
               accrualDuration=4, followupTime=60, 
               rho1=0, rho2=0.5))
#>     time subjects  nevents nevents1 nevents2 ndropouts ndropouts1 ndropouts2
#> 1  5.363      100 50.00056 22.47993 27.52063         0          0          0
#> 2 50.324      100 99.90000 49.90030 49.99970         0          0          0
#>   nfmax nfmax1 nfmax2    uscore    vscore  logRankZ hazardRatioH0
#> 1     0      0      0 -2.089089  3.241192 -1.160391             1
#> 2     0      0      0 -8.253880 10.082599 -2.599393             1
```

It follows that the mean of \\\mathbf{Z}=(Z\_{1,1}, Z\_{1,2},
Z\_{2,2})\\ is \\ \left(\frac{3.179}{\sqrt{12.464}},
\frac{10.545}{\sqrt{22.270}}, \frac{6.608}{\sqrt{6.161}} \right) =
(0.900, 2.234, 2.662) \\ and the covariance matrix of \\\mathbf{Z}\\ is
\\ \left(\begin{array}{ccc} 1 & \sqrt{\frac{12.464}{22.270}} &
\frac{3.241}{\sqrt{12.464\times 6.161}} \\ \sqrt{\frac{12.464}{22.270}}
& 1 & \frac{10.083}{\sqrt{22.270\times 6.161}} \\
\frac{3.241}{\sqrt{12.464\times 6.161}} &
\frac{10.083}{\sqrt{22.270\times 6.161}} & 1 \end{array} \right) =
\left(\begin{array}{ccc} 1 & 0.748 & 0.370 \\ 0.748 & 1 & 0.861 \\ 0.370
& 0.861 & 1 \end{array} \right) \\ Now, we obtain the critical value
\\u_2\\ by solving equation (1).

``` r
mu = c(0.900, 2.234, 2.662)
sigma = matrix(c(1, 0.748, 0.370, 0.748, 1, 0.861, 0.370, 0.861, 1), 3, 3)
u1 = 2.968
alpha = 0.025
f <- function(u2, u1, sigma, alpha) {
  1 - pmvnormr(upper=c(u1, u2, u2), sigma=sigma) - alpha
}
(u2 = uniroot(f, c(1,3), u1, sigma, alpha)$root)
#> [1] 2.136503
```

The power can be estimated by plugging in the mean under the alternative
hypothesis.

``` r
1 - pmvnormr(upper=c(u1, u2, u2), mean=mu, sigma=sigma)
#> [1] 0.7244549
#> attr(,"method")
#> [1] "qmc"
#> attr(,"error")
#> [1] 8.932615e-05
#> attr(,"nsamples")
#> [1] 24576
```

For the simulation study, we use very large critical values for the
FH(0,0) and FH(0,1) log-rank statistics at the interim and final
analyses for each iteration, and then construct the maximum combo test
statistic at the final analysis. Finally, we tally the rejections across
iterations to estimate the power. The same seed should be used to
produce identical simulated data.

``` r
sim1 = lrsim(kMax = 2, informationRates = c(0.5, 1),
             criticalValues = c(6, 6),
             accrualIntensity = 25, 
             piecewiseSurvivalTime = c(0, 1.5),
             lambda1 = c(0.25, 0.125), lambda2 = c(0.25, 0.25),
             n = 100,
             rho1 = 0, rho2 = 0,
             plannedEvents = c(50, 100), 
             maxNumberOfIterations = 10000,
             seed = 314159)

sim2 = lrsim(kMax = 2, informationRates = c(0.5, 1),
             criticalValues = c(6, 6),
             accrualIntensity = 25,
             piecewiseSurvivalTime = c(0, 1.5),
             lambda1 = c(0.25, 0.125), lambda2 = c(0.25, 0.25),
             n = 100,
             rho1 = 0, rho2 = 1, 
             plannedEvents = c(50, 100),
             maxNumberOfIterations = 10000,
             seed = 314159)

w1max = subset(-sim1$sumdata$logRankStatistic, sim1$sumdata$stageNumber==1)

w2max = pmax(-sim1$sumdata$logRankStatistic, -sim2$sumdata$logRankStatistic) 
w2max = subset(w2max, sim1$sumdata$stageNumber==2)

mean((w1max > u1) | (w2max > u2))
#> [1] 0.7256
```

The analytic method yielded a power of about 72%, while the simulation
method yielded a power of about 73%, both are close to that reported by
Prior (2020).

Karrison, Theodore G. 2016. “Versatile Tests for Comparing Survival
Curves Based on Weighted Log-Rank Statistics.” *The Stata Journal:
Promoting Communications on Statistics and Stata* 16 (3): 678–90.

Prior, Thomas J. 2020. “Group Sequential Monitoring Based on the Maximum
of Weighted Log-Rank Statistics with the FlemingHarrington Class of
Weights in Oncology Clinical Trials.” *Statistical Methods in Medical
Research* 29 (12): 3525–32.
