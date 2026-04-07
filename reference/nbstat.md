# Negative Binomial Rate Ratio

Obtains the number of subjects accrued, number of events, number of
dropouts, number of subjects reaching the maximum follow-up, total
exposure, and variance for log rate in each group, rate ratio, variance,
and Wald test statistic of log rate ratio at given calendar times.

## Usage

``` r
nbstat(
  time = NA_real_,
  rateRatioH0 = 1,
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
  fixedFollowup = FALSE,
  nullVariance = FALSE
)
```

## Arguments

- time:

  A vector of calendar times for data cut.

- rateRatioH0:

  Rate ratio under the null hypothesis.

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

- nullVariance:

  Whether to calculate the variance for log rate ratio under the null
  hypothesis.

## Value

A list with two components:

- `resultsUnderH1`: A data frame containing the following variables:

  - `time`: The analysis time since trial start.

  - `subjects`: The number of enrolled subjects.

  - `nevents`: The total number of events.

  - `nevents1`: The number of events in the active treatment group.

  - `nevents2`: The number of events in the control group.

  - `ndropouts`: The total number of dropouts.

  - `ndropouts1`: The number of dropouts in the active treatment group.

  - `ndropouts2`: The number of dropouts in the control group.

  - `nfmax`: The total number of subjects reaching maximum follow-up.

  - `nfmax1`: The number of subjects reaching maximum follow-up in the
    active treatment group.

  - `nfmax2`: The number of subjects reaching maximum follow-up in the
    control group.

  - `exposure`: The total exposure time.

  - `exposure1`: The exposure time for the active treatment group.

  - `exposure2`: The exposure time for the control group.

  - `rateRatio`: The rate ratio of the active treatment group versus the
    control group.

  - `vlogRate1`: The variance for the log rate parameter for the active
    treatment group.

  - `vlogRate2`: The variance for the log rate parameter for the control
    group.

  - `vlogRR`: The variance of log rate ratio.

  - `information`: The information of log rate ratio.

  - `zlogRR`: The Z-statistic for log rate ratio.

- `resultsUnderH0` when `nullVariance = TRUE`: A data frame with the
  following variables:

  - `time`: The analysis time since trial start.

  - `lambda1H0`: The restricted maximum likelihood estimate of the event
    rate for the active treatment group.

  - `lambda2H0`: The restricted maximum likelihood estimate of the event
    rate for the control group.

  - `rateRatioH0`: The rate ratio under H0.

  - `vlogRate1H0`: The variance for the log rate parameter for the
    active treatment group under H0.

  - `vlogRate2H0`: The variance for the log rate parameter for the
    control group under H0.

  - `vlogRRH0`: The variance of log rate ratio under H0.

  - `informationH0`: The information of log rate ratio under H0.

  - `zlogRRH0`: The Z-statistic for log rate ratio with variance
    evaluated under H0.

  - `varianceRatio`: The ratio of the variance under H0 versus the
    variance under H1.

  - `lambda1`: The true event rate for the active treatment group.

  - `lambda2`: The true event rate for the control group.

  - `rateRatio`: The true rate ratio.

- `resultsUnderH0` when `nullVariance = FALSE`: A data frame with the
  following variables:

  - `time`: The analysis time since trial start.

  - `rateRatioH0`: The rate ratio under H0.

  - `varianceRatio`: Equal to 1.

  - `lambda1`: The true event rate for the active treatment group.

  - `lambda2`: The true event rate for the control group.

  - `rateRatio`: The true rate ratio.

## Details

The probability mass function for a negative binomial distribution with
dispersion parameter \\\kappa_i\\ and rate parameter \\\lambda_i\\ is
given by \$\$P(Y\_{ij} = y) =
\frac{\Gamma(y+1/\kappa_i)}{\Gamma(1/\kappa_i) y!} \left(\frac{1}{1 +
\kappa_i \lambda_i t\_{ij}}\right)^{1/\kappa_i} \left(\frac{\kappa_i
\lambda_i t\_{ij}} {1 + \kappa_i \lambda_i t\_{ij}}\right)^{y},\$\$
where \\Y\_{ij}\\ is the event count for subject \\j\\ in treatment
group \\i\\, and \\t\_{ij}\\ is the exposure time for the subject. If
\\\kappa_i=0\\, the negative binomial distribution reduces to the
Poisson distribution.

For treatment group \\i\\, let \\\beta_i = \log(\lambda_i)\\. The
log-likelihood for \\\\(\kappa_i, \beta_i):i=1,2\\\\ can be written as
\$\$l = \sum\_{i=1}^{2}\sum\_{j=1}^{n\_{i}} \\\log \Gamma(y\_{ij} +
1/\kappa_i) - \log \Gamma(1/\kappa_i) + y\_{ij} (\log(\kappa_i) +
\beta_i) - (y\_{ij} + 1/\kappa_i) \log(1+ \kappa_i \exp(\beta_i)
t\_{ij})\\.\$\$ It follows that \$\$\frac{\partial l}{\partial \beta_i}
= \sum\_{j=1}^{n_i} \left\\y\_{ij} - (y\_{ij} + 1/\kappa_i)
\frac{\kappa_i \exp(\beta_i) t\_{ij}} {1 + \kappa_i
\exp(\beta_i)t\_{ij}}\right\\,\$\$ and \$\$-\frac{\partial^2 l}{\partial
\beta_i^2} = \sum\_{j=1}^{n_i} (y\_{ij} + 1/\kappa_i) \frac{\kappa_i
\lambda_i t\_{ij}} {(1 + \kappa_i \lambda_i t\_{ij})^2}.\$\$ The Fisher
information for \\\beta_i\\ is \$\$E\left(-\frac{\partial^2 l}{\partial
\beta_i^2}\right) = n_i E\left(\frac{\lambda_i t\_{ij}} {1 + \kappa_i
\lambda_i t\_{ij}}\right).\$\$ In addition, we can show that
\$\$E\left(-\frac{\partial^2 l} {\partial \beta_i \partial
\kappa_i}\right) = 0.\$\$ Therefore, the variance of \\\hat{\beta}\_i\\
is \$\$Var(\hat{\beta}\_i) = \frac{1}{n_i} \left\\
E\left(\frac{\lambda_i t\_{ij}}{1 + \kappa_i \lambda_i t\_{ij}}\right)
\right\\^{-1}.\$\$

To evaluate the integral, we need to obtain the distribution of the
exposure time, \$\$t\_{ij} = \min(\tau - W\_{ij}, C\_{ij},
T\_{fmax}),\$\$ where \\\tau\\ denotes the calendar time since trial
start, \\W\_{ij}\\ denotes the enrollment time for subject \\j\\ in
treatment group \\i\\, \\C\_{ij}\\ denotes the time to dropout after
enrollment for subject \\j\\ in treatment group \\i\\, and \\T\_{fmax}\\
denotes the maximum follow-up time for all subjects. Therefore,
\$\$P(t\_{ij} \geq t) = P(W\_{ij} \leq \tau - t)P(C\_{ij} \geq t)
I(t\leq T\_{fmax}).\$\$ Let \\H\\ denote the distribution function of
the enrollment time, and \\G_i\\ denote the survival function of the
dropout time for treatment group \\i\\. By the change of variables, we
have \$\$E\left(\frac{\lambda_i t\_{ij}}{1 + \kappa_i \lambda_i t\_{ij}}
\right) = \int\_{0}^{\tau \wedge T\_{fmax}} \frac{\lambda_i}{(1 +
\kappa_i \lambda_i t)^2} H(\tau - t) G_i(t) dt.\$\$ A numerical
integration algorithm for a univariate function can be used to evaluate
the above integral.

For the restricted maximum likelihood (reml) estimate of
\\(\beta_1,\beta_2)\\ subject to the constraint that \\\beta_1 - \beta_2
= \Delta\\, we express the log-likelihood in terms of
\\(\beta_2,\Delta,\kappa_1,\kappa_2)\\, and takes the derivative of the
log-likelihood function with respect to \\\beta_2\\. The resulting score
equation has asymptotic limit \$\$E\left(\frac{\partial l}{\partial
\beta_2}\right) = s_1 + s_2,\$\$ where \$\$s_1 = n r E\left\\\lambda_1
t\_{1j} - \left(\lambda_1t\_{1j} + \frac{1}{\kappa_1}\right)
\frac{\kappa_1 e^{\tilde{\beta}\_2 + \Delta}t\_{1j}}{1 + \kappa_1
e^{\tilde{\beta}\_2 +\Delta}t\_{1j}}\right\\,\$\$ and \$\$s_2 = n (1-r)
E\left\\\lambda_2 t\_{2j} - \left(\lambda_2 t\_{2j} +
\frac{1}{\kappa_2}\right) \frac{\kappa_2 e^{\tilde{\beta}\_2} t\_{2j}}
{1 + \kappa_2 e^{\tilde{\beta}\_2}t\_{2j}}\right\\.\$\$ Here \\r\\ is
the randomization probability for the active treatment group. The
asymptotic limit of the reml of \\\beta_2\\ is the solution
\\\tilde{\beta}\_2\\ to \\E\left(\frac{\partial l}{\partial
\beta_2}\right) = 0.\\

## Author

Kaifeng Lu, <kaifenglu@gmail.com>

## Examples

``` r
# Example 1: Variable follow-up design

nbstat(time = c(1, 1.25, 2, 3, 4),
       accrualIntensity = 1956/1.25,
       kappa1 = 5,
       kappa2 = 5,
       lambda1 = 0.7*0.125,
       lambda2 = 0.125,
       gamma1 = 0,
       gamma2 = 0,
       accrualDuration = 1.25,
       followupTime = 2.75)
#> $resultsUnderH1
#>   time subjects  nevents  nevents1  nevents2 ndropouts ndropouts1 ndropouts2
#> 1 1.00   1564.8  83.1300  34.23000  48.90000         0          0          0
#> 2 1.25   1956.0 129.8906  53.48438  76.40625         0          0          0
#> 3 2.00   1956.0 285.7594 117.66562 168.09375         0          0          0
#> 4 3.00   1956.0 493.5844 203.24062 290.34375         0          0          0
#> 5 4.00   1956.0 701.4094 288.81562 412.59375         0          0          0
#>   nfmax nfmax1 nfmax2 exposure exposure1 exposure2 rateRatio   vlogRate1
#> 1     0      0      0    782.4    391.20    391.20       0.7 0.037481104
#> 2     0      0      0   1222.5    611.25    611.25       0.7 0.025270509
#> 3     0      0      0   2689.5   1344.75   1344.75       0.7 0.013838648
#> 4     0      0      0   4645.5   2322.75   2322.75       0.7 0.010091604
#> 5     0      0      0   6601.5   3300.75   3300.75       0.7 0.008598729
#>     vlogRate2     vlogRR information    zlogRR
#> 1 0.028633294 0.06611440    15.12530 -1.387154
#> 2 0.019585298 0.04485581    22.29366 -1.684082
#> 3 0.011259561 0.02509821    39.84348 -2.251393
#> 4 0.008605162 0.01869677    53.48519 -2.608491
#> 5 0.007555189 0.01615392    61.90449 -2.806297
#> 
#> $resultsUnderH0
#>   time rateRatioH0 varianceRatio lambda1 lambda2 rateRatio
#> 1 1.00           1             1  0.0875   0.125       0.7
#> 2 1.25           1             1  0.0875   0.125       0.7
#> 3 2.00           1             1  0.0875   0.125       0.7
#> 4 3.00           1             1  0.0875   0.125       0.7
#> 5 4.00           1             1  0.0875   0.125       0.7
#> 

# Example 2: Fixed follow-up design

nbstat(time = c(0.5, 1, 1.5, 2),
       accrualIntensity = 220/1.5,
       stratumFraction = c(0.2, 0.8),
       kappa1 = 3,
       kappa2 = 3,
       lambda1 = c(0.5*8.4, 0.6*10.5),
       lambda2 = c(8.4, 10.5),
       gamma1 = 0,
       gamma2 = 0,
       accrualDuration = 1.5,
       followupTime = 0.5,
       fixedFollowup = 1,
       nullVariance = 1)
#> $resultsUnderH1
#>   time  subjects nevents nevents1 nevents2 ndropouts ndropouts1 ndropouts2
#> 1  0.5  73.33333   146.3     53.9     92.4         0          0          0
#> 2  1.0 146.66667   438.9    161.7    277.2         0          0          0
#> 3  1.5 220.00000   731.5    269.5    462.0         0          0          0
#> 4  2.0 220.00000   877.8    323.4    554.4         0          0          0
#>       nfmax    nfmax1    nfmax2  exposure exposure1 exposure2 rateRatio
#> 1   0.00000   0.00000   0.00000  18.33333  9.166667  9.166667 0.5785155
#> 2  73.33333  36.66667  36.66667  55.00000 27.500000 27.500000 0.5785155
#> 3 146.66667  73.33333  73.33333  91.66667 45.833333 45.833333 0.5785155
#> 4 220.00000 110.00000 110.00000 110.00000 55.000000 55.000000 0.5785155
#>    vlogRate1  vlogRate2     vlogRR information    zlogRR
#> 1 0.11098462 0.10035910 0.21134372    4.731629 -1.190482
#> 2 0.05010035 0.04667902 0.09677937   10.332781 -1.759244
#> 3 0.03235374 0.03041238 0.06276613   15.932160 -2.184514
#> 4 0.03044733 0.02909091 0.05953824   16.795928 -2.242949
#> 
#> $resultsUnderH0
#>   time lambda1H0 lambda2H0 rateRatioH0 vlogRate1H0 vlogRate2H0   vlogRRH0
#> 1  0.5  7.930335  7.930335           1  0.10432521  0.10432521 0.20865042
#> 2  1.0  7.930335  7.930335           1  0.04795146  0.04795146 0.09590292
#> 3  1.5  7.930335  7.930335           1  0.03113041  0.03113041 0.06226083
#> 4  2.0  7.930335  7.930335           1  0.02958153  0.02958153 0.05916306
#>   informationH0  zlogRRH0 varianceRatio lambda1 lambda2 rateRatio
#> 1      4.792705 -1.198141     0.9872563 5.80928 10.0417 0.5785155
#> 2     10.427211 -1.767264     0.9909439 5.80928 10.0417 0.5785155
#> 3     16.061463 -2.193360     0.9919495 5.80928 10.0417 0.5785155
#> 4     16.902439 -2.250050     0.9936985 5.80928 10.0417 0.5785155
#> 
```
