# Power and Sample Size Calculation for Non-Proportional Hazards and Beyond

Performs power and sample size calculation for non-proportional hazards
model using the Fleming-Harrington family of weighted log-rank tests.

## Details

For proportional hazards, the power is determined by the total number of
events and the constant hazard ratio along with information rates and
spending functions. For non-proportional hazards, the hazard ratio
varies over time and the calendar time plays a key role in determining
the mean and variance of the log-rank test score statistic. It requires
an iterative algorithm to find the calendar time at which the targeted
number of events will be reached for each interim analysis. The lrstat
package uses the analytic method in Lu (2021) to find the mean and
variance of the weighted log-rank test score statistic at each interim
analysis. In addition, the package approximates the variance and
covariance matrix of the sequentially calculated log-rank test
statistics under the alternative hypothesis with that under the null
hypothesis to take advantage of the independent increments structure in
Tsiatis (1982) applicable for the Fleming-Harrington family of weighted
log-rank tests.

The most useful functions in the package are `lrstat`, `lrpower`,
`lrsamplesize`, and `lrsim`, which calculate the mean and variance of
log-rank test score statistic at a sequence of given calendar times, the
power of the log-rank test, the sample size in terms of accrual duration
and follow-up duration, and the log-rank test simulation, respectively.
The `accrual` function calculates the number of patients accrued at
given calendar times. The `caltime` function finds the calendar times to
reach the targeted number of events. The `exitprob` function calculates
the stagewise exit probabilities for specified boundaries with a varying
mean parameter over time based on an adaptation of the recursive
integration algorithm described in Chapter 19 of Jennison and Turnbull
(2000).

The development of the lrstat package is heavily influenced by the rpact
package. We find their function arguments to be self-explanatory. We
have used the same names whenever appropriate so that users familiar
with the rpact package can learn the lrstat package quickly. However,
there are notable differences:

- lrstat uses direct approximation, while rpact uses the Schoenfeld
  method for log-rank test power and sample size calculation.

- lrstat uses `accrualDuration` to explicitly set the end of accrual
  period, while rpact incorporates the end of accrual period in
  `accrualTime`.

- lrstat considers the trial a failure at the last stage if the log-rank
  test cannot reject the null hypothesis up to this stage and cannot
  stop for futility at an earlier stage.

- the `lrsim` function uses the variance of the log-rank test score
  statistic as the information.

In addition to the log-rank test power and sample size calculations, the
lrstat package can also be used for the following tasks:

- design generic group sequential trials.

- design generic group sequential equivalence trials.

- design adaptive group sequential trials for changes in sample size,
  error spending function, number and spacing or future looks.

- calculate the terminating and repeated confidence intervals for
  standard and adaptive group sequential trials.

- calculate the conditional power for non-proportional hazards with or
  without design changes.

- perform multiplicity adjustment based on graphical approaches using
  weighted Bonferroni tests, Bonferroni mixture of weighted Simes test,
  and Bonferroni mixture of Dunnett test as well as group sequential
  trials with multiple hypotheses.

- perform multiplicity adjustment using stepwise gatekeeping procedures
  for two sequences of hypotheses and the standard or modified mixture
  gatekeeping procedures in the general case.

- design parallel-group trials with the primary endpoint analyzed using
  mixed-model for repeated measures (MMRM).

- design crossover trials to estimate direct treatment effects while
  accounting for carryover effects.

- design one-way repeated measures ANOVA trials.

- design two-way ANOVA trials.

- design Simon's 2-stage trials.

- design modified toxicity probability-2 (mTPI-2) trials.

- design Bayesian optimal interval (BOIN) trials.

- design group sequential trials for negative binomial endpoints with
  censoring.

- design trials using Wilcoxon, Fisher's exact, and McNemar's test.

- calculate Clopper-Pearson confidence interval for single proportions.

- calculate Brookmeyer-Crowley confidence interval for quantiles of
  censored survival data.

- calculate Miettinen & Nurminen confidence interval for stratified risk
  difference, risk ratio, odds ratio, rate difference, and rate ratio.

- perform power and sample size calculation for logistic regression.

- perform power and sample size calculation for Cohen's kappa.

- calculate Hedges' g effect size.

- generate random numbers from truncated piecewise exponential
  distribution.

- perform power and sample size calculations for negative binomial data.

## References

Anastasios A. Tsiatis. Repeated significance testing for a general class
of statistics used in censored survival analysis. J Am Stat Assoc.
1982;77:855-861.

Christopher Jennison, Bruce W. Turnbull. Group Sequential Methods with
Applications to Clinical Trials. Chapman & Hall/CRC: Boca Raton, 2000,
ISBN:0849303168

Kaifeng Lu. Sample size calculation for logrank test and prediction of
number of events over time. Pharm Stat. 2021;20:229-244.

## See also

rpact, gsDesign

## Author

Kaifeng Lu, <kaifenglu@gmail.com>

## Examples

``` r
lrpower(kMax = 2, informationRates = c(0.8, 1),
        criticalValues = c(2.250, 2.025), accrualIntensity = 20,
        piecewiseSurvivalTime = c(0, 6),
        lambda1 = c(0.0533, 0.0309), lambda2 = c(0.0533, 0.0533),
        gamma1 = 0.00427, gamma2 = 0.00427,
        accrualDuration = 22, followupTime = 18)
#>                                                                      
#> Group-sequential design with 2 stages for log-rank test              
#> Overall power: 0.7786, overall significance level (1-sided): 0.025   
#> Maximum # events: 296.6, expected # events: 271.1                    
#> Maximum # dropouts: 28.1, expected # dropouts: 25.4                  
#> Maximum # subjects: 440, expected # subjects: 440                    
#> Maximum information: 73.36, expected information: 67.23              
#> Total study duration: 40, expected study duration: 35.7              
#> Accrual duration: 22, follow-up duration: 18, fixed follow-up: FALSE 
#> Allocation ratio: 1                                                  
#> Alpha spending: Lan-DeMets O'Brien-Fleming, beta spending: None      
#>                                                                      
#>                        Stage 1 Stage 2
#> Information rate       0.800   1.000  
#> Efficacy boundary (Z)  2.250   2.025  
#> Cumulative rejection   0.4286  0.7786 
#> Cumulative alpha spent 0.0122  0.0250 
#> Number of events       237.2   296.6  
#> Number of dropouts     21.8    28.1   
#> Number of subjects     440.0   440.0  
#> Analysis time          29.9    40.0   
#> Efficacy boundary (HR) 0.745   0.789  
#> Efficacy boundary (p)  0.0122  0.0214 
#> Information            59.06   73.36  
#> HR                     0.764   0.722  
```
