# Power and Sample Size Calculation for Non-Proportional Hazards and Beyond

A comprehensive clinical trial design and analysis package with a focus
on non-proportional hazards and weighted log-rank methods, plus broad
support for group sequential design, adaptive design, multiplicity,
dose-finding, and endpoint-specific power and sample size calculations.

## Details

For proportional hazards, power is determined by the total number of
events and the constant hazard ratio together with information rates and
spending functions. For non-proportional hazards, the hazard ratio
varies over time and calendar time determines the mean and variance of
the weighted log-rank score statistic. The package uses the analytic
approach in Lu (2021) and approximates the variance-covariance matrix of
sequential statistics under the alternative by that under the null to
leverage the independent increments structure in Tsiatis (1982) for the
Fleming-Harrington family.

The package capabilities can be grouped as follows:

- **Time-to-event design under proportional and non-proportional
  hazards:** weighted log-rank statistics and design operating
  characteristics (for example `lrstat`, `lrpower`, `lrsamplesize`,
  `lrsim`, `lrschoenfeld`), including equivalence settings, event
  prediction, accrual modeling, and calendar time determination (for
  example `accrual`, `caltime`, `nevent`, `pevent`, `patrisk`,
  `natrisk`, `getDurationFromNevents`, `getNeventsFromHazardRatio`).

- **Group sequential and adaptive designs:** standard, adaptive,
  multi-arm multi-stage (MAMS), and seamless frameworks, including
  boundaries, conditional power, confidence intervals, and exit
  probabilities (for example `getDesign`, `adaptDesign`, `getBound`,
  `getCP`, `getCI`, `getRCI`, `getADCI`, `getADRCI`, `exitprob`, and
  their `_mams` and `_seamless` variants).

- **Fixed and sequential sample size/power across endpoint types:**
  binary, continuous, ordinal/multinomial, paired, crossover/carryover,
  MMRM, repeated measures ANOVA, two-way ANOVA, Wilcoxon, Fisher's
  exact, McNemar, one-sample exact, agreement, logistic regression, and
  negative binomial endpoints (for example `getDesign*`, `power*Exact`,
  `samplesize*Exact`, `kmpower`, `rmpower`, `nbpower`, and associated
  `*samplesize` and `*equiv` functions).

- **Multiplicity and graphical procedures:** weighted Bonferroni and
  graph-based updates, Bonferroni mixtures of weighted Simes and
  Dunnett, and gatekeeping procedures for multiple hypotheses (for
  example `fadjpbon`, `fadjpsim`, `fadjpdun`, `fseqbon`, `fstp2seq`,
  `fstdmix`, `fmodmix`, `fwgtmat`, `updateGraph`).

- **Early-phase and dose-finding tools:** Simon's two-stage, Bayesian
  Simon designs, mTPI-2, and BOIN utilities (for example `simon2stage`,
  `simonBayesAnalysis`, `simonBayesSim`, `mTPI2Table`, `BOINTable`).

- **Estimation, confidence intervals, and supporting utilities:**
  Clopper-Pearson, exact risk ratio/risk difference, Brookmeyer-Crowley
  survival quantiles, Miettinen-Nurminen and REML intervals/statistics
  for stratified measures, Hedges' g, piecewise-exponential distribution
  tools, multivariate normal integration, and model diagnostics (for
  example `ClopperPearsonCI`, `riskRatioExactCI`, `riskDiffExactCI`,
  `survQuantile`, `mnRiskDiffCI`, `remlRiskDiff`, `hedgesg`, `dtpwexp`,
  `ptpwexp`, `qtpwexp`, `rtpwexp`, `pbvnorm`, `pmvnormr`, `qmvnormr`,
  `phregr`, `liferegr`, `residuals_phregr`, `survfit_phregr`,
  `zph_phregr`).

- **Interactive use:** a bundled Shiny interface for exploratory design
  workflows (`runShinyApp_lrstat`).

The development of lrstat was strongly influenced by `rpact`, with
argument naming aligned where possible for ease of adoption. Key
differences include direct approximation (rather than Schoenfeld-based
approximation) for weighted log-rank design calculations, explicit use
of `accrualDuration` to define the accrual end, and treatment of final
stage trial outcome when early futility stopping does not occur.

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
