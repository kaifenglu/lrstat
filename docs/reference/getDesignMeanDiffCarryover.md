# Power and Sample Size for Direct Treatment Effects in Crossover Trials

Obtains the power and sample size for direct treatment effects in
crossover trials accounting or without accounting for carryover effects.

## Usage

``` r
getDesignMeanDiffCarryover(
  beta = NA_real_,
  n = NA_real_,
  trtpair = NA_real_,
  carryover = TRUE,
  meanDiffH0 = 0,
  meanDiff = 0.5,
  stDev = 1,
  corr = 0.5,
  design = NA_real_,
  cumdrop = NA_real_,
  allocationRatioPlanned = NA_real_,
  normalApproximation = FALSE,
  rounding = TRUE,
  alpha = 0.025
)
```

## Arguments

- beta:

  The type II error.

- n:

  The total sample size.

- trtpair:

  The treatment pair of interest to power the study. If not given, it
  defaults to comparing the first treatment to the last treatment.

- carryover:

  Whether to account for carryover effects in the power calculation.
  Defaults to TRUE.

- meanDiffH0:

  The mean difference for the treatment pair of interest under the null
  hypothesis. Defaults to 0.

- meanDiff:

  The mean difference for the treatment pair of interest under the
  alternative hypothesis.

- stDev:

  The standard deviation for within-subject random error.

- corr:

  The intra-subject correlation due to subject random effect.

- design:

  The crossover design represented by a matrix with rows indexing the
  sequences, columns indexing the periods, and matrix entries indicating
  the treatments.

- cumdrop:

  The cumulative dropout rate over periods.

- allocationRatioPlanned:

  Allocation ratio for the sequences. Defaults to equal randomization if
  not provided.

- normalApproximation:

  The type of computation of the p-values. If `TRUE`, the variance is
  assumed to be known, otherwise the calculations are performed with the
  t distribution.

- rounding:

  Whether to round up the sample size. Defaults to TRUE for sample size
  rounding.

- alpha:

  The one-sided significance level. Defaults to 0.025.

## Value

An S3 class `designMeanDiffCarryover` object with the following
components:

- `power`: The power to reject the null hypothesis.

- `alpha`: The one-sided significance level.

- `numberOfSubjects`: The maximum number of subjects.

- `trtpair`: The treatment pair of interest to power the study.

- `carryover`: Whether to account for carryover effects in the power
  calculation.

- `meanDiffH0`: The mean difference for the treatment pair of interest
  under the null hypothesis.

- `meanDiff`: The mean difference for the treatment pair of interest
  under the alternative hypothesis.

- `stDev`: The standard deviation for within-subject random error.

- `corr`: The intra-subject correlation due to subject random effect.

- `design`: The crossover design represented by a matrix with rows
  indexing the sequences, columns indexing the periods, and matrix
  entries indicating the treatments.

- `designMatrix`: The design matrix accounting for intercept, sequence,
  period, direct treatment effects and carryover treatment effects when
  `carryover = TRUE`, or the design matrix accounting for intercept,
  sequence, period, and direct treatment effects when
  `carryover = FALSE`.

- `nseq`: The number of sequences.

- `nprd`: The number of periods.

- `ntrt`: The number of treatments.

- `cumdrop`: The cumulative dropout rate over periods.

- `V_direct_only`: The covariance matrix for direct treatment effects
  without accounting for carryover effects. The treatment comparisons
  for the covariance matrix are for the first \\t-1\\ treatments
  relative to the last treatment.

- `V_direct_carry`: The covariance matrix for direct and carryover
  treatment effects.

- `v_direct_only`: The variance of the direct treatment effect for the
  treatment pair of interest without accounting for carryover effects.

- `v_direct`: The variance of the direct treatment effect for the
  treatment pair of interest accounting for carryover effects.

- `v_carry`: The variance of the carryover treatment effect for the
  treatment pair of interest.

- `releff_direct`: The relative efficiency of the design for estimating
  the direct treatment effect for the treatment pair of interest after
  accounting for carryover effects with respect to that without
  accounting for carryover effects. This is equal to
  `v_direct_only/v_direct`.

- `releff_carry`: The relative efficiency of the design for estimating
  the carryover effect for the treatment pair of interest. This is equal
  to `v_direct_only/v_carry`.

- `half_width`: The half-width of the confidence interval for the direct
  treatment effect for the treatment pair of interest.

- `nu`: Degrees of freedom for the t-test.

- `allocationRatioPlanned`: Allocation ratio for the sequences.

- `normalApproximation`: The type of computation of the p-values. If
  `TRUE`, the variance is assumed to be known, otherwise the
  calculations are performed with the t distribution.

- `rounding`: Whether to round up the sample size.

## Details

The linear mixed-effects model to assess the direct treatment effects in
the presence of carryover treatment effects is given by \$\$y\_{ijk} =
\mu + \alpha_i + b\_{ij} + \gamma_k + \tau\_{d(i,k)} +
\lambda\_{c(i,k-1)} + e\_{ijk}\$\$ \$\$i=1,\ldots,n; j=1,\ldots,r_i; k =
1,\ldots,p; d,c = 1,\ldots,t\$\$ where \\\mu\\ is the general mean,
\\\alpha_i\\ is the effect of the \\i\\th treatment sequence,
\\b\_{ij}\\ is the random effect with variance \\\sigma_b^2\\ for the
\\j\\th subject of the \\i\\th treatment sequence, \\\gamma_k\\ is the
period effect, and \\e\_{ijk}\\ is the random error with variance
\\\sigma^2\\ for the subject in period \\k\\. The direct effect of the
treatment administered in period \\k\\ of sequence \\i\\ is
\\\tau\_{d(i,k)}\\, and \\\lambda\_{c(i,k-1)}\\ is the carryover effect
of the treatment administered in period \\k-1\\ of sequence \\i\\. The
value of the carryover effect for the observed response in the first
period is \\\lambda\_{c(i,0)} = 0\\ since there is no carryover effect
in the first period. The intra-subject correlation due to the subject
random effect is \$\$\rho = \frac{\sigma_b^2}{\sigma_b^2 +
\sigma^2}.\$\$ Therefore, `stDev` = \\\sigma^2\\ and `corr` = \\\rho\\.
By constructing the design matrix \\X\\ for the linear model with a
compound symmetry covariance matrix for the response vector of a
subject, we can obtain \$\$Var(\hat{\beta}) = (X'V^{-1}X)^{-1}.\$\$

The covariance matrix for the direct treatment effects and carryover
treatment effects can be extracted from the appropriate sub-matrices.
The covariance matrix for the direct treatment effects without
accounting for the carryover treatment effects can be obtained by
omitting the carryover effect terms from the model.

The power is for the direct treatment effect for the treatment pair of
interest with or without accounting for carryover effects as determined
by the input parameter `carryover`. The relative efficiency is for the
direct treatment effect for the treatment pair of interest accounting
for carryover effects relative to that without accounting for carryover
effects.

The degrees of freedom for the t-test accounting for carryover effects
can be calculated as the total number of observations minus the number
of subjects minus \\p-1\\ minus \\2(t-1)\\ to account for the subject
effect, period effect, and direct and carryover treatment effects. The
degrees of freedom for the t-test without accounting for carryover
effects is equal to the total number of observations minus the number of
subjects minus \\p-1\\ minus \\t-1\\ to account for the subject effect,
period effect, and direct treatment effects.

## References

Robert O. Kuehl. Design of Experiments: Statistical Principles of
Research Design and Analysis. Brooks/Cole: Pacific Grove, CA. 2000.

## Author

Kaifeng Lu, <kaifenglu@gmail.com>

## Examples

``` r
# Williams design for 4 treatments

(design1 = getDesignMeanDiffCarryover(
  beta = 0.2, n = NA,
  meanDiff = 0.5, stDev = 1,
  design = matrix(c(1, 4, 2, 3,
                    2, 1, 3, 4,
                    3, 2, 4, 1,
                    4, 3, 1, 2),
                  4, 4, byrow = TRUE),
  alpha = 0.025))
#>                                                                   
#> Testing direct treatment effects accounting for carryover effects 
#>      Prd1 Prd2 Prd3 Prd4
#> Seq1    1    4    2    3
#> Seq2    2    1    3    4
#> Seq3    3    2    4    1
#> Seq4    4    3    1    2
#>                                                                                   
#> Treatment comparison of interest: 1 - 4                                           
#> Mean difference under H0: 0, mean difference under H1: 0.5                        
#> Within-subject standard deviation: 1, intra-subject correlation: 0.5              
#> Power: 0.8015, alpha (1-sided): 0.025, CI half-width: 0.35                        
#> Cumulative dropout rates over periods: 0, 0, 0, 0                                 
#> Without accounting for carryover effects, variance for direct effect: 0.029       
#> Accounting for carryover, variance for direct effect: 0.031, for carryover: 0.046 
#> Relative efficiency for direct effect: 0.909, for carryover effect: 0.625         
#> Number of subjects: 70                                                            
#> Sequence allocation ratio: 1 1 1 1                                                
#> Test statistic: t-test                                                            
```
