# Group Sequential Design for Two-Sample Mean Difference From the MMRM Model

Obtains the power and sample size for two-sample mean difference at the
last time point from the mixed-model for repeated measures (MMRM) model.

## Usage

``` r
getDesignMeanDiffMMRM(
  beta = NA_real_,
  meanDiffH0 = 0,
  meanDiff = 0.5,
  k = 1,
  t = NA_real_,
  covar1 = diag(k),
  covar2 = NA_real_,
  accrualTime = 0,
  accrualIntensity = NA_real_,
  piecewiseSurvivalTime = 0,
  gamma1 = 0,
  gamma2 = 0,
  accrualDuration = NA_real_,
  allocationRatioPlanned = 1,
  normalApproximation = TRUE,
  rounding = TRUE,
  kMax = 1L,
  informationRates = NA_real_,
  efficacyStopping = NA_integer_,
  futilityStopping = NA_integer_,
  criticalValues = NA_real_,
  alpha = 0.025,
  typeAlphaSpending = "sfOF",
  parameterAlphaSpending = NA_real_,
  userAlphaSpending = NA_real_,
  futilityBounds = NA_real_,
  futilityCP = NA_real_,
  futilityMeanDiff = NA_real_,
  typeBetaSpending = "none",
  parameterBetaSpending = NA_real_,
  userBetaSpending = NA_real_,
  spendingTime = NA_real_
)
```

## Arguments

- beta:

  The type II error.

- meanDiffH0:

  The mean difference at the last time point under the null hypothesis.
  Defaults to 0.

- meanDiff:

  The mean difference at the last time point under the alternative
  hypothesis.

- k:

  The number of postbaseline time points.

- t:

  The postbaseline time points.

- covar1:

  The covariance matrix for the repeated measures given baseline for the
  active treatment group.

- covar2:

  The covariance matrix for the repeated measures given baseline for the
  control group. If missing, it will be set equal to the covariance
  matrix for the active treatment group.

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

- gamma1:

  The hazard rate for exponential dropout, or a vector of hazard rates
  for piecewise exponential dropout for the active treatment group.

- gamma2:

  The hazard rate for exponential dropout, or a vector of hazard rates
  for piecewise exponential dropout for the control group.

- accrualDuration:

  Duration of the enrollment period.

- allocationRatioPlanned:

  Allocation ratio for the active treatment versus control. Defaults to
  1 for equal randomization.

- normalApproximation:

  The type of computation of the p-values. If `TRUE`, the variance is
  assumed to be known, otherwise the calculations are performed with the
  t distribution. The degrees of freedom for the t-distribution is the
  total effective sample size minus 2. The exact calculation using the t
  distribution is only implemented for the fixed design.

- rounding:

  Whether to round up sample size. Defaults to 1 for sample size
  rounding.

- kMax:

  The maximum number of stages.

- informationRates:

  The information rates. Defaults to `(1:kMax) / kMax` if left
  unspecified.

- efficacyStopping:

  Indicators of whether efficacy stopping is allowed at each stage.
  Defaults to `TRUE` if left unspecified.

- futilityStopping:

  Indicators of whether futility stopping is allowed at each stage.
  Defaults to `TRUE` if left unspecified.

- criticalValues:

  Upper boundaries on the z-test statistic scale for stopping for
  efficacy.

- alpha:

  The significance level. Defaults to 0.025.

- typeAlphaSpending:

  The type of alpha spending. One of the following: `"OF"` for
  O'Brien-Fleming boundaries, `"P"` for Pocock boundaries, `"WT"` for
  Wang & Tsiatis boundaries, `"sfOF"` for O'Brien-Fleming type spending
  function, `"sfP"` for Pocock type spending function, `"sfKD"` for Kim
  & DeMets spending function, `"sfHSD"` for Hwang, Shi & DeCani spending
  function, `"user"` for user defined spending, and `"none"` for no
  early efficacy stopping. Defaults to `"sfOF"`.

- parameterAlphaSpending:

  The parameter value for the alpha spending. Corresponds to \\\Delta\\
  for `"WT"`, \\\rho\\ for `"sfKD"`, and \\\gamma\\ for `"sfHSD"`.

- userAlphaSpending:

  The user defined alpha spending. Cumulative alpha spent up to each
  stage.

- futilityBounds:

  Lower boundaries on the z-test statistic scale for stopping for
  futility at stages `1, ..., kMax-1`. Defaults to `rep(-6, kMax-1)` if
  left unspecified. The futility bounds are non-binding for the
  calculation of critical values.

- futilityCP:

  The futility bounds on the conditional power scale.

- futilityMeanDiff:

  The futility bounds on the mean difference scale.

- typeBetaSpending:

  The type of beta spending. One of the following: `"sfOF"` for
  O'Brien-Fleming type spending function, `"sfP"` for Pocock type
  spending function, `"sfKD"` for Kim & DeMets spending function,
  `"sfHSD"` for Hwang, Shi & DeCani spending function, `"user"` for user
  defined spending, and `"none"` for no early futility stopping.
  Defaults to `"none"`.

- parameterBetaSpending:

  The parameter value for the beta spending. Corresponds to \\\rho\\ for
  `"sfKD"`, and \\\gamma\\ for `"sfHSD"`.

- userBetaSpending:

  The user defined beta spending. Cumulative beta spent up to each
  stage.

- spendingTime:

  A vector of length `kMax` for the error spending time at each
  analysis. Defaults to missing, in which case, it is the same as
  `informationRates`.

## Value

An S3 class `designMeanDiffMMRM` object with three components:

- `overallResults`: A data frame containing the following variables:

  - `overallReject`: The overall rejection probability.

  - `alpha`: The overall significance level.

  - `attainedAlpha`: The attained significance level, which is different
    from the overall significance level in the presence of futility
    stopping.

  - `kMax`: The number of stages.

  - `theta`: The parameter value.

  - `information`: The maximum information.

  - `expectedInformationH1`: The expected information under H1.

  - `expectedInformationH0`: The expected information under H0.

  - `drift`: The drift parameter, equal to `theta*sqrt(information)`.

  - `inflationFactor`: The inflation factor (relative to the fixed
    design).

  - `numberOfSubjects`: The maximum number of subjects.

  - `studyDuration`: The maximum study duration.

  - `expectedNumberOfSubjectsH1`: The expected number of subjects under
    H1.

  - `expectedNumberOfSubjectsH0`: The expected number of subjects under
    H0.

  - `expectedStudyDurationH1`: The expected study duration under H1.

  - `expectedStudyDurationH0`: The expected study duration under H0.

  - `accrualDuration`: The accrual duration.

  - `followupTime`: The follow-up time.

  - `fixedFollowup`: Whether a fixed follow-up design is used.

  - `meanDiffH0`: The mean difference under H0.

  - `meanDiff`: The mean difference under H1.

- `byStageResults`: A data frame containing the following variables:

  - `informationRates`: The information rates.

  - `efficacyBounds`: The efficacy boundaries on the Z-scale.

  - `futilityBounds`: The futility boundaries on the Z-scale.

  - `rejectPerStage`: The probability for efficacy stopping.

  - `futilityPerStage`: The probability for futility stopping.

  - `cumulativeRejection`: The cumulative probability for efficacy
    stopping.

  - `cumulativeFutility`: The cumulative probability for futility
    stopping.

  - `cumulativeAlphaSpent`: The cumulative alpha spent.

  - `efficacyP`: The efficacy boundaries on the p-value scale.

  - `futilityP`: The futility boundaries on the p-value scale.

  - `information`: The cumulative information.

  - `efficacyStopping`: Whether to allow efficacy stopping.

  - `futilityStopping`: Whether to allow futility stopping.

  - `rejectPerStageH0`: The probability for efficacy stopping under H0.

  - `futilityPerStageH0`: The probability for futility stopping under
    H0.

  - `cumulativeRejectionH0`: The cumulative probability for efficacy
    stopping under H0.

  - `cumulativeFutilityH0`: The cumulative probability for futility
    stopping under H0.

  - `efficacyMeanDiff`: The efficacy boundaries on the mean difference
    scale.

  - `futilityMeanDiff`: The futility boundaries on the mean difference
    scale.

  - `numberOfSubjects`: The number of subjects.

  - `numberOfCompleters`: The number of completers.

  - `analysisTime`: The average time since trial start.

- `settings`: A list containing the following input parameters:

  - `typeAlphaSpending`: The type of alpha spending.

  - `parameterAlphaSpending`: The parameter value for alpha spending.

  - `userAlphaSpending`: The user defined alpha spending.

  - `typeBetaSpending`: The type of beta spending.

  - `parameterBetaSpending`: The parameter value for beta spending.

  - `userBetaSpending`: The user defined beta spending.

  - `spendingTime`: The error spending time at each analysis.

  - `allocationRatioPlanned`: The allocation ratio for the active
    treatment versus control.

  - `accrualTime`: A vector that specifies the starting time of
    piecewise Poisson enrollment time intervals.

  - `accrualIntensity`: A vector of accrual intensities. One for each
    accrual time interval.

  - `piecewiseSurvivalTime`: A vector that specifies the starting time
    of piecewise exponential survival time intervals.

  - `gamma1`: The hazard rate for exponential dropout or a vector of
    hazard rates for piecewise exponential dropout for the active
    treatment group.

  - `gamma2`: The hazard rate for exponential dropout or a vector of
    hazard rates for piecewise exponential dropout for the control
    group.

  - `k`: The number of postbaseline time points.

  - `t`: The postbaseline time points.

  - `covar1`: The covariance matrix for the repeated measures given
    baseline for the active treatment group.

  - `covar2`: The covariance matrix for the repeated measures given
    baseline for the control group.

  - `normalApproximation`: The type of computation of the p-values. If
    `TRUE`, the variance is assumed to be known, otherwise the
    calculations are performed with the t distribution.

  - `rounding`: Whether to round up sample size.

## Details

Consider a longitudinal study with two treatment groups. The outcome is
measured at baseline and at \\k\\ postbaseline time points. For each
treatment group, the outcomes are assumed to follow a multivariate
normal distribution. Conditional on baseline, the covariance matrix of
the post-baseline outcomes is denoted by \\\Sigma_1\\ for the active
treatment group and \\\Sigma_2\\ for the control group. Let \\\mu_1\\
and \\\mu_2\\ denote the mean vectors of post-baseline outcomes for the
active and control groups, respectively. We are interested in testing
the null hypothesis \\H_0: \mu\_{1,k} - \mu\_{2,k} = \delta_0\\ against
the alternative \\H_1: \mu\_{1,k} - \mu\_{2,k} = \delta\\.

The study design is based on the information for treatment difference at
the last postbaseline time point. This information is given by \$\$I =
1/\text{Var}(\hat{\mu}\_{1,k} - \hat{\mu}\_{2,k})\$\$ In the presence of
monotone missing data, let \\p\_{g,1},\ldots,p\_{g,k}\\ denote the
proportions of subjects in observed data patterns 1 through \\k\\ for
treatment group \\g=1\\ (active) or 2 (control). A subject in pattern
\\j\\ has complete data up to time \\t_j\\, i.e., the observed outcomes
are \\y\_{i,1},\ldots,y\_{i,j}\\, with missing values for
\\y\_{i,j+1},\ldots,y\_{i,k}\\.

According to Lu et al. (2008), the information matrix for the
post-baseline mean vector in group \\g\\ is \$\$I_g = n \pi_g J_g\$\$
where \\\pi_g\\ is the proportion of subjects in group \\g\\, and
\$\$J_g = \sum\_{j=1}^k p\_{g,j} \left( \begin{array}{cc}
\Sigma\_{g,j}^{-1} & 0 \\ 0 & 0 \end{array}\right)\$\$ Here,
\\\Sigma\_{g,j}\\ is the leading \\j\times j\\ principal submatrix of
\\\Sigma_g\\. It follows that \$\$\text{Var}(\hat{\mu}\_{1,k} -
\hat{\mu}\_{2,k}) = \frac{1}{n}\left(\frac{1}{\pi_1} J_1^{-1}\[k,k\] +
\frac{1}{\pi_2} J_2^{-1}\[k,k\]\right)\$\$

The observed data pattern probabilities depend on the accrual and
dropout distributions. Let \\H(u)\\ denote the cumulative distribution
function of enrollment time \\u\\, \\G_g(t)\\ denote the survival
function of dropout time \\t\\ for treatment group \\g\\, and \\\tau\\
denote the calendar time at interim or final analysis. Then, for
\\j=1,\ldots,k-1\\, the probability that a subject in group \\g\\ falls
into observed data pattern \\j\\ is \$\$p\_{g,j} = H(\tau -
t_j)G_g(t_j) - H(\tau - t\_{j+1})G_g(t\_{j+1})\$\$ For the last pattern
(\\j=k\\, i.e., completers), \$\$p\_{g,k} = H(\tau - t_k)G_g(t_k)\$\$
For the final analysis, \\\tau\\ is the study duration, so \\H(\tau -
t_j) = 1\\ for all \\j\\. Therefore, the pattern probabilities depend
only on the dropout distribution: \$\$p\_{g,j} = G_g(t_j) -
G_g(t\_{j+1}), \quad j=1,\ldots,k-1\$\$ and \$\$p\_{g,k} = G_g(t_k)\$\$

Cumulative dropout probabilities at post-baseline time points can be
used to define a piecewise exponential dropout distribution. Let
\\F_g(t_j)\\ denote the cumulative probability of dropout by time
\\t_j\\ for treatment group \\g\\. The left endpoints of the piecewise
survival time intervals are given by \\t_0=0,t_1,\ldots,t\_{k-1}\\. The
hazard rate in the interval \\(t\_{j-1},t_j\]\\ is given by
\$\$\gamma\_{g,j} = -\log\left(\frac{1 - F_g(t_j)}{1 - F_g(t\_{j-1})}
\right) / (t_j - t\_{j-1}), \quad j=1,\ldots,k\$\$

## References

Kaifeng Lu, Xiaohui Luo, and Pei-Yun Chen. Sample size estimation for
repeated measures analysis in randomized clinical trials with missing
data. The International Journal of Biostatistics 2008; 14(1), Article 9.

## Author

Kaifeng Lu, <kaifenglu@gmail.com>

## Examples

``` r
# function to generate the AR(1) correlation matrix
ar1_cor <- function(n, corr) {
  exponent <- abs(matrix((1:n) - 1, n, n, byrow = TRUE) - ((1:n) - 1))
  corr^exponent
}

(design1 = getDesignMeanDiffMMRM(
  beta = 0.2,
  meanDiffH0 = 0,
  meanDiff = 0.5,
  k = 4,
  t = c(1,2,3,4),
  covar1 = ar1_cor(4, 0.7),
  accrualIntensity = 10,
  gamma1 = 0.02634013,
  gamma2 = 0.02634013,
  accrualDuration = NA,
  allocationRatioPlanned = 1,
  kMax = 3,
  alpha = 0.025,
  typeAlphaSpending = "sfOF"))
#>                                                                                          
#> Group-sequential design with 3 stages for two-sample mean difference from the MMRM model 
#> Mean difference under H0: 0, mean difference under H1: 0.5                               
#> Standard deviation for treatment: 1, standard deviation for control: 1                   
#> Overall power: 0.8021, overall alpha (1-sided): 0.025                                    
#> Drift parameter: 2.827, inflation factor: 1.013                                          
#> Maximum information: 31.97, expected under H1: 27.29, expected under H0: 31.9            
#> Maximum # subjects: 139, expected under H1: 132.5, expected under H0: 138.9              
#> Total study duration: 17.9, expected under H1: 15.6, expected under H0: 17.9             
#> Accrual duration: 13.9, follow-up duration: 4, fixed follow-up: TRUE                     
#> Allocation ratio: 1                                                                      
#> Alpha spending: Lan-DeMets O'Brien-Fleming, beta spending: None                          
#>                                                                                          
#>                               Stage 1 Stage 2 Stage 3
#> Information rate              0.333   0.667   1.000  
#> Efficacy boundary (Z)         3.710   2.511   1.993  
#> Cumulative rejection          0.0189  0.4199  0.8021 
#> Cumulative alpha spent        0.0001  0.0060  0.0250 
#> Number of subjects            80.0    125.5   139.0  
#> Number of completers          40.0    85.5    139.0  
#> Analysis time                 8.0     12.5    17.9   
#> Efficacy boundary (mean diff) 1.137   0.544   0.352  
#> Efficacy boundary (p)         0.0001  0.0060  0.0231 
#> Information                   10.66   21.31   31.97  
```
