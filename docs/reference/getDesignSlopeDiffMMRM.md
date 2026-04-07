# Group Sequential Design for Two-Sample Slope Difference From the MMRM Model

Obtains the power given sample size or obtains the sample size given
power for two-sample slope difference from the growth curve MMRM model.

## Usage

``` r
getDesignSlopeDiffMMRM(
  beta = NA_real_,
  slopeDiffH0 = 0,
  slopeDiff = 0.5,
  stDev = 1,
  stDevIntercept = 1,
  stDevSlope = 1,
  corrInterceptSlope = 0.5,
  w = NA_real_,
  N = NA_real_,
  accrualTime = 0,
  accrualIntensity = NA_real_,
  piecewiseSurvivalTime = 0,
  gamma1 = 0,
  gamma2 = 0,
  accrualDuration = NA_real_,
  followupTime = NA_real_,
  fixedFollowup = FALSE,
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
  typeBetaSpending = "none",
  parameterBetaSpending = NA_real_,
  userBetaSpending = NA_real_,
  spendingTime = NA_real_
)
```

## Arguments

- beta:

  The type II error.

- slopeDiffH0:

  The slope difference under the null hypothesis. Defaults to 0.

- slopeDiff:

  The slope difference under the alternative hypothesis.

- stDev:

  The standard deviation of the residual.

- stDevIntercept:

  The standard deviation of the random intercept.

- stDevSlope:

  The standard deviation of the random slope.

- corrInterceptSlope:

  The correlation between the random intercept and random slope.

- w:

  The number of time units (e.g. weeks) per measurement visit in a
  period. In general, visits are more frequent in the beginning of the
  study and less frequent towards the end.

- N:

  The number of measurement visits in a period. For example,
  `w = c(8, 16)` and `N = c(2, Inf)` means that the response variable
  will be collected at baseline, week 8, week 16, and every 16 weeks
  thereafter.

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

- followupTime:

  Follow-up time for the last enrolled subject.

- fixedFollowup:

  Whether a fixed follow-up design is used. Defaults to `FALSE` for
  variable follow-up.

- allocationRatioPlanned:

  Allocation ratio for the active treatment versus control. Defaults to
  1 for equal randomization.

- normalApproximation:

  The type of computation of the p-values. If `TRUE`, the variance is
  assumed to be known, otherwise the calculations are performed with the
  t distribution. The degrees of freedom for the t-distribution for
  testing the slope difference is calculated using the containment
  method, and is equal to the total number of observations minus two
  times the total number of subjects. The exact calculation using the t
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

An S3 class `designSlopeDiffMMRM` object with three components:

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

  - `slopeDiffH0`: The slope difference under H0.

  - `slopeDiff`: The slope difference under H1.

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

  - `efficacySlopeDiff`: The efficacy boundaries on the slope difference
    scale.

  - `futilitySlopeDiff`: The futility boundaries on the slope difference
    scale.

  - `numberOfSubjects`: The number of subjects.

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

  - `w`: The number of time units per measurement visit in a period.

  - `N`: The number of measurement visits in a period.

  - `stdDev`: The standard deviation of the residual.

  - `G`: The covariance matrix for the random intercept and random
    slope.

  - `normalApproximation`: The type of computation of the p-values. If
    `TRUE`, the variance is assumed to be known, otherwise the
    calculations are performed with the t distribution.

  - `rounding`: Whether to round up sample size.

## Details

We use the following random-effects model to compare two slopes:
\$\$y\_{ij} = \alpha + (\beta + \gamma x_i) t_j + a_i + b_i t_j +
e\_{ij}\$\$ where

- \\\alpha\\: overall intercept common across treatment groups due to
  randomization

- \\\beta\\: slope for the control group

- \\\gamma\\: difference in slopes between the active treatment and
  control groups

- \\x_i\\: treatment indicator for subject \\i\\, 1 for the active
  treatment and 0 for the control

- \\t_j\\: time point \\j\\ for repeated measurements, \\t_1 = 0 \< t_2
  \< \ldots \< t_k\\

- \\(a_i, b_i)\\: random intercept and random slope for subject \\i\\,
  \\Var(a_i) = \sigma_a^2\\, \\Var(b_i) = \sigma_b^2\\, \\Corr(a_i, b_i)
  = \rho\\

- \\e\_{ij}\\: within-subject residual with variance \\\sigma_e^2\\

By accounting for randomization, we improve the efficiency for
estimating the difference in slopes. The model also accommodates
unequally spaced time points and missing data. Specifically, given a
calendar time \\\tau\\ for an interim or final analysis, let \\k\\ be
the number of scheduled time points up to and including \\\tau\\,
subject to the follow-up duration for fixed follow-up designs. Let the
observed time points be \\t_1, t_2, \ldots, t_k\\, where \\t_1 = 0\\
denotes baseline.

For a subject in treatment group \\g\\ with observed data pattern \\j\\,
the design matrix for the fixed effects \\(\alpha, \beta, \gamma)'\\ is
given by \$\$X\_{g,j} = (\bm{1}\_j, \vec{t}\_j, I(g=1)\vec{t}\_j)\$\$
where \\\bm{1}\_j\\ is a \\j\\-vector of ones, and \\\vec{t}\_j =
(t_1,\ldots,t_j)'\\ is the column vector of observed time points. The
design matrix for the random effects \\(a_i, b_i)'\\ is \$\$Z_j =
(\bm{1}\_j, \vec{t}\_j)\$\$ The variance-covariance matrix of the random
effects is \$\$D = \left(\begin{array}{cc} \sigma_a^2 & \rho \sigma_a
\sigma_b \\ \rho \sigma_a \sigma_b & \sigma_b^2 \end{array}\right)\$\$
Therefore, the variance-covariance matrix for the observed data for the
subject is \$\$V\_{j} = Z_j D Z_j' + \sigma_e^2 I_j\$\$ where \\I_j\\ is
the \\j\times j\\ identity matrix. Let \\\pi_g\\ denote the proportion
of subjects in group \\g\\. The information matrix for the fixed effects
is \$\$I = nJ\$\$ where \$\$J = \sum\_{g=1}^{2} \pi_g \sum\_{j=1}^{k}
p\_{g,j} X\_{g,j}' V_j^{-1} X\_{g,j}\$\$ and \\p\_{g,j}\\ is the
proportion of subjects in group \\g\\ with observed data pattern \\j\\.

The variance of the estimator for the slope difference \\\hat{\gamma}\\
is given by \$\$\text{Var}(\hat{\gamma}) = \frac{1}{n} J^{-1}\[3,3\]\$\$
which can be used to calculate the power and sample size for the group
sequential design to detect a slope difference.

## References

Daniel O. Scharfstein, Anastasios A. Tsiatis, and James M. Robins.
Semiparametric efficiency and its implication on the design and analysis
of group-sequential studies. Journal of the American Statistical
Association 1997; 92:1342-1350.

## Author

Kaifeng Lu, <kaifenglu@gmail.com>

## Examples

``` r
(design1 <- getDesignSlopeDiffMMRM(
  beta = 0.2, slopeDiff = log(1.15)/52,
  stDev = sqrt(.182),
  stDevIntercept = sqrt(.238960),
  stDevSlope = sqrt(.000057),
  corrInterceptSlope = .003688/sqrt(.238960*.000057),
  w = 8,
  N = 10000,
  accrualIntensity = 15,
  gamma1 = 1/(4.48*52),
  gamma2 = 1/(4.48*52),
  accrualDuration = NA,
  followupTime = 8,
  alpha = 0.025))
#>                                                                              
#> Fixed design for two-sample slope difference from the MMRM model             
#> Slope difference under H0: 0, slope difference under H1: 0.003               
#> Standard deviation (SD) of within-subject residual: 0.427                    
#> SD of random intercept: 0.489, SD of random slope: 0.008, correlation: 0.999 
#> Overall power: 0.8004, overall alpha (1-sided): 0.025                        
#> Drift parameter: 2.803, inflation factor: 1                                  
#> Information: 1087517.15                                                      
#> Number of subjects: 1013                                                     
#> Study duration: 75.5                                                         
#> Accrual duration: 67.5, follow-up duration: 8, fixed follow-up: FALSE        
#> Allocation ratio: 1                                                          
#>                                                                              
#>                                      
#> Efficacy boundary (Z)          1.960 
#> Efficacy boundary (slope diff) 0.002 
#> Efficacy boundary (p)          0.0250
```
