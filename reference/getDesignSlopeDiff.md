# Group Sequential Design for Two-Sample Slope Difference

Obtains the power given sample size or obtains the sample size given
power for a group sequential design for two-sample slope difference.

## Usage

``` r
getDesignSlopeDiff(
  beta = NA_real_,
  n = NA_real_,
  slopeDiffH0 = 0,
  slopeDiff = 0.5,
  stDev = 1,
  stDevCovariate = 1,
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
  futilitySlopeDiff = NA_real_,
  typeBetaSpending = "none",
  parameterBetaSpending = NA_real_,
  userBetaSpending = NA_real_,
  spendingTime = NA_real_
)
```

## Arguments

- beta:

  The type II error.

- n:

  The total sample size.

- slopeDiffH0:

  The slope difference under the null hypothesis. Defaults to 0.

- slopeDiff:

  The slope difference under the alternative hypothesis.

- stDev:

  The standard deviation of the residual.

- stDevCovariate:

  The standard deviation of the covariate.

- allocationRatioPlanned:

  Allocation ratio for the active treatment versus control. Defaults to
  1 for equal randomization.

- normalApproximation:

  The type of computation of the p-values. If `TRUE`, the variance is
  assumed to be known, otherwise the calculations are performed with the
  t distribution. The exact calculation using the t distribution is only
  implemented for the fixed design.

- rounding:

  Whether to round up sample size. Defaults to 1 for sample size
  rounding.

- kMax:

  The maximum number of stages.

- informationRates:

  The information rates. Fixed prior to the trial. Defaults to
  `(1:kMax) / kMax` if left unspecified.

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

- futilitySlopeDiff:

  The futility bounds on the slope difference scale.

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

An S3 class `designSlopeDiff` object with three components:

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

  - `expectedNumberOfSubjectsH1`: The expected number of subjects under
    H1.

  - `expectedNumberOfSubjectsH0`: The expected number of subjects under
    H0.

  - `slopeDiffH0`: The slope difference under the null hypothesis.

  - `slopeDiff`: The slope difference under the alternative hypothesis.

  - `stDev`: The standard deviation of the residual.

  - `stDevCovariate`: The standard deviation of the covariate.

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

- `settings`: A list containing the following input parameters:

  - `typeAlphaSpending`: The type of alpha spending.

  - `parameterAlphaSpending`: The parameter value for alpha spending.

  - `userAlphaSpending`: The user defined alpha spending.

  - `typeBetaSpending`: The type of beta spending.

  - `parameterBetaSpending`: The parameter value for beta spending.

  - `userBetaSpending`: The user defined beta spending.

  - `spendingTime`: The error spending time at each analysis.

  - `allocationRatioPlanned`: Allocation ratio for the active treatment
    versus control.

  - `normalApproximation`: The type of computation of the p-values. If
    `TRUE`, the variance is assumed to be known, otherwise the
    calculations are performed with the t distribution.

  - `rounding`: Whether to round up sample size.

## Details

We assume a simple linear regression of the form \$\$y\_{g,i} =
\alpha_g + \beta_g x\_{g,i} + \epsilon\_{g,i}\$\$ for treatment group
\\g\\, where \\\epsilon\_{g,i}\\ is the residual error for subject \\i\\
in group \\g\\, which is assumed to be normally distributed with mean 0
and standard deviation \\\sigma\_\epsilon\\. The covariate \\x\_{g,i}\\
is assumed to be normally distributed with mean 0 and standard deviation
\\\sigma_x\\. Since \$\$\hat{\beta}\_g = \frac{\sum\_{i=1}^{n_g}
(x\_{g,i}-\bar{x}\_g) y\_{g,i}}
{\sum\_{i=1}^{n_g}(x\_{g,i}-\bar{x}\_g)^2}\$\$ where \\n_g\\ is the
sample size for group \\g\\, it follows that \$\$\hat{\beta}\_g \sim
N(\beta_g,
\frac{\sigma\_\epsilon^2}{\sum\_{i=1}^{n_g}(x\_{g,i}-\bar{x}\_g)^2}).\$\$
The slope difference is defined as \$\$\hat{\beta}\_1 -
\hat{\beta}\_2\$\$ where \\\hat{\beta}\_1\\ and \\\hat{\beta}\_2\\ are
the estimated slopes for treatment groups 1 and 2, respectively. Since
the variance of \\\hat{\beta}\_g\\ is \$\$\frac{\sigma\_\epsilon^2}{n_g
\sigma_x^2}\$\$ we have \$\$\hat{\beta}\_1 - \hat{\beta}\_2 \sim
N(\beta_1 - \beta_2, (n_1^{-1} +
n_2^{-1})\sigma\_\epsilon^2/\sigma_x^2)\$\$ which can be used to
calculate the power and sample size for the group sequential design.

## Author

Kaifeng Lu, <kaifenglu@gmail.com>

## Examples

``` r
(design1 <- getDesignSlopeDiff(
  beta = 0.1, n = NA, slopeDiff = -0.5,
  stDev = 10, stDevCovariate = 6,
  normalApproximation = FALSE, alpha = 0.025))
#>                                                                        
#> Fixed design for two-sample slope difference                           
#> Slope difference under H0: 0, slope difference under H1: -0.5          
#> Standard deviation of residual: 10, standard deviation of covariate: 6 
#> Overall power: 0.9, overall alpha (1-sided): 0.025                     
#> Drift parameter: 3.248, inflation factor: 1                            
#> Information: 42.21                                                     
#> Number of subjects: 469                                                
#> Allocation ratio: 1                                                    
#>                                                                        
#>                                      
#> Efficacy boundary (t)          1.965 
#> Efficacy boundary (slope diff) -0.302
#> Efficacy boundary (p)          0.0250
```
