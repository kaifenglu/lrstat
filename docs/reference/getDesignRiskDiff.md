# Group Sequential Design for Two-Sample Risk Difference

Obtains the power given sample size or obtains the sample size given
power for a group sequential design for two-sample risk difference.

## Usage

``` r
getDesignRiskDiff(
  beta = NA_real_,
  n = NA_real_,
  riskDiffH0 = 0,
  pi1 = NA_real_,
  pi2 = NA_real_,
  nullVariance = TRUE,
  allocationRatioPlanned = 1,
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

- n:

  The total sample size.

- riskDiffH0:

  The risk difference under the null hypothesis. Defaults to 0.

- pi1:

  The assumed probability for the active treatment group.

- pi2:

  The assumed probability for the control group.

- nullVariance:

  Whether to use the variance under the null or the empirical variance
  under the alternative.

- allocationRatioPlanned:

  Allocation ratio for the active treatment versus control. Defaults to
  1 for equal randomization.

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

An S3 class `designRiskDiff` object with three components:

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

  - `riskDiffH0`: The risk difference under the null hypothesis.

  - `pi1`: The assumed probability for the active treatment group.

  - `pi2`: The assumed probability for the control group.

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

  - `efficacyRiskDiff`: The efficacy boundaries on the risk difference
    scale.

  - `futilityRiskDiff`: The futility boundaries on the risk difference
    scale.

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

  - `numberOfSubjects`: The number of subjects.

- `settings`: A list containing the following input parameters:

  - `typeAlphaSpending`: The type of alpha spending.

  - `parameterAlphaSpending`: The parameter value for alpha spending.

  - `userAlphaSpending`: The user defined alpha spending.

  - `typeBetaSpending`: The type of beta spending.

  - `parameterBetaSpending`: The parameter value for beta spending.

  - `userBetaSpending`: The user defined beta spending.

  - `spendingTime`: The error spending time at each analysis.

  - `varianceRatio`: The ratio of the variance under H0 to the variance
    under H1.

  - `nullVariance`: Whether to use the variance under the null or the
    empirical variance under the alternative.

  - `allocationRatioPlanned`: Allocation ratio for the active treatment
    versus control.

  - `rounding`: Whether to round up sample size.

## Details

Consider a group sequential design for two-sample risk difference. The
parameter of interest is \$\$\theta = \pi_1 - \pi_2\$\$ where \\\pi_1\\
is the response probability for the active treatment group and \\\pi_2\\
is the response probability for the control group. The variance of
\\\hat{\theta}\\ can be obtained from the binomial distributions as
follows: \$\$Var(\hat{\theta}) = \frac{1}{n} \\
\frac{\pi_1(1-\pi_1)}{r} + \frac{\pi_2(1-\pi_2)}{1-r} \\\$\$ where \\n\\
is the total number of subjects and \\r\\ is the randomization
probability for the active treatment group. When `nullVariance = TRUE`,
the variance is computed under the null hypothesis. In this case, the
values of \\\pi_1\\ and \\\pi_2\\ in the variance formula are replaced
with their restricted maximum likelihood counterparts, subject to the
constraint \$\$\pi_1 - \pi_2 = \theta_0\$\$

## Author

Kaifeng Lu, <kaifenglu@gmail.com>

## Examples

``` r
(design1 <- getDesignRiskDiff(
  beta = 0.2, n = NA, pi1 = 0.1, pi2 = 0.15,
  kMax = 3, alpha = 0.025, typeAlphaSpending = "sfOF",
  nullVariance = FALSE))
#>                                                                                        
#> Group-sequential design with 3 stages for two-sample risk difference                   
#> Risk difference under H0: 0, proportion on treatment: 0.1, proportion on control: 0.15 
#> Overall power: 0.8002, overall alpha (1-sided): 0.025                                  
#> Drift parameter: 2.82, inflation factor: 1.013                                         
#> Maximum information: 3181.61, expected under H1: 2718.8, expected under H0: 3175.08    
#> Maximum # subjects: 1384, expected under H1: 1182.7, expected under H0: 1381.2         
#> Allocation ratio: 1, variance of standardized test statistic: under H1                 
#> Alpha spending: Lan-DeMets O'Brien-Fleming, beta spending: None                        
#>                                                                                        
#>                               Stage 1 Stage 2 Stage 3
#> Information rate              0.333   0.667   1.000  
#> Efficacy boundary (Z)         3.712   2.511   1.993  
#> Cumulative rejection          0.0186  0.4181  0.8002 
#> Cumulative alpha spent        0.0001  0.0061  0.0250 
#> Number of subjects            461.0   923.0   1384.0 
#> Efficacy boundary (risk diff) -0.114  -0.055  -0.035 
#> Efficacy boundary (p)         0.0001  0.0060  0.0231 
#> Information                   1059.77 2121.84 3181.61
```
