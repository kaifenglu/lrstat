# Group Sequential Design for Two-Sample Wilcoxon Test

Obtains the power given sample size or obtains the sample size given
power for a group sequential design for two-sample Wilcoxon test.

## Usage

``` r
getDesignWilcoxon(
  beta = NA_real_,
  n = NA_real_,
  pLarger = 0.6,
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

- pLarger:

  The probability that a randomly chosen sample from the treatment group
  is larger than a randomly chosen sample from the control group under
  the alternative hypothesis.

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

An S3 class `designWilcoxon` object with three components:

- `overallResults`: A data frame containing the following variables:

  - `overallReject`: The overall rejection probability.

  - `alpha`: The overall significance level.

  - `attainedAlpha`: The attained significance level, which is different
    from the overall significance level in the presence of futility
    stopping..

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

  - `pLarger`: The probability that a randomly chosen sample from the
    treatment group is larger than a randomly chosen sample from the
    control group under the alternative hypothesis.

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

  - `efficacyPLarger`: The efficacy boundaries on the proportion of
    pairs of samples from the two treatment groups with the sample from
    the treatment group greater than that from the control group.

  - `futilityPLarger`: The futility boundaries on the proportion of
    pairs of samples from the two treatment groups with the sample from
    the treatment group greater than that from the control group.

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

  - `rounding`: Whether to round up sample size.

## Details

The Mann-Whitney U test is a non-parametric test for the difference in
location between two independent groups. It is also known as the
Wilcoxon rank-sum test. The test is based on the ranks of the data
rather than the actual values, making it robust to outliers and
non-normal distributions. The test statistic is the number of times a
randomly chosen sample from the treatment group is larger than a
randomly chosen sample from the control group, i.e., \$\$W\_{XY} =
\sum\_{i=1}^{n_1}\sum\_{j=1}^{n_2} I(X_i \> Y_j)\$\$ where \\X_i\\ and
\\Y_j\\ are the samples from the treatment and control groups,
respectively. The test is often used when the data do not meet the
assumptions of the t-test, such as non-normality or unequal variances.
The test is also applicable to ordinal data. The test is one-sided,
meaning that it only tests whether the treatment group is larger than
the control group. Asymptotically, \$\$\frac{W\_{XY} - n_1
n_2/2}{\sqrt{n_1 n_2 (n+1)/12}} \sim N(0,1) \quad \text{under} H_0\$\$
where \\n_1\\ and \\n_2\\ are the sample sizes of the treatment and
control groups, respectively, and \\n=n_1+n_2\\. Let \\\theta = P(X \>
Y)\\, and \\\hat{\theta} = \frac{1}{nm}W\_{XY}\\. It follows that
\$\$\sqrt{n}(\hat{\theta} - 1/2) \sim N\left(0, \frac{1}{12
r(1-r)}\right) \quad \text{under} H_0\$\$ where \\r = n_1/(n_1+n_2)\\ is
the randomization probability for the active treatment group. This
formulation allows for group sequential testing with futility stopping
and efficacy stopping.

## Author

Kaifeng Lu, <kaifenglu@gmail.com>

## Examples

``` r
# Example 1: fixed design
(design1 <- getDesignWilcoxon(
  beta = 0.1, n = NA,
  pLarger = pnorm((8 - 2)/sqrt(2*25^2)), alpha = 0.025))
#>                                                                              
#> Fixed design for two-sample Wilcoxon test                                    
#> Probability of observations in treatment larger than those in control: 0.567 
#> Overall power: 0.9002, overall alpha (1-sided): 0.025                        
#> Drift parameter: 3.243, inflation factor: 1                                  
#> Information: 2316                                                            
#> Number of subjects: 772                                                      
#> Allocation ratio: 1                                                          
#>                                                                              
#>                                   
#> Efficacy boundary (Z)       1.960 
#> Efficacy boundary (pLarger) 0.5407
#> Efficacy boundary (p)       0.0250

# Example 2: group sequential design
(design2 <- getDesignWilcoxon(
  beta = 0.1, n = NA,
  pLarger = pnorm((8 - 2)/sqrt(2*25^2)), alpha = 0.025,
  kMax = 3, typeAlphaSpending = "sfOF"))
#>                                                                                   
#> Group-sequential design with 3 stages for two-sample Wilcoxon test                
#> Probability of observations in treatment larger than those in control: 0.567      
#> Overall power: 0.9001, overall alpha (1-sided): 0.025                             
#> Drift parameter: 3.261, inflation factor: 1.012                                   
#> Maximum information: 2343, expected under H1: 1878.94, expected under H0: 2338.19 
#> Maximum # subjects: 781, expected under H1: 626.3, expected under H0: 779.4       
#> Allocation ratio: 1                                                               
#> Alpha spending: Lan-DeMets O'Brien-Fleming, beta spending: None                   
#>                                                                                   
#>                             Stage 1 Stage 2 Stage 3
#> Information rate            0.333   0.667   1.000  
#> Efficacy boundary (Z)       3.713   2.510   1.993  
#> Cumulative rejection        0.0335  0.5613  0.9001 
#> Cumulative alpha spent      0.0001  0.0061  0.0250 
#> Number of subjects          260.0   521.0   781.0  
#> Efficacy boundary (pLarger) 0.6329  0.5635  0.5412 
#> Efficacy boundary (p)       0.0001  0.0060  0.0231 
#> Information                 780.00  1563.00 2343.00
```
