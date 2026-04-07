# Group Sequential Design for McNemar's Test for Paired Proportions

Obtains the power given sample size or obtains the sample size given
power for a group sequential design for McNemar's test for paired
proportions.

## Usage

``` r
getDesignPairedPropMcNemar(
  beta = NA_real_,
  n = NA_real_,
  pDiscordant = NA_real_,
  riskDiff = NA_real_,
  nullVariance = TRUE,
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

- pDiscordant:

  The proportion of discordant pairs (xi = pi01 + pi10).

- riskDiff:

  The risk difference between the active and control treatments (delta =
  pi_t - pi_c = pi01 - pi10)

- nullVariance:

  Whether to use the variance under the null or the variance under the
  alternative.

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
  for `"WT"`, \\rho\\ for `"sfKD"`, and \\\gamma\\ for `"sfHSD"`.

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

An S3 class `designPairedPropMcNemar` object with three components:

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

  - `pDiscordant`: The proportion of discordant pairs (xi = pi01 +
    pi10).

  - `riskDiff`: The risk difference between the active and control
    treatments (delta = pi_t - pi_c = pi01 - pi10)

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

  - `efficacyRiskDiff`: The efficacy boundaries on the risk difference
    scale.

  - `futilityRiskDiff`: The futility boundaries on the risk difference
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

  - `varianceRatio`: The ratio of the variance under H0 to the variance
    under H1.

  - `rounding`: Whether to round up sample size.

## Details

Consider a group sequential design for McNemar's test for paired
proportions. The table below shows joint probabilities for each cell
(\\\pi\_{ij}\\ where \\i\\ is for control group and \\j\\ is for
experimental group), with marginal totals.

|  |  |  |  |
|----|----|----|----|
|  | Experimental: No Response | Experimental: Response | Row Total |
| Control: No Response | \\\pi\_{00}\\ | \\\pi\_{01}\\ | \\1-\pi_c\\ |
| Control: Response | \\\pi\_{10}\\ | \\\pi\_{11}\\ | \\\pi_c\\ |
| Column Total | \\1-\pi_t\\ | \\\pi_t\\ | 1 |

The parameters \\\pi\_{01}\\ and \\\pi\_{10}\\ are the discordant pairs
(i.e., \\\pi\_{01} + \pi\_{10} = \xi\\) and the risk difference is
\\\pi\_{01} - \pi\_{10} = \delta\\. The parameter \\\pi_t\\ is the
proportion of experimental group response, and \\\pi_c\\ is the
proportion of control group response. The parameter of interest is
\$\$\theta = \pi_t - \pi_c = \pi\_{01} - \pi\_{10} = \delta\$\$ The
variance of \\\hat{\theta}\\ can be obtained from the multinomial
distribution as follows: \$\$Var(\hat{\theta}) = \frac{1}{n} \\
\pi\_{01}(1-\pi\_{01}) + \pi\_{10}(1-\pi\_{10}) + 2\pi\_{01}\pi\_{10}
\\\$\$ which can be simplified to \$\$Var(\hat{\theta}) = \frac{1}{n}
(\xi - \delta^2)\$\$ Here, \\n\\ is the total number of treatment pairs.
This is the unconditional variance, which is used for the overall
design.

## Author

Kaifeng Lu, <kaifenglu@gmail.com>

## Examples

``` r
# Example 1: fixed design
(design1 <- getDesignPairedPropMcNemar(
  beta = 0.1, n = NA, pDiscordant = 0.16, riskDiff = 0.1,
  alpha = 0.025))
#>                                                            
#> Fixed design for McNemar's test                            
#> Proportion of discordant pairs: 0.16, risk difference: 0.1 
#> Overall power: 0.9001, overall alpha (1-sided): 0.025      
#> Drift parameter: 3.307, inflation factor: 1                
#> Information: 1093.33                                       
#> Number of subjects: 164                                    
#> Variance of standardized test statistic: under H0          
#>                                                            
#>                                     
#> Efficacy boundary (Z)         1.960 
#> Efficacy boundary (risk diff) 0.0612
#> Efficacy boundary (p)         0.0250

# Example 2: group sequential design
(design2 <- getDesignPairedPropMcNemar(
  beta = 0.1, n = NA, pDiscordant = 0.16, riskDiff = 0.1,
  alpha = 0.025, kMax = 3, typeAlphaSpending = "sfOF"))
#>                                                                                        
#> Group-sequential design with 3 stages for McNemar's test                               
#> Proportion of discordant pairs: 0.16, risk difference: 0.1                             
#> Overall power: 0.9016, overall alpha (1-sided): 0.025                                  
#> Drift parameter: 3.337, inflation factor: 1.013                                        
#> Maximum information: 1113.33, expected under H1: 897.63, expected under H0: 1111.06    
#> Maximum # subjects: 167, expected under H1: 134.6, expected under H0: 166.7            
#> Variance of standardized test statistic: under H0                                      
#> Alpha spending: Lan-DeMets O'Brien-Fleming, Alpha spending: Lan-DeMets O'Brien-Fleming 
#>                                                                                        
#>                               Stage 1 Stage 2 Stage 3
#> Information rate              0.335   0.665   1.000  
#> Efficacy boundary (Z)         3.698   2.516   1.993  
#> Cumulative rejection          0.0296  0.5487  0.9016 
#> Cumulative alpha spent        0.0001  0.0060  0.0250 
#> Number of subjects            56.0    111.0   167.0  
#> Efficacy boundary (risk diff) 0.198   0.096   0.062  
#> Efficacy boundary (p)         0.0001  0.0059  0.0231 
#> Information                   373.33  740.00  1113.33
```
