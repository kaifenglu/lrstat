# Adaptive Multi-Arm Multi-Stage Design

Calculates the conditional power for specified incremental information,
given the interim results, parameter value, data-dependent changes in
treatment selection, the error spending function, and the number and
spacing of interim looks. Conversely, calculates the incremental
information required to attain a specified conditional power, given the
interim results, parameter value, data-dependent changes in treatment
selection, the error spending function, and the number and spacing of
interim looks.

## Usage

``` r
adaptDesign_mams(
  betaNew = NA_real_,
  INew = NA_real_,
  M = NA_integer_,
  r = 1,
  corr_known = TRUE,
  L = NA_integer_,
  zL = NA_real_,
  theta = NA_real_,
  IMax = NA_real_,
  kMax = NA_integer_,
  informationRates = NA_real_,
  efficacyStopping = NA_integer_,
  futilityStopping = NA_integer_,
  criticalValues = NULL,
  alpha = 0.025,
  typeAlphaSpending = "sfOF",
  parameterAlphaSpending = NA_real_,
  userAlphaSpending = NA_real_,
  futilityBounds = NULL,
  futilityCP = NULL,
  futilityTheta = NULL,
  typeBetaSpending = "none",
  parameterBetaSpending = NA_real_,
  spendingTime = NA_real_,
  MullerSchafer = FALSE,
  MNew = NA_integer_,
  selected = NA_integer_,
  rNew = 1,
  kNew = NA_integer_,
  informationRatesNew = NA_real_,
  efficacyStoppingNew = NA_integer_,
  futilityStoppingNew = NA_integer_,
  typeAlphaSpendingNew = "sfOF",
  parameterAlphaSpendingNew = NA_real_,
  futilityBoundsInt = NULL,
  futilityCPInt = NULL,
  futilityThetaInt = NULL,
  typeBetaSpendingNew = "none",
  parameterBetaSpendingNew = NA_real_,
  userBetaSpendingNew = NA_real_,
  spendingTimeNew = NA_real_
)
```

## Arguments

- betaNew:

  The type II error for the secondary trial.

- INew:

  The maximum information for any active arm versus the common control
  in the secondary trial. Either `betaNew` or `INew` should be provided,
  while the other must be missing.

- M:

  Number of active treatment arms in the primary trial.

- r:

  Randomization ratio of each active arm to the common control in the
  primary trial.

- corr_known:

  Logical. If `TRUE`, the correlation between Wald statistics is derived
  from the randomization ratio \\r\\ as \\r / (r + 1)\\. If `FALSE`, a
  conservative correlation of 0 is assumed.

- L:

  The interim adaptation look of the primary trial.

- zL:

  The z-test statistics at the interim adaptation look of the primary
  trial.

- theta:

  A vector of length \\M\\ representing the assumed treatment effects
  for each active arm versus the common control. The global null is
  \\\theta_i = 0\\ for all \\i\\, and alternatives are one-sided:
  \\\theta_i \> 0\\ for at least one \\i = 1, \ldots, M\\.

- IMax:

  Maximum information for any active arm versus the common control for
  the primary trial. Must be provided.

- kMax:

  The maximum number of stages of the primary trial.

- informationRates:

  The information rates of the primary trial.

- efficacyStopping:

  Indicators of whether efficacy stopping is allowed at each stage of
  the primary trial. Defaults to `TRUE` if left unspecified.

- futilityStopping:

  Indicators of whether futility stopping is allowed at each stage of
  the primary trial.

- criticalValues:

  The matrix of by-level upper boundaries on the max z-test statistic
  scale for efficacy stopping for the primary trial. The first column is
  for level `M`, the second column is for level `M - 1`, and so on, with
  the last column for level 1. If left unspecified, the critical values
  will be computed based on the specified alpha spending function.

- alpha:

  The significance level of the primary trial. Defaults to 0.025.

- typeAlphaSpending:

  The type of alpha spending for the primary trial. One of the
  following: `"OF"` for O'Brien-Fleming boundaries, `"P"` for Pocock
  boundaries, `"WT"` for Wang & Tsiatis boundaries, `"sfOF"` for
  O'Brien-Fleming type spending function, `"sfP"` for Pocock type
  spending function, `"sfKD"` for Kim & DeMets spending function,
  `"sfHSD"` for Hwang, Shi & DeCani spending function, `"user"` for user
  defined spending, and `"none"` for no early efficacy stopping.
  Defaults to `"sfOF"`.

- parameterAlphaSpending:

  The parameter value of alpha spending for the primary trial.
  Corresponds to \\\Delta\\ for `"WT"`, \\\rho\\ for `"sfKD"`, and
  \\\gamma\\ for `"sfHSD"`.

- userAlphaSpending:

  The user-defined alpha spending for the primary trial. Represents the
  cumulative alpha spent up to each stage.

- futilityBounds:

  The futility boundaries on the max-z statistic scale for the primary
  trial. Defaults to `rep(-8, kMax-1)` if left unspecified.

- futilityCP:

  The conditional power-based futility bounds for the primary trial.

- futilityTheta:

  The parameter value-based futility bounds for the primary trial.

- typeBetaSpending:

  The type of beta spending for the primary trial. One of the following:
  `"sfOF"` for O'Brien-Fleming type spending function, `"sfP"` for
  Pocock type spending function, `"sfKD"` for Kim & DeMets spending
  function, `"sfHSD"` for Hwang, Shi & DeCani spending function, and
  `"none"` for no early futility stopping. Defaults to `"none"`.

- parameterBetaSpending:

  The parameter value of beta spending for the primary trial.
  Corresponds to \\\rho\\ for `"sfKD"`, and \\\gamma\\ for `"sfHSD"`.

- spendingTime:

  The error spending time of the primary trial. Defaults to missing, in
  which case it is assumed to be the same as `informationRates`.

- MullerSchafer:

  Whether to use the Muller and Schafer (2001) method for trial
  adaptation.

- MNew:

  Number of active treatment arms in the secondary trial.

- selected:

  The indices of the selected active treatment arms for the secondary
  trial.

- rNew:

  Randomization ratio of each active arm to the common control in the
  secondary trial.

- kNew:

  The number of looks of the secondary trial.

- informationRatesNew:

  The spacing of looks of the secondary trial.

- efficacyStoppingNew:

  The indicators of whether efficacy stopping is allowed at each look of
  the secondary trial. Defaults to `TRUE` if left unspecified.

- futilityStoppingNew:

  The indicators of whether futility stopping is allowed at each look of
  the secondary trial. Defaults to `TRUE` if left unspecified.

- typeAlphaSpendingNew:

  The type of alpha spending for the secondary trial. One of the
  following: `"OF"` for O'Brien-Fleming boundaries, `"sfOF"` for
  O'Brien-Fleming type spending function, `"sfP"` for Pocock type
  spending function, `"sfKD"` for Kim & DeMets spending function,
  `"sfHSD"` for Hwang, Shi & DeCani spending function, and `"none"` for
  no early efficacy stopping. Defaults to `"sfOF"`.

- parameterAlphaSpendingNew:

  The parameter value of alpha spending for the secondary trial.
  Corresponds to \\\rho\\ for `"sfKD"`, and \\\gamma\\ for `"sfHSD"`.

- futilityBoundsInt:

  The futility boundaries on the max-z statistic scale for new stages of
  the integrated trial.

- futilityCPInt:

  The conditional power-based futility bounds for new stages of the
  integrated trial.

- futilityThetaInt:

  The parameter value-based futility bounds for the new stages of the
  integrated trial.

- typeBetaSpendingNew:

  The type of beta spending for the secondary trial. One of the
  following: `"sfOF"` for O'Brien-Fleming type spending function,
  `"sfP"` for Pocock type spending function, `"sfKD"` for Kim & DeMets
  spending function, `"sfHSD"` for Hwang, Shi & DeCani spending
  function, `"user"` for user defined spending, and `"none"` for no
  early futility stopping. Defaults to `"none"`.

- parameterBetaSpendingNew:

  The parameter value of beta spending for the secondary trial.
  Corresponds to \\\rho\\ for `"sfKD"`, and \\\gamma\\ for `"sfHSD"`.

- userBetaSpendingNew:

  The user-defined cumulative beta spending. Represents the cumulative
  beta spent up to each stage of the secondary trial.

- spendingTimeNew:

  The error spending time of the secondary trial. Defaults to missing,
  in which case it is assumed to be the same as `informationRatesNew`.

## Value

An `adaptDesign_mams` object with three list components:

- `primaryTrial`: A list of selected information for the primary trial,
  including `M`, `r`, `corr_known`, `L`, `zL`, `theta`,
  `maxInformation`, `kMax`, `informationRates`, `efficacyBounds`,
  `futilityBounds`, `information`, `alpha`, `conditionalAlpha`,
  `conditionalPower`, `MullerSchafer`, and `byLevelBounds`, where
  `byLevelBounds` is a data frame with columns `level`, `stage`, and
  `efficacyBounds`, representing the efficacy bounds for each
  combination of the number of active arms and the stage of analysis in
  the primary trial.

- `secondaryTrial`: A list of selected information for the secondary
  trial, including `overallReject`, `alpha`, `M`, `r`, `selected`,
  `corr_known`, `kMax`, `maxInformation`, `informationRates`,
  `cumulativeRejection`, `cumulativeAlphaSpent`, `information`,
  `typeAlphaSpending`, `parameterAlphaSpending`, `typeBetaSpending`,
  `parameterBetaSpending`, `userBetaSpending`, `spendingTime`, and
  `byHypothesisBounds`, where `byHypothesisBounds` is a data frame with
  columns `hypothesis`, `stage`, `efficacyBounds`, and `futilityBounds`,
  representing the efficacy and futility bounds for each hypothesis and
  each stage of analysis in the secondary trial.

- `integratedTrial`: A list of selected information for the integrated
  trial, including `M`, `r`, `corr_known`, `MNew`, `rNew`, `selected`,
  `L`, `zL`, `theta`, `maxInformation`, `kMax`, `informationRates`,
  `efficacyBounds`, `futilityBounds`, `information`, and
  `byIntersectionBounds`, where `byIntersectionBounds` is a data frame
  with columns `intersectionHypothesis`, `stage`, and `efficacyBounds`,
  representing the efficacy bounds for each intersection hypothesis and
  each stage of analysis in the integrated trial.

## References

Ping Gao, Yingqiu Li. Adaptive multiple comparison sequential design
(AMCSD) for clinical trials. Journal of Biopharmaceutical Statistics,
2024, 34(3), 424-440.

## See also

[`getDesign_mams`](https://kaifenglu.github.io/lrstat/reference/getDesign_mams.md)

## Author

Kaifeng Lu, <kaifenglu@gmail.com>

## Examples

``` r
# Two active treatment arms are compared with a common control in a
# two-look time-to-event design using O'Brien–Fleming–type alpha spending.
# Suppose each active arm has a true hazard ratio of 0.75 versus control,
# and the total number of events across all three arms at the final analysis
# is 486. This corresponds to approximately 324 events for each active arm
# versus the common control. Under these assumptions, the trial has about
# 80% power to detect the treatment effect in at least one active arm.

(des1 <- getDesign_mams(
  IMax = 324 / 4, theta = c(-log(0.75), -log(0.75)),
  M = 2, r = 1, kMax = 2, informationRates = c(1/2, 1),
  alpha = 0.025, typeAlphaSpending = "OF"))
#>                                                                            
#> Multi-arm multi-stage design                                               
#> Overall power: 0.7994, overall alpha (1-sided): 0.025                      
#> Number of active arms: 2                                                   
#> Randomization ratio of each active vs. control: 1                          
#> Using correlation for critical value calculation: TRUE                     
#> Max information for pairwise comparion: 81                                 
#> Number of looks: 2                                                         
#> Expected information under H1: 74.53, expected information under H0: 80.93 
#> Alpha spending: O'Brien-Fleming, beta spending: None                       
#>                                                                            
#>                           Stage 1 Stage 2
#> Information rate          0.500   1.000  
#> Efficacy boundary (Z)     3.142   2.222  
#> Cumulative rejection      0.1598  0.7994 
#> Cumulative alpha spent    0.0016  0.0250 
#> Efficacy boundary (theta) 0.494   0.247  
#> Efficacy boundary (p)     0.0008  0.0132 
#> Information               40.50   81.00  
#> 
#> By level critical values
#>   Level Stage Boundary (Z)
#> 1     2     1        3.142
#> 2     2     2        2.222
#> 3     1     1        2.797
#> 4     1     2        1.977

# Now assume that, at the interim analysis, the observed hazard ratios for
# the two active arms versus control are 0.91 and 0.78, respectively. Using
# the rule “drop any arm with an observed hazard ratio > 0.9”, arm 1 is
# dropped. We then aim to achieve 80% conditional power to detect a hazard
# ratio of 0.78 for the remaining arm at the final look. The analysis below
# indicates that the required total number of events for arm 2 versus control
# at the final analysis should be increased from 324 to 535.

(des2 <- adaptDesign_mams(
  betaNew = 0.2, M = 2, r = 1, corr_known = FALSE,
  L = 1, zL = c(-log(0.91), -log(0.78)) * sqrt(324 / 4 / 2),
  theta = c(-log(0.91), -log(0.78)),
  IMax = 324 / 4, kMax = 2, informationRates = c(1/2, 1),
  alpha = 0.025, typeAlphaSpending = "OF",
  MNew = 1, selected = 2, rNew = 1))
#>                                                            
#> Primary trial:                                             
#> Multi-arm multi-stage design                               
#> Number of active arms: 2                                   
#> Randomization ratio of each active vs. control: 1          
#> Using correlation for critical value calculation: FALSE    
#> Max information for pairwise comparion: 81                 
#> Number of looks: 2                                         
#> Interim adaptation look: 1, z-statistic value: 0.6, 1.581  
#> theta: 0.094, 0.248                                        
#> Conditional type I error: 0.0597, conditional power: 0.496 
#> Muller & Schafer method for secondary trial: FALSE         
#>                                                            
#>                       Stage 1 Stage 2
#> Information rate      0.500   1.000  
#> Efficacy boundary (Z) 3.179   2.248  
#> Information           40.50   81.00  
#> 
#> By level critical values for primary trial
#>   Level Stage Boundary (Z)
#> 1     2     1        3.179
#> 2     2     2        2.248
#> 3     1     1        2.797
#> 4     1     2        1.977
#>                                                                  
#> Secondary trial:                                                 
#> Multi-arm multi-stage design                                     
#> Number of selected active arms: 1, selected active arms: 2       
#> Randomization ratio of each active vs. control: 1                
#> Maximum information: 93.22                                       
#> Overall power: 0.8, overall significance level (1-sided): 0.0597 
#>                                                                  
#>                        Stage 1
#> Information rate       1.000  
#> Cumulative rejection   0.8000 
#> Cumulative alpha spent 0.0597 
#> Information            93.22  
#> 
#> By hypothesis critical values for secondary trial
#>   Hypothesis Stage Efficacy boundary (Z) Futility boundary (Z)
#> 1          2     1                 1.557                 1.557
#>                                                            
#> Integrated trial:                                          
#> Adaptive multi-arm multi-stage design                      
#> Number of active arms before adaptation: 2                 
#> Number of selected active arms: 1, selected active arms: 2 
#> Total number of looks: 2                                   
#> Interim adaptation look: 1, z-statistic value: 0.6, 1.581  
#>                                                            
#>                     Stage 1 Stage 2
#> Information rate    0.303   1.000  
#> Efficacy bounds (Z) 3.179   2.170  
#> Information         40.50   133.72 
```
