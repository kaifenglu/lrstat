# Power and Sample Size for Phase 2/3 Seamless Design

Computes either the maximum information and stopping boundaries for a
phase 2/3 seamless design, or the achieved power when the maximum
information and stopping boundaries are provided. Both efficacy and
futility stopping can be incorporated.

## Usage

``` r
getDesign_seamless(
  beta = NA_real_,
  IMax = NA_real_,
  theta = NA_real_,
  M = NA_integer_,
  r = 1,
  corr_known = TRUE,
  K = 1L,
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
  userBetaSpending = NA_real_,
  spendingTime = NA_real_
)
```

## Arguments

- beta:

  Type II error rate. Provide either `beta` or `IMax`; the other should
  be missing.

- IMax:

  Maximum information for any active arm versus the common control.
  Provide either `IMax` or `beta`; the other should be missing.

- theta:

  A vector of length \\M\\ representing the true treatment effects for
  each active arm versus the common control. The global null is
  \\\theta_i = 0\\ for all \\i\\, and alternatives are one-sided:
  \\\theta_i \> 0\\ for at least one \\i = 1, \ldots, M\\.

- M:

  Number of active treatment arms in Phase 2.

- r:

  Randomization ratio of each active arm to the common control in Phase
  2.

- corr_known:

  Logical. If `TRUE`, the correlation between Wald statistics in Phase 2
  is derived from the randomization ratio \\r\\ as \\r / (r + 1)\\. If
  `FALSE`, a conservative correlation of 0 is used.

- K:

  Number of sequential looks in Phase 3.

- informationRates:

  A numeric vector of information rates fixed before the trial. If
  unspecified, defaults to \\(1:(K+1)) / (K+1)\\.

- efficacyStopping:

  Indicators of whether efficacy stopping is allowed at each stage.
  Defaults to `TRUE` if left unspecified.

- futilityStopping:

  Indicators of whether futility stopping is allowed at each stage.
  Defaults to `TRUE` if left unspecified.

- criticalValues:

  The upper boundaries on the max-Z statistic scale for Phase 2 and the
  Z statistics for the selected arm in Phase 3. If missing, boundaries
  will be computed based on the specified alpha spending function.

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

  A numeric vector of length \\K\\ specifying futility boundaries on the
  max-Z scale at the end of Phase 2 and on the Z scale for the \\K - 1\\
  analyses in Phase 3. The final analysis uses the efficacy boundary as
  the futility boundary.

- futilityCP:

  A numeric vector of length \\K\\ specifying futility boundaries on the
  conditional power scale.

- futilityTheta:

  A numeric vector of length \\K\\ specifying futility boundaries on the
  parameter scale.

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

  A numeric vector of length \\K+1\\ specifying the error spending time
  at each analysis. Values must be strictly increasing and end at 1. If
  omitted, defaults to `informationRates`.

## Value

An S3 object of class `seamless` with the following components:

- `overallResults`: A data frame containing:

  - `overallReject`: Overall probability of rejecting the null
    hypothesis.

  - `alpha`: Overall significance level.

  - `attainedAlpha`: The attained significance level, which may differ
    from `alpha` in the presence of futility stopping.

  - `M`: Number of active arms in Phase 2.

  - `r`: Randomization ratio per active arm versus control in Phase 2.

  - `corr_known`: Whether the phase-2 correlation was assumed known.

  - `K`: Number of looks in Phase 3.

  - `information`: Maximum information for any active arm versus
    control.

  - `expectedInformationH1`: Expected information under the alternative.

  - `expectedInformationH0`: Expected information under the null.

  - `informationOverall`: Maximum information for the overall study.

  - `expectedInformationH1`: Expected information under the alternative
    for the overall study.

  - `expectedInformationH0`: Expected information under the null for the
    overall study.

- `byStageResults`: A data frame containing:

  - `informationRates`: Information rates at each analysis.

  - `efficacyBounds`: Efficacy boundaries on the Z scale.

  - `futilityBounds`: Futility boundaries on the Z scale.

  - `rejectPerStage`: Probability of efficacy stopping at each stage.

  - `futilityPerStage`: Probability of futility stopping at each stage.

  - `cumulativeRejection`: Cumulative probability of efficacy stopping.

  - `cumulativeFutility`: Cumulative probability of futility stopping.

  - `cumulativeAlphaSpent`: Cumulative alpha spent.

  - `efficacyTheta`: Efficacy boundaries on the parameter scale.

  - `futilityTheta`: Futility boundaries on the parameter scale.

  - `efficacyP`: Efficacy boundaries on the p-value scale.

  - `futilityP`: Futility boundaries on the p-value scale.

  - `information`: Cumulative information at each analysis.

  - `informationOverall`: Cumulative information for the overall study
    at each analysis.

  - `efficacyStopping`: Indicator of whether efficacy stopping is
    permitted.

  - `futilityStopping`: Indicator of whether futility stopping is
    permitted.

  - `rejectPerStageH0`: Probability of efficacy stopping under the
    global null.

  - `futilityPerStageH0`: Probability of futility stopping under the
    global null.

  - `cumulativeRejectionH0`: Cumulative probability of efficacy stopping
    under the global null.

  - `cumulativeFutilityH0`: Cumulative probability of futility stopping
    under the global null.

- `byArmResults`: A data frame containing:

  - `theta`: Parameter values for the active arms.

  - `selectAsBest`: Probability an arm is selected as best at the end of
    Phase 2.

  - `powerByArm`: Probability of rejecting the null for each arm by
    trial end.

  - `condPowerByArm`: Conditional power for each arm given it was
    selected as the best at the end of Phase 2.

- `settings`: A list of input settings:

  - `typeAlphaSpending`: Type of alpha spending function.

  - `parameterAlphaSpending`: Parameter value for the chosen alpha
    spending function.

  - `userAlphaSpending`: User-specified alpha spending values.

  - `typeBetaSpending`: Type of beta spending function.

  - `parameterBetaSpending`: Parameter value for the chosen beta
    spending function.

  - `userBetaSpending`: User-specified beta spending values.

  - `spendingTime`: Error-spending times at each analysis.

## Details

If `corr_known` is `FALSE`, critical boundaries are computed assuming
independence among the Phase-2 Wald statistics (a conservative
assumption). Power calculations, however, use the correlation implied by
the randomization ratio \\r\\.

Futility boundaries may be supplied directly on the Z scale, derived
from conditional power, derived from parameter values, or computed from
a beta spending function.

## References

Ping Gao, Yingqiu Li. Adaptive two-stage seamless sequential design for
clinical trials. Journal of Biopharmaceutical Statistics, 2025, 35(4),
565-587.

## Author

Kaifeng Lu, <kaifenglu@gmail.com>

## Examples

``` r
# Example 1: obtain the maximum information given power with no futility
(design1 <- getDesign_seamless(
  beta = 0.1, theta = c(0.3, 0.5), M = 2, r = 1.0,
  K = 2, informationRates = seq(1, 3)/3,
  alpha = 0.025, typeAlphaSpending = "OF"))
#>                                                                             
#> Phase 2/3 seamless group-sequential design                                  
#> Overall power: 0.9, overall alpha (1-sided): 0.025                          
#> Number of active arms in phase 2: 2                                         
#> Randomization ratio of each active vs. control: 1                           
#> Using correlation for critical value calculation: TRUE                      
#> Number of looks in phase 3: 2                                               
#> Max information for pairwise comparion: 54.67                               
#> Expected information under H1: 42.39, expected information under H0: 54.54  
#> Max information for oveall study: 63.78                                     
#> Expected overall info under H1: 51.5, expected overall info under H0: 63.65 
#> Alpha spending: O'Brien-Fleming, beta spending: None                        
#>                                                                             
#>                               Stage 1 Stage 2 Stage 3
#> Information rate              0.333   0.667   1.000  
#> Efficacy boundary (Z)         3.777   2.670   2.180  
#> Cumulative rejection          0.0541  0.6198  0.9000 
#> Cumulative alpha spent        0.0002  0.0066  0.0250 
#> Efficacy boundary (theta)     0.885   0.442   0.295  
#> Efficacy boundary (p)         0.0001  0.0038  0.0146 
#> Information for pairwise comp 18.22   36.44   54.67  
#> Information for overall study 27.33   45.55   63.78  
#> 
#>                           Arm 1  Arm 2 
#> Treatment effect (theta)  0.300  0.500 
#> Being the best in phase 2 0.1966 0.8034
#> Power                     0.1353 0.7647
#> Conditional power         0.6883 0.9518

# Example 2: obtain power given the maximum information and a futility rule
(design2 <- getDesign_seamless(
  IMax = 110/(2*1^2), theta = c(0.3, 0.5), M = 2, r = 1.0,
  K = 2, informationRates = seq(1, 3)/3,
  alpha = 0.025, typeAlphaSpending = "OF",
  futilityBounds = c(0.0, 0.5)))
#>                                                                              
#> Phase 2/3 seamless group-sequential design                                   
#> Overall power: 0.898, overall alpha (1-sided): 0.025                         
#> Number of active arms in phase 2: 2                                          
#> Randomization ratio of each active vs. control: 1                            
#> Using correlation for critical value calculation: TRUE                       
#> Number of looks in phase 3: 2                                                
#> Max information for pairwise comparion: 55                                   
#> Expected information under H1: 42.17, expected information under H0: 37.22   
#> Max information for oveall study: 64.17                                      
#> Expected overall info under H1: 51.34, expected overall info under H0: 46.39 
#> Alpha spending: O'Brien-Fleming, beta spending: None                         
#>                                                                              
#>                               Stage 1 Stage 2 Stage 3
#> Information rate              0.333   0.667   1.000  
#> Efficacy boundary (Z)         3.777   2.670   2.180  
#> Futility boundary (Z)         0.000   0.500   2.180  
#> Cumulative rejection          0.0548  0.6229  0.8980 
#> Cumulative futility           0.0078  0.0140  0.1020 
#> Cumulative alpha spent        0.0002  0.0066  0.0250 
#> Efficacy boundary (theta)     0.882   0.441   0.294  
#> Futility boundary (theta)     0.000   0.083   0.294  
#> Efficacy boundary (p)         0.0001  0.0038  0.0146 
#> Futility boundary (p)         0.5000  0.3085  0.0146 
#> Information for pairwise comp 18.33   36.67   55.00  
#> Information for overall study 27.50   45.83   64.17  
#> Cumulative rejection under H0 0.0002  0.0066  0.0244 
#> Cumulative futility under H0  0.3333  0.6297  0.9756 
#> 
#>                           Arm 1  Arm 2 
#> Treatment effect (theta)  0.300  0.500 
#> Being the best in phase 2 0.1959 0.8041
#> Power                     0.1349 0.7631
#> Conditional power         0.6887 0.9490

# Example 3: derive futility boundaries using beta spending
(design3 <- getDesign_seamless(
  beta = 0.1, theta = c(0.3, 0.5), M = 2, r = 1.0,
  K = 2, informationRates = seq(1, 3)/3,
  alpha = 0.025, typeAlphaSpending = "OF",
  typeBetaSpending = "sfOF"))
#>                                                                              
#> Phase 2/3 seamless group-sequential design                                   
#> Overall power: 0.9, overall alpha (1-sided): 0.025                           
#> Number of active arms in phase 2: 2                                          
#> Randomization ratio of each active vs. control: 1                            
#> Using correlation for critical value calculation: TRUE                       
#> Number of looks in phase 3: 2                                                
#> Max information for pairwise comparion: 57.14                                
#> Expected information under H1: 42.83, expected information under H0: 35.57   
#> Max information for oveall study: 66.66                                      
#> Expected overall info under H1: 52.35, expected overall info under H0: 45.09 
#> Alpha spending: O'Brien-Fleming, beta spending: Lan-DeMets O'Brien-Fleming   
#>                                                                              
#>                               Stage 1 Stage 2 Stage 3
#> Information rate              0.333   0.667   1.000  
#> Efficacy boundary (Z)         3.777   2.670   2.180  
#> Futility boundary (Z)         -0.149  1.262   2.180  
#> Cumulative rejection          0.0594  0.6437  0.9000 
#> Cumulative futility           0.0044  0.0440  0.1000 
#> Cumulative alpha spent        0.0002  0.0066  0.0250 
#> Efficacy boundary (theta)     0.865   0.433   0.288  
#> Futility boundary (theta)     -0.034  0.205   0.288  
#> Efficacy boundary (p)         0.0001  0.0038  0.0146 
#> Futility boundary (p)         0.5591  0.1034  0.0146 
#> Information for pairwise comp 19.05   38.09   57.14  
#> Information for overall study 28.57   47.61   66.66  
#> Cumulative rejection under H0 0.0002  0.0066  0.0230 
#> Cumulative futility under H0  0.2763  0.8494  0.9770 
#> 
#>                           Arm 1  Arm 2 
#> Treatment effect (theta)  0.300  0.500 
#> Being the best in phase 2 0.1914 0.8086
#> Power                     0.1325 0.7675
#> Conditional power         0.6924 0.9491
```
