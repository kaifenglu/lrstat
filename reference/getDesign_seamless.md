# Power and Sample Size for a Phase 2/3 Seamless Design

Computes either the maximum information and stopping boundaries for a
phase 2/3 seamless design, or the achieved power when the maximum
information and stopping boundaries are provided.

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
  criticalValues = NA_real_,
  alpha = 0.025,
  typeAlphaSpending = "sfOF",
  parameterAlphaSpending = NA_real_,
  userAlphaSpending = NA_real_,
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

- criticalValues:

  The upper boundaries on the max z-test statistic scale for Phase 2 and
  the z-test statistics for the selected arm in Phase 3 for the primary
  trial. If missing, boundaries will be computed based on the specified
  alpha spending function.

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

- spendingTime:

  A numeric vector of length \\K+1\\ specifying the error spending time
  at each analysis. Values must be strictly increasing and ends at 1. If
  omitted, defaults to `informationRates`.

## Value

An S3 object of class `seamless` with the following components:

- `overallResults`: A data frame containing:

  - `overallReject`: Overall probability of rejecting the null
    hypothesis.

  - `alpha`: Overall significance level.

  - `M`: Number of active arms in phase 2.

  - `r`: Randomization ratio per active arm versus control in phase 2.

  - `corr_known`: Whether the phase-2 correlation was assumed known.

  - `K`: Number of stages in phase 3.

  - `information`: Maximum information for any active arm versus
    control.

- `byStageResults`: A data frame containing:

  - `informationRates`: Information rates at each analysis.

  - `efficacyBounds`: Efficacy boundaries on the Z-scale.

  - `rejectPerStage`: Probability of efficacy stopping at each stage.

  - `cumulativeRejection`: Cumulative probability of efficacy stopping.

  - `cumulativeAlphaSpent`: Cumulative alpha spent.

  - `efficacyTheta`: Efficacy boundaries on the parameter scale.

  - `efficacyP`: Efficacy boundaries on the p-value scale.

  - `information`: Cumulative information for any active arm versus
    control at each analysis.

  - `efficacyStopping`: Indicator of whether efficacy stopping is
    permitted.

- `byArmResults`: A data frame containing:

  - `theta`: Parameter values for the active arms.

  - `selectAsBest`: Probability an arm is selected as best at the end of
    phase 2.

  - `powerByArm`: Probability of rejecting the null for each arm by
    trial end.

  - `condPowerByArm`: Conditional power for each arm given it was
    selected as the best at the end of phase 2.

- `settings`: A list of input settings:

  - `typeAlphaSpending`: Type of alpha spending function.

  - `parameterAlphaSpending`: Parameter value for the chosen alpha
    spending function.

  - `userAlphaSpending`: User-specified alpha spending values.

  - `spendingTime`: Error-spending times at each analysis.

## Details

If `corr_known` is `FALSE`, critical boundaries are computed assuming
independence among the phase-2 Wald statistics (a conservative
assumption). Power calculations, however, use the correlation implied by
the randomization ratio \\r\\.

## References

Ping Gao, Yingqiu Li. Adaptive two-stage seamless sequential design for
clinical trials. Journal of Biopharmaceutical Statistics, 2025, 35(4),
565-587.

## Author

Kaifeng Lu, <kaifenglu@gmail.com>

## Examples

``` r
# Example 1: obtain the maximum information given power
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
#> Max information for pairwise comparion: 54.6651641304325 
#> Number of looks in phase 3: 2                            
#> Alpha spending: O'Brien-Fleming                          
#>                                                          
#>                           Stage 1 Stage 2 Stage 3
#> Information rate          0.333   0.667   1.000  
#> Efficacy boundary (Z)     3.777   2.670   2.180  
#> Cumulative rejection      0.0541  0.6198  0.9000 
#> Cumulative alpha spent    0.0002  0.0066  0.0250 
#> Efficacy boundary (theta) 0.885   0.442   0.295  
#> Efficacy boundary (p)     0.0001  0.0038  0.0146 
#> Information               18.22   36.44   54.67  
#> 
#>                           Arm 1  Arm 2 
#> Treatment effect (theta)  0.300  0.500 
#> Being the best in phase 2 0.1966 0.8033
#> Power                     0.1353 0.7646
#> Conditional power         0.6883 0.9518

# Example 2: obtain power given the maximum information
(design2 <- getDesign_seamless(
  IMax = 110/(2*1^2), theta = c(0.3, 0.5), M = 2, r = 1.0,
  K = 2, informationRates = seq(1, 3)/3,
  alpha = 0.025, typeAlphaSpending = "OF"))
#>                                                        
#> Phase 2/3 seamless group-sequential design             
#> Overall power: 0.9016, overall alpha (1-sided): 0.025  
#> Number of active arms in phase 2: 2                    
#> Randomization ratio of each active vs. control: 1      
#> Using correlation for critical value calculation: TRUE 
#> Max information for pairwise comparion: 55             
#> Number of looks in phase 3: 2                          
#> Alpha spending: O'Brien-Fleming                        
#>                                                        
#>                           Stage 1 Stage 2 Stage 3
#> Information rate          0.333   0.667   1.000  
#> Efficacy boundary (Z)     3.777   2.670   2.180  
#> Cumulative rejection      0.0548  0.6231  0.9016 
#> Cumulative alpha spent    0.0002  0.0066  0.0250 
#> Efficacy boundary (theta) 0.882   0.441   0.294  
#> Efficacy boundary (p)     0.0001  0.0038  0.0146 
#> Information               18.33   36.67   55.00  
#> 
#>                           Arm 1  Arm 2 
#> Treatment effect (theta)  0.300  0.500 
#> Being the best in phase 2 0.1959 0.8040
#> Power                     0.1354 0.7662
#> Conditional power         0.6909 0.9529
```
