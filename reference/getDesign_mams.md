# Power and Sample Size for a Multi-Arm Multi-Stage Design

Computes either the maximum information and stopping boundaries for a
multi-arm multi-stage design, or the achieved power when the maximum
information and stopping boundaries are provided.

## Usage

``` r
getDesign_mams(
  beta = NA_real_,
  IMax = NA_real_,
  theta = NA_real_,
  M = NA_integer_,
  r = 1,
  corr_known = TRUE,
  kMax = 1L,
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

  Number of active treatment arms.

- r:

  Randomization ratio of each active arm to the common control.

- corr_known:

  Logical. If `TRUE`, the correlation between Wald statistics is derived
  from the randomization ratio \\r\\ as \\r / (r + 1)\\. If `FALSE`, a
  conservative correlation of 0 is used.

- kMax:

  Number of sequential looks.

- informationRates:

  A numeric vector of information rates fixed before the trial. If
  unspecified, defaults to \\(1:kMax) / kMax\\.

- efficacyStopping:

  Indicators of whether efficacy stopping is allowed at each stage.
  Defaults to `TRUE` if left unspecified.

- futilityStopping:

  Indicators of whether futility stopping is allowed at each stage.
  Defaults to `TRUE` if left unspecified.

- criticalValues:

  The matrix of by-level upper boundaries on the max z-test statistic
  scale for efficacy stopping. The first column is for level `M`, the
  second column is for level `M - 1`, and so on, with the last column
  for level 1. If left unspecified, the critical values will be computed
  based on the specified alpha spending function.

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

  A numeric vector of length `kMax - 1` specifying the futility
  boundaries on the max z-test statistic scale for futility stopping.

- futilityCP:

  A numeric vector of length `kMax - 1` specifying the futility
  boundaries on the conditional power scale for futility stopping.

- futilityTheta:

  A numeric vector of length `kMax - 1` specifying the futility
  boundaries on the parameter scale for futility stopping.

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

  A numeric vector of length `kMax` specifying the error spending time
  at each analysis. Values must be strictly increasing and ends at 1. If
  omitted, defaults to `informationRates`.

## Value

An S3 object of class `mams` with the following components:

- `overallResults`: A data frame containing:

  - `overallReject`: Overall probability of rejecting the global null
    hypothesis.

  - `alpha`: Overall significance level.

  - `attainedAlpha`: The attained significance level, which is different
    from the overall significance level in the presence of futility
    stopping.

  - `M`: Number of active arms.

  - `r`: Randomization ratio per active arm versus control.

  - `corr_known`: Whether the correlation among Wald statistics was
    assumed known.

  - `kMax`: Number of stages.

  - `information`: Maximum information for any active arm versus
    control.

  - `expectedInformationH1`: The expected information under H1.

  - `expectedInformationH0`: The expected information under H0.

- `byStageResults`: A data frame containing:

  - `informationRates`: Information rates at each analysis.

  - `efficacyBounds`: Efficacy boundaries on the max Z-scale.

  - `futilityBounds`: Futility boundaries on the max Z-scale.

  - `rejectPerStage`: Probability of efficacy stopping at each stage.

  - `futilityPerStage`: Probability of futility stopping at each stage.

  - `cumulativeRejection`: Cumulative probability of efficacy stopping.

  - `cumulativeFutility`: Cumulative probability of futility stopping.

  - `cumulativeAlphaSpent`: Cumulative alpha spent.

  - `efficacyTheta`: Efficacy boundaries on the parameter scale.

  - `futilityTheta`: Futility boundaries on the parameter scale.

  - `efficacyP`: Efficacy boundaries on the p-value scale.

  - `futilityP`: Futility boundaries on the p-value scale.

  - `information`: Cumulative information for any active arm versus
    control at each analysis.

  - `efficacyStopping`: Indicator of whether efficacy stopping is
    permitted at each stage.

  - `futilityStopping`: Indicator of whether futility stopping is
    permitted at each stage.

  - `rejectPerStageH0`: The probability for efficacy stopping under H0.

  - `futilityPerStageH0`: The probability for futility stopping under
    H0.

  - `cumulativeRejectionH0`: The cumulative probability for efficacy
    stopping under H0.

  - `cumulativeFutilityH0`: The cumulative probability for futility
    stopping under H0.

- `settings`: A list of input settings:

  - `typeAlphaSpending`: The type of alpha spending.

  - `parameterAlphaSpending`: The parameter value for the chosen alpha
    spending function.

  - `userAlphaSpending`: The user-specified alpha spending values.

  - `typeBetaSpending`: The type of beta spending.

  - `parameterBetaSpending`: The parameter value for the chosen beta
    spending function.

  - `userBetaSpending`: The user-specified beta spending values.

  - `spendingTime`: The error-spending time at each analysis.

- `byLevelBounds`: A data frame containing the efficacy boundaries for
  each level of testing (i.e., number of active arms remaining) and each
  stage. Columns include:

  - `level`: Number of active arms remaining (1 to \\M\\).

  - `stage`: Stage index (1 to `kMax`).

  - `efficacyBounds`: Efficacy boundaries on the max Z-scale for the
    given level and stage.

## Details

If `corr_known` is `FALSE`, critical boundaries are computed assuming
independence among the Wald statistics in each stage (a conservative
assumption). Power calculations, however, use the correlation implied by
the randomization ratio \\r\\.

## References

Ping Gao, Yingqiu Li. Adaptive multiple comparison sequential design
(AMCSD) for clinical trials. Journal of Biopharmaceutical Statistics,
2024, 34(3), 424-440.

## Author

Kaifeng Lu, <kaifenglu@gmail.com>

## Examples

``` r
# Example 1: obtain the maximum information given power
(design1 <- getDesign_mams(
  beta = 0.1, theta = c(0.3, 0.5), M = 2, r = 1.0,
  kMax = 3, informationRates = seq(1, 3)/3,
  alpha = 0.025, typeAlphaSpending = "OF"))
#>                                                                             
#> Multi-arm multi-stage design                                                
#> Overall power: 0.9, overall alpha (1-sided): 0.025                          
#> Number of active arms: 2                                                    
#> Randomization ratio of each active vs. control: 1                           
#> Using correlation for critical value calculation: TRUE                      
#> Number of looks: 3                                                          
#> Max information for pairwise comparion: 47.15                               
#> Max information for overall study: 70.72                                    
#> Expected information under H1: 38.07, expected information under H0: 47.05  
#> Expected overall info under H1: 57.1, expected overall info under H0: 70.58 
#> Alpha spending: O'Brien-Fleming, beta spending: None                        
#>                                                                             
#>                               Stage 1 Stage 2 Stage 3
#> Information rate              0.333   0.667   1.000  
#> Efficacy boundary (Z)         3.887   2.748   2.244  
#> Cumulative rejection          0.0308  0.5470  0.9000 
#> Cumulative alpha spent        0.0001  0.0058  0.0250 
#> Efficacy boundary (theta)     0.980   0.490   0.327  
#> Efficacy boundary (p)         0.0001  0.0030  0.0124 
#> Information for pairwise comp 15.72   31.43   47.15  
#> Information for overall study 23.57   47.15   70.72  
#> 
#> By level critical values
#>   Level Stage Boundary (Z)
#> 1     2     1        3.887
#> 2     2     2        2.748
#> 3     2     3        2.244
#> 4     1     1        3.471
#> 5     1     2        2.454
#> 6     1     3        2.004

# Example 2: obtain power given the maximum information
(design2 <- getDesign_mams(
  IMax = 110/(2*1^2), theta = c(0.3, 0.5), M = 2, r = 1.0,
  kMax = 3, informationRates = seq(1, 3)/3,
  alpha = 0.025, typeAlphaSpending = "OF"))
#>                                                                             
#> Multi-arm multi-stage design                                                
#> Overall power: 0.9399, overall alpha (1-sided): 0.025                       
#> Number of active arms: 2                                                    
#> Randomization ratio of each active vs. control: 1                           
#> Using correlation for critical value calculation: TRUE                      
#> Number of looks: 3                                                          
#> Max information for pairwise comparion: 55                                  
#> Max information for overall study: 82.5                                     
#> Expected information under H1: 42.6, expected information under H0: 54.89   
#> Expected overall info under H1: 63.9, expected overall info under H0: 82.34 
#> Alpha spending: O'Brien-Fleming, beta spending: None                        
#>                                                                             
#>                               Stage 1 Stage 2 Stage 3
#> Information rate              0.333   0.667   1.000  
#> Efficacy boundary (Z)         3.887   2.748   2.244  
#> Cumulative rejection          0.0433  0.6329  0.9399 
#> Cumulative alpha spent        0.0001  0.0058  0.0250 
#> Efficacy boundary (theta)     0.908   0.454   0.303  
#> Efficacy boundary (p)         0.0001  0.0030  0.0124 
#> Information for pairwise comp 18.33   36.67   55.00  
#> Information for overall study 27.50   55.00   82.50  
#> 
#> By level critical values
#>   Level Stage Boundary (Z)
#> 1     2     1        3.887
#> 2     2     2        2.748
#> 3     2     3        2.244
#> 4     1     1        3.471
#> 5     1     2        2.454
#> 6     1     3        2.004
```
