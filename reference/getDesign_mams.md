# Power and Sample Size for Multi-Arm Multi-Stage Design

Computes either the maximum information and stopping boundaries for a
generic multiple comparison sequential design, or the achieved power
when the maximum information and stopping boundaries are provided.

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

  Number of active treatment arms.

- r:

  Randomization ratio of each active arm to the common control.

- corr_known:

  Logical. If `TRUE`, the correlation between Wald statistics is derived
  from the randomization ratio `r` as \\r / (r + 1)\\. If `FALSE`, a
  conservative correlation of 0 is used.

- kMax:

  Number of sequential looks.

- informationRates:

  A numeric vector of information rates fixed before the trial. If
  unspecified, defaults to \\(1:kMax) / kMax\\.

- efficacyStopping:

  Indicators of whether efficacy stopping is allowed at each stage.
  Defaults to `TRUE` if left unspecified.

- criticalValues:

  The upper boundaries on the max z-test statistic scale for efficacy
  stopping. If missing, boundaries will be computed based on the
  specified information rates and alpha spending function.

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

  A numeric vector of length \\kMax\\ specifying the error spending time
  at each analysis. Values must be strictly increasing and ends at 1. If
  omitted, defaults to `informationRates`.

## Value

An S3 object of class `mams` with these components:

- `overallResults`: A data frame containing:

  - `overallReject`: Overall probability of rejecting the null
    hypothesis.

  - `alpha`: Overall significance level.

  - `M`: Number of active arms.

  - `r`: Randomization ratio per active arm versus control.

  - `corr_known`: Whether the phase-2 correlation was assumed known.

  - `kMax`: Number of stages.

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

- `settings`: A list of input settings:

  - `typeAlphaSpending`: Type of alpha spending function.

  - `parameterAlphaSpending`: Parameter value for the chosen alpha
    spending function.

  - `userAlphaSpending`: User-specified alpha spending values.

  - `spendingTime`: Error-spending times at each analysis.

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
#> Max information for pairwise comparion: 47.15          
#> Number of looks: 3                                     
#> Alpha spending: O'Brien-Fleming                        
#>                                                        
#>                           Stage 1 Stage 2 Stage 3
#> Information rate          0.333   0.667   1.000  
#> Efficacy boundary (Z)     3.887   2.748   2.244  
#> Cumulative rejection      0.0308  0.5470  0.9000 
#> Cumulative alpha spent    0.0001  0.0058  0.0250 
#> Efficacy boundary (theta) 0.980   0.490   0.327  
#> Efficacy boundary (p)     0.0001  0.0030  0.0124 
#> Information               15.72   31.43   47.15  
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
#> Max information for pairwise comparion: 55             
#> Number of looks: 3                                     
#> Alpha spending: O'Brien-Fleming                        
#>                                                        
#>                           Stage 1 Stage 2 Stage 3
#> Information rate          0.333   0.667   1.000  
#> Efficacy boundary (Z)     3.887   2.748   2.244  
#> Cumulative rejection      0.0433  0.6329  0.9399 
#> Cumulative alpha spent    0.0001  0.0058  0.0250 
#> Efficacy boundary (theta) 0.908   0.454   0.303  
#> Efficacy boundary (p)     0.0001  0.0030  0.0124 
#> Information               18.33   36.67   55.00  
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
