# Adaptive Two-Stage Seamless Sequential Design

Calculates the conditional power for specified incremental information,
given the interim results, parameter value, data-dependent changes in
the error spending function, and the number and spacing of interim
looks. Conversely, calculates the incremental information required to
attain a specified conditional power, given the interim results,
parameter value, data-dependent changes in the error spending function,
and the number and spacing of interim looks.

## Usage

``` r
adaptDesign_seamless(
  betaNew = NA_real_,
  INew = NA_real_,
  M = NA_integer_,
  r = NA_real_,
  corr_known = TRUE,
  L = NA_integer_,
  zL = NA_real_,
  theta = NA_real_,
  IMax = NA_real_,
  K = NA_integer_,
  informationRates = NA_real_,
  efficacyStopping = NA_integer_,
  criticalValues = NA_real_,
  alpha = 0.025,
  typeAlphaSpending = "sfOF",
  parameterAlphaSpending = NA_real_,
  userAlphaSpending = NA_real_,
  spendingTime = NA_real_,
  MullerSchafer = FALSE,
  kNew = NA_integer_,
  informationRatesNew = NA_real_,
  efficacyStoppingNew = NA_integer_,
  typeAlphaSpendingNew = "sfOF",
  parameterAlphaSpendingNew = NA_real_,
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

  Number of active treatment arms in Phase 2.

- r:

  Randomization ratio of each active arm to the common control in Phase
  2.

- corr_known:

  Logical. If `TRUE`, the correlation between Wald statistics in Phase 2
  is derived from the randomization ratio `r` as \\r / (r + 1)\\. If
  `FALSE`, a conservative correlation of 0 is used.

- L:

  The interim adaptation look in Phase 3.

- zL:

  The z-test statistic at the interim adaptation look of Phase 3.

- theta:

  The treatment effect for the selected arm versus the common control.

- IMax:

  Maximum information for any active arm versus the common control for
  the original trial. Must be provided.

- K:

  Number of sequential looks in Phase 3.

- informationRates:

  A numeric vector of information rates fixed before the trial. If
  unspecified, defaults to \\(1:(K+1)) / (K+1)\\.

- efficacyStopping:

  Indicators of whether efficacy stopping is allowed at each stage of
  the primary trial. Defaults to `TRUE` if left unspecified.

- criticalValues:

  The upper boundaries on the max z-test statistic scale for Phase 2 and
  the z-test statistics for the selected arm in Phase 3 for the primary
  trial. If missing, boundaries will be computed based on the specified
  alpha spending function.

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

- spendingTime:

  The error spending time of the primary trial. Defaults to missing, in
  which case it is assumed to be the same as `informationRates`.

- MullerSchafer:

  Whether to use the Muller and Schafer (2001) method for trial
  adaptation.

- kNew:

  The number of looks of the secondary trial.

- informationRatesNew:

  The spacing of looks of the secondary trial.

- efficacyStoppingNew:

  The indicators of whether efficacy stopping is allowed at each look of
  the secondary trial. Defaults to `TRUE` if left unspecified.

- typeAlphaSpendingNew:

  The type of alpha spending for the secondary trial. One of the
  following: `"OF"` for O'Brien-Fleming boundaries, `"P"` for Pocock
  boundaries, `"WT"` for Wang & Tsiatis boundaries, `"sfOF"` for
  O'Brien-Fleming type spending function, `"sfP"` for Pocock type
  spending function, `"sfKD"` for Kim & DeMets spending function,
  `"sfHSD"` for Hwang, Shi & DeCani spending function, and `"none"` for
  no early efficacy stopping. Defaults to `"sfOF"`.

- parameterAlphaSpendingNew:

  The parameter value of alpha spending for the secondary trial.
  Corresponds to \\\Delta\\ for `"WT"`, \\\rho\\ for `"sfKD"`, and
  \\\gamma\\ for `"sfHSD"`.

- spendingTimeNew:

  The error spending time of the secondary trial. Defaults to missing,
  in which case it is assumed to be the same as `informationRatesNew`.

## Value

An `adaptDesign_seamless` object with three list components:

- `primaryTrial`: A list of selected information for the primary trial,
  including `M`, `r`, `corr_known`, `K`, `L`, `zL`, `theta`,
  `maxInformation`, `kMax`, `informationRates`, `efficacyBounds`,
  `information`, `alpha`, `conditionalAlpha`, `conditionalPower`, and
  `MullerSchafer`.

- `secondaryTrial`: A `design` object for the secondary trial.

- `integratedTrial`: A list of selected information for the integrated
  trial, including `M`, `r`, `corr_known`, `K`, `L`, `zL`, `theta`,
  `maxInformation`, `kMax`, `informationRates`, `efficacyBounds`, and
  `information`.

## References

Ping Gao, Yingqiu Li. Adaptive two-stage seamless sequential design for
clinical trials. Journal of Biopharmaceutical Statistics, 2025, 35(4),
565-587.

## See also

[`getDesign_seamless`](https://kaifenglu.github.io/lrstat/reference/getDesign_seamless.md)

## Author

Kaifeng Lu, <kaifenglu@gmail.com>

## Examples

``` r
(des1 <- adaptDesign_seamless(
  betaNew = 0.1, M = 2, r = 1, corr_known = FALSE,
  L = 1, zL = -log(0.67) * sqrt(80 / 4), theta = -log(0.691),
  IMax = 120 / 4, K = 2, informationRates = c(1/3, 2/3, 1),
  alpha = 0.025, typeAlphaSpending = "OF", kNew = 1))
#>                                                                 
#> Primary trial:                                                  
#> Phase 2/3 seamless group-sequential design                      
#> Number of active arms in phase 2: 2                             
#> Randomization ratio of each active vs. control: 1               
#> Using correlation for critical value calculation: FALSE         
#> Max information for pairwise comparion: 30                      
#> Number of looks in phase 3: 2                                   
#> Interim adaptation look in Phase 3: 1, z-statistic value: 1.791 
#> Conditional type I error: 0.0935, conditional power: 0.44       
#> Muller & Schafer method for secondary trial: FALSE              
#>                                                                 
#>                       Stage 1 Stage 2 Stage 3
#> Information rate      0.333   0.667   1.000  
#> Efficacy boundary (Z) 3.852   2.724   2.224  
#> Information           10.00   20.00   30.00  
#>                                                                  
#> Secondary trial:                                                 
#> Fixed design                                                     
#> theta: 0.37, maximum information: 49.51                          
#> Overall power: 0.9, overall significance level (1-sided): 0.0935 
#> Drift parameter: 2.601, inflation factor: 1                      
#>                                                                  
#>                                 
#> Efficacy boundary (Z)     1.319 
#> Efficacy boundary (theta) 0.187 
#> Efficacy boundary (p)     0.0935
#>                                                                 
#> Integrated trial:                                               
#> Adaptive Phase 2/3 seamless design                              
#> Total number of looks in Phase 3: 2                             
#> Interim adaptation look in Phase 3: 1, z-statistic value: 1.791 
#>                                                                 
#>                     Stage 1 Stage 2 Stage 3
#> Information rate    0.144   0.288   1.000  
#> Efficacy bounds (Z) 3.852   2.724   2.074  
#> Information         10.00   20.00   69.51  
```
