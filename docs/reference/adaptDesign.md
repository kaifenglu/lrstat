# Adaptive Design at an Interim Look

Calculates the conditional power for specified incremental information,
given the interim results, parameter value, data-dependent changes in
the error spending function, and the number and spacing of interim
looks. Conversely, calculates the incremental information required to
attain a specified conditional power, given the interim results,
parameter value, data-dependent changes in the error spending function,
and the number and spacing of interim looks.

## Usage

``` r
adaptDesign(
  betaNew = NA_real_,
  INew = NA_real_,
  L = NA_integer_,
  zL = NA_real_,
  theta = NA_real_,
  IMax = NA_real_,
  kMax = NA_integer_,
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
  spendingTime = NA_real_,
  MullerSchafer = FALSE,
  kNew = NA_integer_,
  informationRatesNew = NA_real_,
  efficacyStoppingNew = NA_integer_,
  futilityStoppingNew = NA_integer_,
  typeAlphaSpendingNew = "sfOF",
  parameterAlphaSpendingNew = NA_real_,
  typeBetaSpendingNew = "none",
  parameterBetaSpendingNew = NA_real_,
  userBetaSpendingNew = NA_real_,
  spendingTimeNew = NA_real_,
  varianceRatio = 1
)
```

## Arguments

- betaNew:

  The type II error for the secondary trial.

- INew:

  The maximum information of the secondary trial. Either `betaNew` or
  `INew` should be provided, while the other must be missing.

- L:

  The interim adaptation look of the primary trial.

- zL:

  The z-test statistic at the interim adaptation look of the primary
  trial.

- theta:

  The parameter value.

- IMax:

  The maximum information of the primary trial. Must be provided if
  `futilityBounds` is missing and `typeBetaSpending` is not equal to
  "none", or if conditional power calculation is desired.

- kMax:

  The maximum number of stages of the primary trial.

- informationRates:

  The information rates of the primary trial.

- efficacyStopping:

  Indicators of whether efficacy stopping is allowed at each stage of
  the primary trial. Defaults to `TRUE` if left unspecified.

- futilityStopping:

  Indicators of whether futility stopping is allowed at each stage of
  the primary trial. Defaults to `TRUE` if left unspecified.

- criticalValues:

  The upper boundaries on the z-test statistic scale for efficacy
  stopping for the primary trial.

- alpha:

  The significance level of the primary trial. Defaults to 0.025.

- typeAlphaSpending:

  The type of alpha spending for the primary trial. One of the
  following: "OF" for O'Brien-Fleming boundaries, "P" for Pocock
  boundaries, "WT" for Wang & Tsiatis boundaries, "sfOF" for
  O'Brien-Fleming type spending function, "sfP" for Pocock type spending
  function, "sfKD" for Kim & DeMets spending function, "sfHSD" for
  Hwang, Shi & DeCani spending function, "user" for user defined
  spending, and "none" for no early efficacy stopping. Defaults to
  "sfOF".

- parameterAlphaSpending:

  The parameter value of alpha spending for the primary trial.
  Corresponds to Delta for "WT", rho for "sfKD", and gamma for "sfHSD".

- userAlphaSpending:

  The user-defined alpha spending for the primary trial. Represents the
  cumulative alpha spent up to each stage.

- futilityBounds:

  The lower boundaries on the z-test statistic scale for futility
  stopping for the primary trial. Defaults to `rep(-6, kMax-1)` if left
  unspecified.

- typeBetaSpending:

  The type of beta spending for the primary trial. One of the following:
  "sfOF" for O'Brien-Fleming type spending function, "sfP" for Pocock
  type spending function, "sfKD" for Kim & DeMets spending function,
  "sfHSD" for Hwang, Shi & DeCani spending function, and "none" for no
  early futility stopping. Defaults to "none".

- parameterBetaSpending:

  The parameter value of beta spending for the primary trial.
  Corresponds to rho for "sfKD", and gamma for "sfHSD".

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

- futilityStoppingNew:

  The indicators of whether futility stopping is allowed at each look of
  the secondary trial. Defaults to `TRUE` if left unspecified.

- typeAlphaSpendingNew:

  The type of alpha spending for the secondary trial. One of the
  following: "OF" for O'Brien-Fleming boundaries, "P" for Pocock
  boundaries, "WT" for Wang & Tsiatis boundaries, "sfOF" for
  O'Brien-Fleming type spending function, "sfP" for Pocock type spending
  function, "sfKD" for Kim & DeMets spending function, "sfHSD" for
  Hwang, Shi & DeCani spending function, and "none" for no early
  efficacy stopping. Defaults to "sfOF".

- parameterAlphaSpendingNew:

  The parameter value of alpha spending for the secondary trial.
  Corresponds to Delta for "WT", rho for "sfKD", and gamma for "sfHSD".

- typeBetaSpendingNew:

  The type of beta spending for the secondary trial. One of the
  following: "sfOF" for O'Brien-Fleming type spending function, "sfP"
  for Pocock type spending function, "sfKD" for Kim & DeMets spending
  function, "sfHSD" for Hwang, Shi & DeCani spending function, "user"
  for user defined spending, and "none" for no early futility stopping.
  Defaults to "none".

- parameterBetaSpendingNew:

  The parameter value of beta spending for the secondary trial.
  Corresponds to rho for "sfKD", and gamma for "sfHSD".

- userBetaSpendingNew:

  The user-defined cumulative beta spending. Represents the cumulative
  beta spent up to each stage of the secondary trial.

- spendingTimeNew:

  The error spending time of the secondary trial. Defaults to missing,
  in which case it is assumed to be the same as `informationRatesNew`.

- varianceRatio:

  The ratio of the variance under H0 to the variance under H1.

## Value

An `adaptDesign` object with two list components:

- `primaryTrial`: A list of selected information for the primary trial,
  including `L`, `zL`, `theta`, `kMax`, `informationRates`,
  `efficacyBounds`, `futilityBounds`, and `MullerSchafer`.

- `secondaryTrial`: A `design` object for the secondary trial.

## References

Lu Chi, H. M. James Hung, and Sue-Jane Wang. Modification of sample size
in group sequential clinical trials. Biometrics 1999;55:853-857.

Hans-Helge Muller and Helmut Schafer. Adaptive group sequential designs
for clinical trials: Combining the advantages of adaptive and of
classical group sequential approaches. Biometrics 2001;57:886-891.

## See also

[`getDesign`](https://github.com/kaifenglu/lrstat/reference/getDesign.md)

## Author

Kaifeng Lu, <kaifenglu@gmail.com>

## Examples

``` r
# original group sequential design with 90% power to detect delta = 6
delta = 6
sigma = 17
n = 282
(des1 = getDesign(IMax = n/(4*sigma^2), theta = delta, kMax = 3,
                  alpha = 0.05, typeAlphaSpending = "sfHSD",
                  parameterAlphaSpending = -4))
#>                                                                          
#> Group-sequential design with 3 stages                                    
#> theta: 6, maximum information: 0.24                                      
#> Overall power: 0.9029, overall alpha (1-sided): 0.05                     
#> Drift parameter: 2.963, inflation factor: 1.014                          
#> Expected information under H1: 0.19, expected information under H0: 0.24 
#> Alpha spending: HSD(gamma = -4), beta spending: None                     
#>                                                                          
#>                           Stage 1 Stage 2 Stage 3
#> Information rate          0.333   0.667   1.000  
#> Efficacy boundary (Z)     2.794   2.289   1.680  
#> Cumulative rejection      0.1395  0.5588  0.9029 
#> Cumulative alpha spent    0.0026  0.0125  0.0500 
#> Efficacy boundary (theta) 9.797   5.676   3.401  
#> Efficacy boundary (p)     0.0026  0.0110  0.0465 
#> Information               0.08    0.16    0.24   

# interim look results
L = 1
n1 = n/3
delta1 = 4.5
sigma1 = 20
zL = delta1/sqrt(4/n1*sigma1^2)

t = des1$byStageResults$informationRates

# conditional power with sample size increase
(des2 = adaptDesign(
  betaNew = NA, INew = 420/(4*sigma1^2),
  L = L, zL = zL, theta = delta1,
  IMax = n/(4*sigma1^2), kMax = 3, informationRates = t,
  alpha = 0.05, typeAlphaSpending = "sfHSD",
  parameterAlphaSpending = -4))
#>                                                      
#> Primary trial:                                       
#> Group-sequential design with 3 stages                
#> Interim adaptation look: 1, z-statistic value: 1.091 
#> Conditional type I error: 0.1033                     
#> Conditional power: 0.606, predictive power: 0.5623   
#> Muller & Schafer method for secondary trial: FALSE   
#>                                                      
#>                       Stage 1 Stage 2 Stage 3
#> Information rate      0.333   0.667   1.000  
#> Efficacy boundary (Z) 2.794   2.289   1.680  
#>                                                                          
#> Secondary trial:                                                         
#> Group-sequential design with 2 stages                                    
#> theta: 4.5, maximum information: 0.26                                    
#> Overall power: 0.8486, overall significance level (1-sided): 0.1033      
#> Drift parameter: 2.306, inflation factor: 1.011                          
#> Expected information under H1: 0.22, Expected information under H0: 0.26 
#>                                                                          
#>                           Stage 1 Stage 2
#> Information rate          0.500   1.000  
#> Efficacy boundary (Z)     2.146   1.286  
#> Cumulative rejection      0.3029  0.8486 
#> Cumulative alpha spent    0.0159  0.1033 
#> Efficacy boundary (theta) 5.925   2.510  
#> Efficacy boundary (p)     0.0159  0.0992 
#> Information               0.13    0.26   

# Muller & Schafer (2001) method to design the secondary trial:
# 3-look gamma(-2) spending with 84% power at delta = 4.5 and sigma = 20
(des2 = adaptDesign(
  betaNew = 0.16, INew = NA,
  L = L, zL = zL, theta = delta1,
  IMax = n/(4*sigma1^2), kMax = 3, informationRates = t,
  alpha = 0.05, typeAlphaSpending = "sfHSD",
  parameterAlphaSpending = -4,
  MullerSchafer = TRUE,
  kNew = 3, typeAlphaSpendingNew = "sfHSD",
  parameterAlphaSpendingNew = -2))
#>                                                      
#> Primary trial:                                       
#> Group-sequential design with 3 stages                
#> Interim adaptation look: 1, z-statistic value: 1.091 
#> Conditional type I error: 0.1033                     
#> Conditional power: 0.606, predictive power: 0.5623   
#> Muller & Schafer method for secondary trial: TRUE    
#>                                                      
#>                       Stage 1 Stage 2 Stage 3
#> Information rate      0.333   0.667   1.000  
#> Efficacy boundary (Z) 2.794   2.289   1.680  
#>                                                                         
#> Secondary trial:                                                        
#> Group-sequential design with 3 stages                                   
#> theta: 4.5, maximum information: 0.26                                   
#> Overall power: 0.84, overall significance level (1-sided): 0.1033       
#> Drift parameter: 2.303, inflation factor: 1.041                         
#> Expected information under H1: 0.2, Expected information under H0: 0.26 
#>                                                                         
#>                           Stage 1 Stage 2 Stage 3
#> Information rate          0.333   0.667   1.000  
#> Efficacy boundary (Z)     2.162   1.781   1.351  
#> Cumulative rejection      0.2028  0.5559  0.8400 
#> Cumulative alpha spent    0.0153  0.0452  0.1033 
#> Efficacy boundary (theta) 7.314   4.261   2.640  
#> Efficacy boundary (p)     0.0153  0.0375  0.0883 
#> Information               0.09    0.17    0.26   

# incremental sample size for sigma = 20
(nNew = 4*sigma1^2*des2$secondaryTrial$overallResults$information)
#> [1] 419.2447
```
