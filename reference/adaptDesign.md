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
  kNew = NA_integer_,
  informationRatesNew = NA_real_,
  efficacyStoppingNew = NA_integer_,
  futilityStoppingNew = NA_integer_,
  typeAlphaSpendingNew = "sfOF",
  parameterAlphaSpendingNew = NA_real_,
  futilityBoundsNew = NULL,
  futilityCPNew = NULL,
  futilityThetaNew = NULL,
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

  The assumed parameter value.

- IMax:

  The maximum information of the primary trial. Must be provided.

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
  stopping for the primary trial. If missing, boundaries will be
  computed based on the specified alpha spending function.

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

  The lower boundaries on the z-test statistic scale for futility
  stopping for the primary trial. Defaults to `rep(-6, kMax-1)` if left
  unspecified.

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

- futilityBoundsNew:

  The lower boundaries on the z-test statistic scale for futility
  stopping for the secondary trial. Defaults to `rep(-6, kNew-1)` if
  left unspecified.

- futilityCPNew:

  The conditional power-based futility bounds for the secondary trial.

- futilityThetaNew:

  The parameter value-based futility bounds for the secondary trial.

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

- varianceRatio:

  The ratio of the variance under H0 to the variance under H1.

## Value

An `adaptDesign` object with three list components:

- `primaryTrial`: A list of selected information for the primary trial,
  including `L`, `zL`, `theta`, `maxInformation`, `kMax`,
  `informationRates`, `efficacyBounds`, `futilityBounds`, `information`,
  `alpha`, `conditionalAlpha`, `conditionalPower`, `predictivePower`,
  and and `MullerSchafer`.

- `secondaryTrial`: A `design` object for the secondary trial.

- `integratedTrial`: A list of selected information for the integrated
  trial, including `L`, `zL`, `theta`, `maxInformation`, `kMax`,
  `informationRates`, `efficacyBounds`, `futilityBounds`, and
  `information`.

## References

Lu Chi, H. M. James Hung, and Sue-Jane Wang. Modification of sample size
in group sequential clinical trials. Biometrics 1999;55:853-857.

Hans-Helge Muller and Helmut Schafer. Adaptive group sequential designs
for clinical trials: Combining the advantages of adaptive and of
classical group sequential approaches. Biometrics 2001;57:886-891.

## See also

[`getDesign`](https://kaifenglu.github.io/lrstat/reference/getDesign.md)

## Author

Kaifeng Lu, <kaifenglu@gmail.com>

## Examples

``` r
# two-arm randomized clinical trial with a normally distributed endpoint
# 90% power to detect mean difference of 15 with a standard deviation of 50
# Design the Stage I Trial with 3 looks and Lan-DeMets O'Brien-Fleming type
# spending function
delta <- 15
sigma <- 50

(des1 <- getDesignMeanDiff(
  beta = 0.1, meanDiff = delta, stDev = sigma,
  kMax = 3, alpha = 0.025, typeAlphaSpending = "sfOF"
))
#>                                                                                   
#> Group-sequential design with 3 stages for two-sample mean difference              
#> Mean difference under H0: 0, mean difference under H1: 15, standard deviation: 50 
#> Overall power: 0.9003, overall alpha (1-sided): 0.025                             
#> Drift parameter: 3.262, inflation factor: 1.012                                   
#> Maximum information: 0.05, expected under H1: 0.04, expected under H0: 0.05       
#> Maximum # subjects: 473, expected under H1: 379.2, expected under H0: 472         
#> Allocation ratio: 1                                                               
#> Alpha spending: Lan-DeMets O'Brien-Fleming, beta spending: None                   
#>                                                                                   
#>                               Stage 1 Stage 2 Stage 3
#> Information rate              0.334   0.666   1.000  
#> Efficacy boundary (Z)         3.706   2.513   1.993  
#> Cumulative rejection          0.0343  0.5596  0.9003 
#> Cumulative alpha spent        0.0001  0.0060  0.0250 
#> Number of subjects            158.0   315.0   473.0  
#> Efficacy boundary (mean diff) 29.484  14.159  9.163  
#> Efficacy boundary (p)         0.0001  0.0060  0.0231 
#> Information                   0.02    0.03    0.05   

s1 <- des1$byStageResults$informationRates
b1 <- des1$byStageResults$efficacyBounds
n <- des1$overallResults$numberOfSubjects

# Monitoring the Stage I Trial
L <- 1
nL <- des1$byStageResults$numberOfSubjects[L]
deltahat <- 8
sigmahat <- 55
sedeltahat <- sigmahat * sqrt( 4 / nL)
zL <- deltahat / sedeltahat

# Making an Adaptive Change: Stage I to Stage II
# revised clinically meaningful difference downward to 10 power the study
# retain the standard deviation at the design stage
# Muller & Schafer (2001) method to design the secondary trial
# with 2 looks and Lan-DeMets Pocock type spending function
# re-estimate sample size to reach 90% conditional power
deltaNew <- 10

(des2 <- adaptDesign(
  betaNew = 0.1, L = L, zL = zL, theta = deltaNew,
  IMax = n / (4 * sigma^2), kMax = 3, informationRates = s1,
  alpha = 0.025, typeAlphaSpending = "sfOF",
  MullerSchafer = TRUE, kNew = 2, typeAlphaSpendingNew = "sfP"
))
#>                                                      
#> Primary trial:                                       
#> Group-sequential design with 3 stages                
#> Interim adaptation look: 1, z-statistic value: 0.914 
#> Conditional type I error: 0.0378                     
#> Conditional power: 0.496, predictive power: 0.3877   
#> Muller & Schafer method for secondary trial: TRUE    
#>                                                      
#>                       Stage 1 Stage 2 Stage 3
#> Information rate      0.334   0.666   1.000  
#> Efficacy boundary (Z) 3.706   2.513   1.993  
#> Information           0.02    0.03    0.05   
#>                                                                         
#> Secondary trial:                                                        
#> Group-sequential design with 2 stages                                   
#> theta: 10, maximum information: 0.1                                     
#> Overall power: 0.9, overall significance level (1-sided): 0.0378        
#> Drift parameter: 3.227, inflation factor: 1.113                         
#> Expected information under H1: 0.07, Expected information under H0: 0.1 
#>                                                                         
#>                           Stage 1 Stage 2
#> Information rate          0.500   1.000  
#> Efficacy boundary (Z)     1.988   2.018  
#> Cumulative rejection      0.6157  0.9000 
#> Cumulative alpha spent    0.0234  0.0378 
#> Efficacy boundary (theta) 8.711   6.253  
#> Efficacy boundary (p)     0.0234  0.0218 
#> Information               0.05    0.10   
#>                                                      
#> Integrated trial:                                    
#> Group-sequential design with 3 stages                
#> Interim adaptation look: 1, z-statistic value: 0.914 
#>                                                      
#>                       Stage 1 Stage 2 Stage 3
#> Information rate      0.132   0.566   1.000  
#> Efficacy boundary (Z) 3.706   2.182   2.212  
#> Information           0.02    0.07    0.12   

INew <- des2$secondaryTrial$overallResults$information
(nNew <- ceiling(INew * 4 * sigma^2))
#> [1] 1042
(nTotal <- nL + nNew)
#> [1] 1200
```
