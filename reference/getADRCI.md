# Repeated Confidence Interval After Adaptation

Obtains the repeated p-value, conservative point estimate, and repeated
confidence interval for an adaptive group sequential trial.

## Usage

``` r
getADRCI(
  L = NA_integer_,
  zL = NA_real_,
  IMax = NA_real_,
  kMax = NA_integer_,
  informationRates = NA_real_,
  efficacyStopping = NA_integer_,
  criticalValues = NA_real_,
  alpha = 0.025,
  typeAlphaSpending = "sfOF",
  parameterAlphaSpending = NA_real_,
  spendingTime = NA_real_,
  L2 = NA_integer_,
  zL2 = NA_real_,
  INew = NA_real_,
  MullerSchafer = 0L,
  informationRatesNew = NA_real_,
  efficacyStoppingNew = NA_integer_,
  typeAlphaSpendingNew = "sfOF",
  parameterAlphaSpendingNew = NA_real_,
  spendingTimeNew = NA_real_
)
```

## Arguments

- L:

  The interim adaptation look of the primary trial.

- zL:

  The z-test statistic at the interim adaptation look of the primary
  trial.

- IMax:

  The maximum information of the primary trial.

- kMax:

  The maximum number of stages of the primary trial.

- informationRates:

  The information rates of the primary trial.

- efficacyStopping:

  Indicators of whether efficacy stopping is allowed at each stage of
  the primary trial. Defaults to true if left unspecified.

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
  Hwang, Shi & DeCani spending function, and "none" for no early
  efficacy stopping. Defaults to "sfOF".

- parameterAlphaSpending:

  The parameter value of alpha spending for the primary trial.
  Corresponds to Delta for "WT", rho for "sfKD", and gamma for "sfHSD".

- spendingTime:

  The error spending time of the primary trial. Defaults to missing, in
  which case, it is the same as `informationRates`.

- L2:

  The look of interest in the secondary trial.

- zL2:

  The z-test statistic at the look of the secondary trial.

- INew:

  The maximum information of the secondary trial.

- MullerSchafer:

  Whether to use the Muller and Schafer (2001) method for trial
  adaptation.

- informationRatesNew:

  The spacing of looks of the secondary trial.

- efficacyStoppingNew:

  The indicators of whether efficacy stopping is allowed at each look of
  the secondary trial up to look `L2`. Defaults to true if left
  unspecified.

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

- spendingTimeNew:

  The error spending time of the secondary trial. up to look `L2`.
  Defaults to missing, in which case, it is the same as
  `informationRatesNew`.

## Value

A data frame with the following variables:

- `pvalue`: Repeated p-value for rejecting the null hypothesis.

- `thetahat`: Point estimate of the parameter.

- `cilevel`: Confidence interval level.

- `lower`: Lower bound of repeated confidence interval.

- `upper`: Upper bound of repeated confidence interval.

## References

Cyrus R. Mehta, Peter Bauer, Martin Posch and Werner Brannath. Repeated
confidence intervals for adaptive group sequential trials. Stat Med.
2007;26:5422–5433.

## See also

[`adaptDesign`](https://kaifenglu.github.io/lrstat/reference/adaptDesign.md)

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

# Muller & Schafer (2001) method to design the secondary trial:
des2 = adaptDesign(
  betaNew = 0.2, L = L, zL = zL, theta = 5,
  kMax = 3, informationRates = t,
  alpha = 0.05, typeAlphaSpending = "sfHSD",
  parameterAlphaSpending = -4,
  MullerSchafer = TRUE,
  kNew = 3, typeAlphaSpendingNew = "sfHSD",
  parameterAlphaSpendingNew = -2)

n2 = ceiling(des2$secondaryTrial$overallResults$information*4*20^2)
ns = round(n2*(1:3)/3)
(des2 = adaptDesign(
  INew = n2/(4*20^2), L = L, zL = zL, theta = 5,
  kMax = 3, informationRates = t,
  alpha = 0.05, typeAlphaSpending = "sfHSD",
  parameterAlphaSpending = -4,
  MullerSchafer = TRUE,
  kNew = 3, informationRatesNew = ns/n2,
  typeAlphaSpendingNew = "sfHSD",
  parameterAlphaSpendingNew = -2))
#>                                                      
#> Primary trial:                                       
#> Group-sequential design with 3 stages                
#> Interim adaptation look: 1, z-statistic value: 1.091 
#> Conditional type I error: 0.1033                     
#> Muller & Schafer method for secondary trial: TRUE    
#>                                                      
#>                       Stage 1 Stage 2 Stage 3
#> Information rate      0.333   0.667   1.000  
#> Efficacy boundary (Z) 2.794   2.289   1.680  
#>                                                                          
#> Secondary trial:                                                         
#> Group-sequential design with 3 stages                                    
#> theta: 5, maximum information: 0.18                                      
#> Overall power: 0.8004, overall significance level (1-sided): 0.1033      
#> Drift parameter: 2.151, inflation factor: 1.043                          
#> Expected information under H1: 0.14, Expected information under H0: 0.18 
#>                                                                          
#>                           Stage 1 Stage 2 Stage 3
#> Information rate          0.334   0.666   1.000  
#> Efficacy boundary (Z)     2.160   1.783   1.351  
#> Cumulative rejection      0.1798  0.5057  0.8004 
#> Cumulative alpha spent    0.0154  0.0450  0.1033 
#> Efficacy boundary (theta) 8.683   5.081   3.142  
#> Efficacy boundary (p)     0.0154  0.0373  0.0883 
#> Information               0.06    0.12    0.18   

# termination at the second look of the secondary trial
L2 = 2
delta2 = 6.86
sigma2 = 21.77
zL2 = delta2/sqrt(4/197*sigma2^2)

t2 = des2$secondaryTrial$byStageResults$informationRates[1:L2]

# repeated confidence interval
getADRCI(L = L, zL = zL,
         IMax = n/(4*sigma1^2), kMax = 3,
         informationRates = t,
         alpha = 0.05, typeAlphaSpending = "sfHSD",
         parameterAlphaSpending = -4,
         L2 = L2, zL2 = zL2,
         INew = n2/(4*sigma2^2),
         MullerSchafer = TRUE,
         informationRatesNew = t2,
         typeAlphaSpendingNew = "sfHSD",
         parameterAlphaSpendingNew = -2)
#>       pvalue thetahat cilevel     lower    upper
#> 1 0.02051599  4.28911     0.9 0.9048009 11.32133
```
