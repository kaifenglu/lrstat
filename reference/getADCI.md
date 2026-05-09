# Confidence Interval After Adaptation

Obtains the p-value, median unbiased point estimate, and confidence
interval after the end of an adaptive trial.

## Usage

``` r
getADCI(
  L = NA_integer_,
  zL = NA_real_,
  IMax = NA_real_,
  kMax = NA_integer_,
  informationRates = NA_real_,
  efficacyStopping = NA_integer_,
  criticalValues = NA_real_,
  alpha = 0.25,
  typeAlphaSpending = "sfOF",
  parameterAlphaSpending = NA_real_,
  spendingTime = NA_real_,
  MullerSchafer = FALSE,
  Lc = NA_integer_,
  zLc = NA_real_,
  INew = NA_real_,
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
  following: `"OF"` for O'Brien-Fleming boundaries, `"P"` for Pocock
  boundaries, `"WT"` for Wang & Tsiatis boundaries, `"sfOF"` for
  O'Brien-Fleming type spending function, `"sfP"` for Pocock type
  spending function, `"sfKD"` for Kim & DeMets spending function,
  `"sfHSD"` for Hwang, Shi & DeCani spending function, and `"none"` for
  no early efficacy stopping. Defaults to `"sfOF"`.

- parameterAlphaSpending:

  The parameter value of alpha spending for the primary trial.
  Corresponds to \\\Delta\\ for `"WT"`, \\\rho\\ for `"sfKD"`, and
  \\\gamma\\ for `"sfHSD"`.

- spendingTime:

  The error spending time of the primary trial. Defaults to missing, in
  which case, it is the same as `informationRates`.

- MullerSchafer:

  Whether to use the Muller and Schafer (2001) method for trial
  adaptation.

- Lc:

  The termination look of the integrated trial.

- zLc:

  The z-test statistic at the termination look of the integrated trial.

- INew:

  The maximum information of the secondary trial.

- informationRatesNew:

  The spacing of looks of the secondary trial up to look `L2`.

- efficacyStoppingNew:

  The indicators of whether efficacy stopping is allowed at each look of
  the secondary trial up to look `L2`. Defaults to true if left
  unspecified.

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

  The error spending time of the secondary trial up to look `L2`.
  Defaults to missing, in which case, it is the same as
  `informationRatesNew`.

## Value

A data frame with the following variables:

- `pvalue`: p-value for rejecting the null hypothesis.

- `thetahat`: Median unbiased point estimate of the parameter.

- `cilevel`: Confidence interval level.

- `lower`: Lower bound of confidence interval.

- `upper`: Upper bound of confidence interval.

## Details

If `typeAlphaSpendingNew` is `"OF"`, `"P"`, `"WT"`, or `"none"`, then
`informationRatesNew`, `efficacyStoppingNew`, and `spendingTimeNew` must
be of full length `kNew`, and `informationRatesNew` and
`spendingTimeNew` must end with 1.

## References

Ping Gao, Lingyun Liu and Cyrus Mehta. Exact inference for adaptive
group sequential designs. Stat Med. 2013;32(23):3991-4005.

## See also

[`adaptDesign`](https://kaifenglu.github.io/lrstat/reference/adaptDesign.md)

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
#> Max information: 0.05                                
#> Interim adaptation look: 1, z-statistic value: 0.914 
#> theta: 10                                            
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
#> Maximum information: 0.1                                         
#> Overall power: 0.9, overall significance level (1-sided): 0.0378 
#>                                                                  
#>                        Stage 1 Stage 2
#> Information rate       0.500   1.000  
#> Efficacy boundary (Z)  1.988   2.018  
#> Cumulative rejection   0.6157  0.9000 
#> Cumulative alpha spent 0.0234  0.0378 
#> Information            0.05    0.10   
#>                                                      
#> Integrated trial:                                    
#> Group-sequential design with 3 stages                
#> Maximum information: 0.12                            
#> Interim adaptation look: 1, z-statistic value: 0.914 
#>                                                      
#>                       Stage 1 Stage 2 Stage 3
#> Information rate      0.132   0.566   1.000  
#> Efficacy boundary (Z) 3.706   2.182   2.212  
#> Information           0.02    0.07    0.12   

INew <- des2$secondaryTrial$maxInformation
(nNew <- ceiling(INew * 4 * sigma^2))
#> [1] 1042
(nTotal <- nL + nNew)
#> [1] 1200

# Monitoring the Integrated Trial
s2 <- des2$secondaryTrial$informationRates

Lc <- 2
deltahatc <- 9.5
sigmahatc <- 52.759
L2 <- Lc - L
nL2 <-  nNew * s2[L2]
nc <- nL + nL2
sedeltahatc <- sigmahatc * sqrt(4 / nc)
zLc <- deltahatc / sedeltahatc
zL2 <- (zLc * sqrt(nc) - zL * sqrt(nL)) / sqrt(nL2)

getADCI(
  L = L, zL = zL, IMax = n / (4 * sigmahatc^2), kMax = 3,
  informationRates = s1, alpha = 0.025, typeAlphaSpending = "sfOF",
  MullerSchafer = TRUE, Lc = Lc, zLc = zLc,
  INew = nNew / (4 * sigmahatc^2), informationRatesNew = s2,
  typeAlphaSpendingNew = "sfP")
#>       pvalue thetahat cilevel    lower    upper
#> 1 0.01237657 9.349566    0.95 1.216479 17.36533
```
