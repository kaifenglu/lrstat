# Group Sequential Design for Equivalence in Paired Mean Ratio

Obtains the power given sample size or obtains the sample size given
power for a group sequential design for equivalence in paired mean
ratio.

## Usage

``` r
getDesignPairedMeanRatioEquiv(
  beta = NA_real_,
  n = NA_real_,
  pairedRatioLower = NA_real_,
  pairedRatioUpper = NA_real_,
  pairedRatio = 1,
  CV = 1,
  normalApproximation = TRUE,
  rounding = TRUE,
  kMax = 1L,
  informationRates = NA_real_,
  alpha = 0.05,
  typeAlphaSpending = "sfOF",
  parameterAlphaSpending = NA_real_,
  userAlphaSpending = NA_real_,
  spendingTime = NA_real_
)
```

## Arguments

- beta:

  The type II error.

- n:

  The total sample size.

- pairedRatioLower:

  The lower equivalence limit of paired ratio.

- pairedRatioUpper:

  The upper equivalence limit of paired ratio.

- pairedRatio:

  The paired ratio under the alternative hypothesis.

- CV:

  The coefficient of variation for paired ratio.

- normalApproximation:

  The type of computation of the p-values. If `TRUE`, the variance is
  assumed to be known, otherwise the calculations are performed with the
  t distribution. The exact calculation using the t distribution is only
  implemented for the fixed design.

- rounding:

  Whether to round up sample size. Defaults to 1 for sample size
  rounding.

- kMax:

  The maximum number of stages.

- informationRates:

  The information rates. Fixed prior to the trial. Defaults to
  `(1:kMax) / kMax` if left unspecified.

- alpha:

  The significance level for each of the two one-sided tests. Defaults
  to 0.05.

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

  A vector of length `kMax` for the error spending time at each
  analysis. Defaults to missing, in which case, it is the same as
  `informationRates`.

## Value

An S3 class `designPairedMeanRatioEquiv` object with three components:

- `overallResults`: A data frame containing the following variables:

  - `overallReject`: The overall rejection probability.

  - `alpha`: The significance level for each of the two one-sided tests.
    Defaults to 0.05.

  - `attainedAlpha`: The attained significance level under H0.

  - `kMax`: The number of stages.

  - `information`: The maximum information.

  - `expectedInformationH1`: The expected information under H1.

  - `expectedInformationH0`: The expected information under H0.

  - `numberOfSubjects`: The maximum number of subjects.

  - `expectedNumberOfSubjectsH1`: The expected number of subjects under
    H1.

  - `expectedNumberOfSubjectsH0`: The expected number of subjects under
    H0.

  - `pairedRatioLower`: The lower equivalence limit of paired ratio.

  - `pairedRatioUpper`: The upper equivalence limit of paired ratio.

  - `pairedRatio`: The paired ratio under the alternative hypothesis.

  - `CV`: The coefficient of variation for paired ratios.

- `byStageResults`: A data frame containing the following variables:

  - `informationRates`: The information rates.

  - `efficacyBounds`: The efficacy boundaries on the Z-scale for each of
    the two one-sided tests.

  - `rejectPerStage`: The probability for efficacy stopping.

  - `cumulativeRejection`: The cumulative probability for efficacy
    stopping.

  - `cumulativeAlphaSpent`: The cumulative alpha for each of the two
    one-sided tests.

  - `cumulativeAttainedAlpha`: The cumulative alpha attained under H0.

  - `efficacyP`: The efficacy bounds on the p-value scale for each of
    the two one-sided tests.

  - `information`: The cumulative information.

  - `numberOfSubjects`: The number of subjects.

  - `efficacyPairedRatioLower`: The efficacy boundaries on the paired
    ratio scale for the one-sided null hypothesis on the lower
    equivalence limit.

  - `efficacyPairedRatioUpper`: The efficacy boundaries on the paired
    ratio scale for the one-sided null hypothesis on the upper
    equivalence limit.

- `settings`: A list containing the following input parameters:

  - `typeAlphaSpending`: The type of alpha spending.

  - `parameterAlphaSpending`: The parameter value for alpha spending.

  - `userAlphaSpending`: The user defined alpha spending.

  - `spendingTime`: The error spending time at each analysis.

  - `normalApproximation`: The type of computation of the p-values. If
    `TRUE`, the variance is assumed to be known, otherwise the
    calculations are performed with the t distribution. The exact
    calculation using the t distribution is only implemented for the
    fixed design.

  - `rounding`: Whether to round up sample size.

## Author

Kaifeng Lu, <kaifenglu@gmail.com>

## Examples

``` r
# Example 1: group sequential trial power calculation
(design1 <- getDesignPairedMeanRatioEquiv(
  beta = 0.1, n = NA, pairedRatioLower = 0.8, pairedRatioUpper = 1.25,
  pairedRatio = 1, CV = 0.35,
  kMax = 4, alpha = 0.05, typeAlphaSpending = "sfOF"))
#>                                                                                  
#> Group-sequential design with 4 stages for equivalence in paired mean ratio       
#> Lower limit for paired ratio: 0.8, upper limit for paired ratio: 1.25            
#> Paired ratio under H1: 1, coefficient of variation for paired ratio: 0.35        
#> Overall power: 0.9029, overall alpha: 0.05, attained alpha: 0.05                 
#> Maximum information: 224.99, expected under H1: 190.4, expected under H0: 223.68 
#> Maximum # subjects: 26, expected under H1: 22, expected under H0: 25.8           
#> Alpha spending: Lan-DeMets O'Brien-Fleming                                       
#>                                                                                  
#>                                         Stage 1 Stage 2 Stage 3 Stage 4
#> Information rate                        0.231   0.500   0.769   1.000  
#> Boundary for each 1-sided test (Z)      3.916   2.539   1.983   1.726  
#> Cumulative rejection                    0.0000  0.0001  0.6661  0.9029 
#> Cumulative alpha for each 1-sided test  0.0000  0.0056  0.0254  0.0500 
#> Cumulative alpha attained under H0      0.0000  0.0000  0.0254  0.0500 
#> Number of subjects                      6.0     13.0    20.0    26.0   
#> Boundary for lower limit (paired ratio) 1.378   1.016   0.930   0.898  
#> Boundary for upper limit (paired ratio) 0.726   0.984   1.075   1.114  
#> Boundary for each 1-sided test (p)      0.0000  0.0056  0.0237  0.0422 
#> Information                             51.92   112.50  173.07  224.99 

# Example 2: sample size calculation for t-test
(design2 <- getDesignPairedMeanRatioEquiv(
  beta = 0.1, n = NA, pairedRatioLower = 0.8, pairedRatioUpper = 1.25,
  pairedRatio = 1, CV = 0.35,
  normalApproximation = FALSE, alpha = 0.05))
#>                                                                           
#> Fixed design for equivalence in paired mean ratio                         
#> Lower limit for paired ratio: 0.8, upper limit for paired ratio: 1.25     
#> Paired ratio under H1: 1, coefficient of variation for paired ratio: 0.35 
#> Overall power: 0.9061, overall alpha: 0.05, attained alpha: 0.05          
#> Information: 233.65                                                       
#> Number of subjects: 27                                                    
#>                                                                           
#>                                               
#> Boundary for each 1-sided test (t)      1.706 
#> Boundary for lower limit (paired ratio) 0.894 
#> Boundary for upper limit (paired ratio) 1.118 
#> Boundary for each 1-sided test (p)      0.0500
```
