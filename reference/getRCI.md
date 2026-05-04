# Repeated Confidence Interval for Group Sequential Design

Obtains the repeated confidence interval for a group sequential trial.

## Usage

``` r
getRCI(
  L = NA_integer_,
  zL = NA_real_,
  IMax = NA_real_,
  informationRates = NA_real_,
  efficacyStopping = NA_integer_,
  criticalValues = NA_real_,
  alpha = 0.025,
  typeAlphaSpending = "sfOF",
  parameterAlphaSpending = NA_real_,
  spendingTime = NA_real_
)
```

## Arguments

- L:

  The look of interest.

- zL:

  The z-test statistic at the look.

- IMax:

  The maximum information of the trial.

- informationRates:

  The information rates up to look `L`.

- efficacyStopping:

  Indicators of whether efficacy stopping is allowed at each stage up to
  look `L`. Defaults to true if left unspecified.

- criticalValues:

  The upper boundaries on the z-test statistic scale for efficacy
  stopping up to look `L`.

- alpha:

  The significance level. Defaults to 0.025.

- typeAlphaSpending:

  The type of alpha spending for the trial. One of the following: `"OF"`
  for O'Brien-Fleming boundaries, `"P"` for Pocock boundaries, `"WT"`
  for Wang & Tsiatis boundaries, `"sfOF"` for O'Brien-Fleming type
  spending function, `"sfP"` for Pocock type spending function, `"sfKD"`
  for Kim & DeMets spending function, `"sfHSD"` for Hwang, Shi & DeCani
  spending function, and `"none"` for no early efficacy stopping.
  Defaults to `"sfOF"`.

- parameterAlphaSpending:

  The parameter value for the alpha spending. Corresponds to \\\Delta\\
  for `"WT"`, \\\rho\\ for `"sfKD"`, and \\\gamma\\ for `"sfHSD"`.

- spendingTime:

  The error spending time up to look `L`. Defaults to missing, in which
  case, it is the same as `informationRates`.

## Value

A data frame with the following components:

- `pvalue`: Repeated p-value for rejecting the null hypothesis.

- `thetahat`: Point estimate of the parameter.

- `cilevel`: Confidence interval level.

- `lower`: Lower bound of repeated confidence interval.

- `upper`: Upper bound of repeated confidence interval.

If `typeAlphaSpending` is `"OF"`, `"P"`, `"WT"`, or `"none"`, then
`informationRates`, `efficacyStopping`, and `spendingTime` must be of
full length `kMax`, and `informationRates` and `spendingTime` must end
with 1.

## References

Christopher Jennison and Bruce W. Turnbull. Interim analyses: the
repeated confidence interval approach (with discussion). J R Stat Soc
Series B. 1989;51:305-361.

## Author

Kaifeng Lu, <kaifenglu@gmail.com>

## Examples

``` r
# group sequential design with 90% power to detect delta = 6
delta <- 6
sigma <- 17
n <- 282
(des1 <- getDesign(IMax = n/(4*sigma^2), theta = delta, kMax = 3,
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

# results at the second look
L <- 2
n1 <- n*2/3
delta1 <- 7
sigma1 <- 20
zL <- delta1/sqrt(4/n1*sigma1^2)

# repeated confidence interval
getRCI(L = L, zL = zL, IMax = n/(4*sigma1^2),
       informationRates = c(1/3, 2/3), alpha = 0.05,
       typeAlphaSpending = "sfHSD", parameterAlphaSpending = -4)
#>      pvalue thetahat cilevel    lower    upper
#> 1 0.0374177        7     0.9 0.322283 13.67772
```
