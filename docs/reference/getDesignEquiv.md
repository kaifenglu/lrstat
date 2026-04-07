# Power and Sample Size for a Generic Group Sequential Equivalence Design

Obtains the maximum information and stopping boundaries for a generic
group sequential equivalence design assuming a constant treatment
effect, or obtains the power given the maximum information and stopping
boundaries.

## Usage

``` r
getDesignEquiv(
  beta = NA_real_,
  IMax = NA_real_,
  thetaLower = NA_real_,
  thetaUpper = NA_real_,
  theta = 0,
  kMax = 1L,
  informationRates = NA_real_,
  criticalValues = NA_real_,
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

- IMax:

  The maximum information. Either `beta` or `IMax` should be provided
  while the other one should be missing.

- thetaLower:

  The parameter value at the lower equivalence limit.

- thetaUpper:

  The parameter value at the upper equivalence limit.

- theta:

  The parameter value under the alternative hypothesis.

- kMax:

  The maximum number of stages.

- informationRates:

  The information rates. Fixed prior to the trial. Defaults to
  `(1:kMax) / kMax` if left unspecified.

- criticalValues:

  Upper boundaries on the z-test statistic scale for stopping for
  efficacy.

- alpha:

  The significance level for each of the two one-sided tests, e.g.,
  0.05.

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
  for `"WT"`, \\rho\\ for `"sfKD"`, and \\\gamma\\ for `"sfHSD"`.

- userAlphaSpending:

  The user defined alpha spending. Cumulative alpha spent up to each
  stage.

- spendingTime:

  A vector of length `kMax` for the error spending time at each
  analysis. Defaults to missing, in which case, it is the same as
  `informationRates`.

## Value

An S3 class `designEquiv` object with three components:

- `overallResults`: A data frame containing the following variables:

  - `overallReject`: The overall rejection probability.

  - `alpha`: The overall significance level.

  - `attainedAlphaH10`: The attained significance level under H10.

  - `attainedAlphaH20`: The attained significance level under H20.

  - `kMax`: The number of stages.

  - `thetaLower`: The parameter value at the lower equivalence limit.

  - `thetaUpper`: The parameter value at the upper equivalence limit.

  - `theta`: The parameter value under the alternative hypothesis.

  - `information`: The maximum information.

  - `expectedInformationH1`: The expected information under H1.

  - `expectedInformationH10`: The expected information under H10.

  - `expectedInformationH20`: The expected information under H20.

- `byStageResults`: A data frame containing the following variables:

  - `informationRates`: The information rates.

  - `efficacyBounds`: The efficacy boundaries on the Z-scale for each of
    the two one-sided tests.

  - `rejectPerStage`: The probability for efficacy stopping.

  - `cumulativeRejection`: The cumulative probability for efficacy
    stopping.

  - `cumulativeAlphaSpent`: The cumulative alpha for each of the two
    one-sided tests.

  - `cumulativeAttainedAlphaH10`: The cumulative probability for
    efficacy stopping under H10.

  - `cumulativeAttainedAlphaH20`: The cumulative probability for
    efficacy stopping under H20.

  - `efficacyThetaLower`: The efficacy boundaries on the parameter scale
    for the one-sided null hypothesis at the lower equivalence limit.

  - `efficacyThetaUpper`: The efficacy boundaries on the parameter scale
    for the one-sided null hypothesis at the upper equivalence limit.

  - `efficacyP`: The efficacy bounds on the p-value scale for each of
    the two one-sided tests.

  - `information`: The cumulative information.

- `settings`: A list containing the following components:

  - `typeAlphaSpending`: The type of alpha spending.

  - `parameterAlphaSpending`: The parameter value for alpha spending.

  - `userAlphaSpending`: The user defined alpha spending.

  - `spendingTime`: The error spending time at each analysis.

## Details

Consider the equivalence design with two one-sided hypotheses:
\$\$H\_{10}: \theta \leq \theta\_{10},\$\$ \$\$H\_{20}: \theta \geq
\theta\_{20}.\$\$ We reject \\H\_{10}\\ at or before look \\k\\ if
\$\$Z\_{1j} = (\hat{\theta}\_j - \theta\_{10})\sqrt{I_j} \geq b_j\$\$
for some \\j=1,\ldots,k\\, where \\\\b_j:j=1,\ldots,K\\\\ are the
critical values associated with the specified alpha-spending function,
and \\I_j\\ is the information for \\\theta\\ (inverse variance of
\\\hat{\theta}\\) at the \\j\\th look. For example, for estimating the
risk difference \\\theta = \pi_1 - \pi_2\\, \$\$I_j = \left\\\frac{\pi_1
(1-\pi_1)}{n\_{1j}} + \frac{\pi_2(1-\pi_2)}{n\_{2j}}\right\\^{-1}.\$\$
It follows that \$\$(Z\_{1j} \geq b_j) = (Z_j \geq b_j +
\theta\_{10}\sqrt{I_j}),\$\$ where \\Z_j = \hat{\theta}\_j \sqrt{I_j}\\.

Similarly, we reject \\H\_{20}\\ at or before look \\k\\ if \$\$Z\_{2j}
= (\hat{\theta}\_j - \theta\_{20})\sqrt{I_j} \leq -b_j\$\$ for some
\\j=1,\ldots,k\\. We have \$\$(Z\_{2j} \leq -b_j) = (Z_j \leq - b_j +
\theta\_{20}\sqrt{I_j}).\$\$

Let \\l_j = b_j + \theta\_{10}\sqrt{I_j}\\, and \\u_j = -b_j +
\theta\_{20}\sqrt{I_j}\\. The cumulative probability to reject \\H_0 =
H\_{10} \cup H\_{20}\\ at or before look \\k\\ under the alternative
hypothesis \\H_1\\ is given by \$\$P\_\theta\left(\cup\_{j=1}^{k}
(Z\_{1j} \geq b_j) \cap \cup\_{j=1}^{k} (Z\_{2j} \leq -b_j)\right) =
p_1 + p_2 - p\_{12},\$\$ where \$\$p_1 = P\_\theta\left(\cup\_{j=1}^{k}
(Z\_{1j} \geq b_j)\right) = P\_\theta\left(\cup\_{j=1}^{k} (Z_j \geq
l_j)\right),\$\$ \$\$p_2 = P\_\theta\left(\cup\_{j=1}^{k} (Z\_{2j} \leq
-b_j)\right) = P\_\theta\left(\cup\_{j=1}^{k} (Z_j \leq u_j)\right),\$\$
and \$\$p\_{12} = P\_\theta\left(\cup\_{j=1}^{k} (Z_j \geq l_j) \cup
(Z_j \leq u_j)\right).\$\$ Of note, both \\p_1\\ and \\p_2\\ can be
evaluated using one-sided exit probabilities for group sequential
designs. If there exists \\j\leq k\\ such that \\l_j \leq u_j\\, then
\\p\_{12} = 1\\. Otherwise, \\p\_{12}\\ can be evaluated using two-sided
exit probabilities for group sequential designs.

Since the equivalent hypothesis is tested using two one-sided tests, the
type I error is controlled. To evaluate the attained type I error of the
equivalence trial under \\H\_{10}\\ (or \\H\_{20}\\), we simply fix the
control group parameters, update the active treatment group parameters
according to the null hypothesis, and use the parameters in the power
calculation outlined above.

## Author

Kaifeng Lu, <kaifenglu@gmail.com>

## Examples

``` r
# Example 1: obtain the maximum information given power
(design1 <- getDesignEquiv(
  beta = 0.2, thetaLower = log(0.8), thetaUpper = log(1.25),
  kMax = 2, informationRates = c(0.5, 1),
  alpha = 0.05, typeAlphaSpending = "sfOF"))
#>                                                                                     
#> Group-sequential design with 2 stages for equivalence test                          
#> Lower equivalence limit: -0.223, upper equivalence limit: 0.223, parameter value: 0 
#> Maximum information: 173.23                                                         
#> Overall power: 0.8, overall alpha: 0.05, attained under H10: 0.05, under H20: 0.05  
#> Expected information under H1: 173.23, under H10: 173.23, under H20: 173.23         
#> Alpha spending: Lan-DeMets O'Brien-Fleming                                          
#>                                                                                     
#>                                        Stage 1 Stage 2
#> Information rate                       0.500   1.000  
#> Boundary for each 1-sided test (Z)     2.538   1.662  
#> Cumulative rejection                   0.0000  0.8000 
#> Cumulative alpha for each 1-sided test 0.0056  0.0500 
#> Cumulative alpha attained under H10    0.0000  0.0500 
#> Cumulative alpha attained under H20    0.0000  0.0500 
#> Boundary for lower limit (theta)       0.050   -0.097 
#> Boundary for upper limit (theta)       -0.050  0.097  
#> Boundary for each 1-sided test (p)     0.0056  0.0482 
#> Information                            86.61   173.23 


# Example 2: obtain power given the maximum information
(design2 <- getDesignEquiv(
  IMax = 72.5, thetaLower = log(0.7), thetaUpper = -log(0.7),
  kMax = 3, informationRates = c(0.5, 0.75, 1),
  alpha = 0.05, typeAlphaSpending = "sfOF"))
#>                                                                                       
#> Group-sequential design with 3 stages for equivalence test                            
#> Lower equivalence limit: -0.357, upper equivalence limit: 0.357, parameter value: 0   
#> Maximum information: 72.5                                                             
#> Overall power: 0.8239, overall alpha: 0.05, attained under H10: 0.05, under H20: 0.05 
#> Expected information under H1: 63.96, under H10: 72.08, under H20: 72.08              
#> Alpha spending: Lan-DeMets O'Brien-Fleming                                            
#>                                                                                       
#>                                        Stage 1 Stage 2 Stage 3
#> Information rate                       0.500   0.750   1.000  
#> Boundary for each 1-sided test (Z)     2.538   2.016   1.720  
#> Cumulative rejection                   0.0000  0.4710  0.8239 
#> Cumulative alpha for each 1-sided test 0.0056  0.0236  0.0500 
#> Cumulative alpha attained under H10    0.0000  0.0231  0.0500 
#> Cumulative alpha attained under H20    0.0000  0.0231  0.0500 
#> Boundary for lower limit (theta)       0.065   -0.083  -0.155 
#> Boundary for upper limit (theta)       -0.065  0.083   0.155  
#> Boundary for each 1-sided test (p)     0.0056  0.0219  0.0427 
#> Information                            36.25   54.38   72.50  
```
