# Power and Sample Size for Exact Unconditional Test for Equivalence in Risk Difference

Obtains the power given sample size or obtains the sample size given
power for exact unconditional test of equivalence in risk difference.

## Usage

``` r
getDesignRiskDiffExactEquiv(
  beta = NA_real_,
  n = NA_real_,
  riskDiffLower = NA_real_,
  riskDiffUpper = NA_real_,
  pi1 = NA_real_,
  pi2 = NA_real_,
  allocationRatioPlanned = 1,
  alpha = 0.05
)
```

## Arguments

- beta:

  The type II error.

- n:

  The total sample size.

- riskDiffLower:

  The lower equivalence limit of risk difference.

- riskDiffUpper:

  The upper equivalence limit of risk difference.

- pi1:

  The assumed probability for the active treatment group.

- pi2:

  The assumed probability for the control group.

- allocationRatioPlanned:

  Allocation ratio for the active treatment versus control. Defaults to
  1 for equal randomization.

- alpha:

  The significance level for each of the two one-sided tests. Defaults
  to 0.05.

## Value

A data frame with the following variables:

- `alpha`: The specified significance level for each of the two
  one-sided tests.

- `attainedAlpha`: The attained significance level.

- `power`: The power.

- `n`: The sample size.

- `riskDiffLower`: The lower equivalence limit of risk difference.

- `riskDiffUpper`: The upper equivalence limit of risk difference.

- `pi1`: The assumed probability for the active treatment group.

- `pi2`: The assumed probability for the control group.

- `riskDiff`: The risk difference.

- `allocationRatioPlanned`: Allocation ratio for the active treatment
  versus control.

- `zstatRiskDiffLower`: The efficacy boundaries on the z-test statistic
  scale for the one-sided null hypothesis on the lower equivalence
  limit.

- `zstatRiskDiffUpper`: The efficacy boundaries on the z-test statistic
  scale for the one-sided null hypothesis on the upper equivalence
  limit.

## Author

Kaifeng Lu, <kaifenglu@gmail.com>

## Examples

``` r
getDesignRiskDiffExactEquiv(
  n = 200, riskDiffLower = -0.2, riskDiffUpper = 0.2,
  pi1 = 0.775, pi2 = 0.775, alpha = 0.05)
#>   alpha attainedAlphaH10 attainedAlphaH20     power   n riskDiffLower
#> 1  0.05        0.0495724        0.0495724 0.9146859 200          -0.2
#>   riskDiffUpper   pi1   pi2 riskDiff allocationRatioPlanned zstatRiskDiffLower
#> 1           0.2 0.775 0.775        0                      1           1.669677
#>   zstatRiskDiffUpper
#> 1          -1.669677
```
