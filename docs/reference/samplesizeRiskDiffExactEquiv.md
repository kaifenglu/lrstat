# Sample Size for Exact Unconditional Test of Equivalence in Risk Difference

Obtains the sample size given power for exact unconditional test of
equivalence in risk difference.

## Usage

``` r
samplesizeRiskDiffExactEquiv(
  beta = NA_real_,
  riskDiffLower = NA_real_,
  riskDiffUpper = NA_real_,
  pi1 = NA_real_,
  pi2 = NA_real_,
  allocationRatioPlanned = 1,
  alpha = 0.05,
  calculateAttainedAlpha = TRUE,
  max_n_search = 1000L,
  window = 10L
)
```

## Arguments

- beta:

  The type II error.

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

- calculateAttainedAlpha:

  Whether to calculate the attained alpha.

- max_n_search:

  The maximum sample size to search up to. If no sample size up to this
  value satisfies the windowed power criterion, an error is thrown.

- window:

  The number of consecutive sample sizes that must all satisfy the power
  criterion to confirm the found sample size. This is to mitigate
  non-monotonicity of power in sample size for the exact test.

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
samplesizeRiskDiffExactEquiv(
  beta = 0.2, riskDiffLower = -0.3, riskDiffUpper = 0.3,
  pi1 = 0.9, pi2 = 0.9, alpha = 0.05)
#>   alpha attainedAlphaH10 attainedAlphaH20    power  n riskDiffLower
#> 1  0.05       0.04239754       0.04239754 0.800266 47          -0.3
#>   riskDiffUpper pi1 pi2 riskDiff allocationRatioPlanned zstatRiskDiffLower
#> 1           0.3 0.9 0.9        0                      1           1.789338
#>   zstatRiskDiffUpper
#> 1          -1.789338
```
