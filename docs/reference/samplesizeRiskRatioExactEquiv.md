# Sample Size for Exact Unconditional Test of Equivalence in Risk Ratio

Obtains the sample size given power for exact unconditional test of
equivalence in risk ratio.

## Usage

``` r
samplesizeRiskRatioExactEquiv(
  beta = NA_real_,
  riskRatioLower = NA_real_,
  riskRatioUpper = NA_real_,
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

- riskRatioLower:

  The lower equivalence limit of risk ratio.

- riskRatioUpper:

  The upper equivalence limit of risk ratio.

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

- `riskRatioLower`: The lower equivalence limit of risk ratio.

- `riskRatioUpper`: The upper equivalence limit of risk ratio.

- `pi1`: The assumed probability for the active treatment group.

- `pi2`: The assumed probability for the control group.

- `riskRatio`: The risk ratio.

- `allocationRatioPlanned`: Allocation ratio for the active treatment
  versus control.

- `zstatRiskRatioLower`: The efficacy boundaries on the z-test statistic
  scale for the one-sided null hypothesis on the lower equivalence
  limit.

- `zstatRiskRatioUpper`: The efficacy boundaries on the z-test statistic
  scale for the one-sided null hypothesis on the upper equivalence
  limit.

## Author

Kaifeng Lu, <kaifenglu@gmail.com>

## Examples

``` r
samplesizeRiskRatioExactEquiv(
  beta = 0.2, riskRatioLower = 0.7, riskRatioUpper = 1/0.7,
  pi1 = 0.95, pi2 = 0.95, alpha = 0.05)
#>   alpha attainedAlphaH10 attainedAlphaH20   power  n riskRatioLower
#> 1  0.05       0.04622368       0.04622368 0.82716 37            0.7
#>   riskRatioUpper  pi1  pi2 riskRatio allocationRatioPlanned zstatRiskRatioLower
#> 1       1.428571 0.95 0.95         1                      1            1.753383
#>   zstatRiskRatioUpper
#> 1           -1.753383
```
