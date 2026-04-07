# Power for Exact Unconditional Test of Equivalence in Risk Ratio

Obtains the power given sample size for exact unconditional test of
equivalence in risk ratio.

## Usage

``` r
powerRiskRatioExactEquiv(
  n = NA_integer_,
  riskRatioLower = NA_real_,
  riskRatioUpper = NA_real_,
  pi1 = NA_real_,
  pi2 = NA_real_,
  allocationRatioPlanned = 1,
  alpha = 0.05,
  calculateAttainedAlpha = TRUE
)
```

## Arguments

- n:

  The total sample size.

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

## Value

A data frame with the following variables:

- `alpha`: The specified significance level for each of the two
  one-sided tests.

- `attainedAlphaH10`: The attained significance level under H10 if
  requested.

- `attainedAlphaH20`: The attained significance level under H20 if
  requested.

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
powerRiskRatioExactEquiv(
  n = 200, riskRatioLower = 0.8, riskRatioUpper = 1.25,
  pi1 = 0.775, pi2 = 0.775, alpha = 0.05)
#>   alpha attainedAlphaH10 attainedAlphaH20     power   n riskRatioLower
#> 1  0.05       0.04691224       0.04691224 0.7514153 200            0.8
#>   riskRatioUpper   pi1   pi2 riskRatio allocationRatioPlanned
#> 1           1.25 0.775 0.775         1                      1
#>   zstatRiskRatioLower zstatRiskRatioUpper
#> 1            1.689398           -1.689398
```
