# Power for Exact Unconditional Test of Risk Ratio

Obtains the power given sample size for exact unconditional test of risk
ratio.

## Usage

``` r
powerRiskRatioExact(
  n = NA_integer_,
  riskRatioH0 = 1,
  pi1 = NA_real_,
  pi2 = NA_real_,
  allocationRatioPlanned = 1,
  alpha = 0.025,
  calculateAttainedAlpha = TRUE
)
```

## Arguments

- n:

  The total sample size.

- riskRatioH0:

  The risk ratio under the null hypothesis. Defaults to 1.

- pi1:

  The assumed probability for the active treatment group.

- pi2:

  The assumed probability for the control group.

- allocationRatioPlanned:

  Allocation ratio for the active treatment versus control. Defaults to
  1 for equal randomization.

- alpha:

  The one-sided significance level. Defaults to 0.025.

- calculateAttainedAlpha:

  Whether to calculate the attained alpha.

## Value

A data frame with the following variables:

- `alpha`: The specified one-sided significance level.

- `attainedAlpha`: The attained one-sided significance level if
  requested.

- `power`: The power.

- `n`: The sample size.

- `riskRatioH0`: The risk ratio under the null hypothesis.

- `pi1`: The assumed probability for the active treatment group.

- `pi2`: The assumed probability for the control group.

- `allocationRatioPlanned`: Allocation ratio for the active treatment
  versus control.

- `zstatRiskRatioBound`: The critical value on the scale of score test
  statistic for risk ratio.

- `pi2star`: The response probability in the control group at which the
  critical value of the test statistic is attained.

## Author

Kaifeng Lu, <kaifenglu@gmail.com>

## Examples

``` r
powerRiskRatioExact(n = 68, pi1 = 0.6, pi2 = 0.25, alpha = 0.05)
#>   alpha attainedAlpha     power  n riskRatioH0 pi1  pi2 allocationRatioPlanned
#> 1  0.05    0.04876955 0.9052923 68           1 0.6 0.25                      1
#>   zstatRiskRatioBound   pi2star
#> 1            1.715316 0.8662844
```
