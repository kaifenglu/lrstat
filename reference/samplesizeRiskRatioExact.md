# Sample Size for Exact Unconditional Test of Risk Ratio

Obtains the sample size given power for exact unconditional test of risk
ratio.

## Usage

``` r
samplesizeRiskRatioExact(
  beta = NA_real_,
  riskRatioH0 = 1,
  pi1 = NA_real_,
  pi2 = NA_real_,
  allocationRatioPlanned = 1,
  alpha = 0.025,
  calculateAttainedAlpha = TRUE,
  max_n_search = 1000L,
  window = 10L
)
```

## Arguments

- beta:

  The type II error.

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

  The one-sided significance level.

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

- `alpha`: The specified one-sided significance level.

- `attainedAlpha`: The attained one-sided significance level.

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
samplesizeRiskRatioExact(beta = 0.2, riskRatioH0 = 0.8,
                         pi1 = 0.95, pi2 = 0.95, alpha = 0.05)
#>   alpha attainedAlpha     power  n riskRatioH0  pi1  pi2 allocationRatioPlanned
#> 1  0.05    0.04790715 0.8092953 51         0.8 0.95 0.95                      1
#>   zstatRiskRatioBound   pi2star
#> 1            1.713737 0.5113777
```
