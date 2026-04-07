# Power and Sample Size for Exact Unconditional Test for Risk Ratio

Obtains the power given sample size or obtains the sample size given
power for exact unconditional test of risk ratio.

## Usage

``` r
getDesignRiskRatioExact(
  beta = NA_real_,
  n = NA_real_,
  riskRatioH0 = 1,
  pi1 = NA_real_,
  pi2 = NA_real_,
  allocationRatioPlanned = 1,
  alpha = 0.025
)
```

## Arguments

- beta:

  The type II error.

- n:

  The total sample size.

- riskRatioH0:

  The risk ratio under the null hypothesis. Defaults to 0.

- pi1:

  The assumed probability for the active treatment group.

- pi2:

  The assumed probability for the control group.

- allocationRatioPlanned:

  Allocation ratio for the active treatment versus control. Defaults to
  1 for equal randomization.

- alpha:

  The one-sided significance level. Defaults to 0.025.

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
# Non-inferiority test

getDesignRiskRatioExact(beta = 0.2, riskRatioH0 = 0.7,
                        pi1 = 0.95, pi2 = 0.95, alpha = 0.025)
#>   alpha attainedAlpha     power  n riskRatioH0  pi1  pi2 allocationRatioPlanned
#> 1 0.025    0.02124311 0.8573103 39         0.7 0.95 0.95                      1
#>   zstatRiskRatioBound   pi2star
#> 1            2.116876 0.1023319

```
