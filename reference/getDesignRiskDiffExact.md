# Power and Sample Size for Exact Unconditional Test for Risk Difference

Obtains the power given sample size or obtains the sample size given
power for exact unconditional test of risk difference.

## Usage

``` r
getDesignRiskDiffExact(
  beta = NA_real_,
  n = NA_real_,
  riskDiffH0 = 0,
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

- riskDiffH0:

  The risk difference under the null hypothesis. Defaults to 0.

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

- `riskDiffH0`: The risk difference under the null hypothesis.

- `pi1`: The assumed probability for the active treatment group.

- `pi2`: The assumed probability for the control group.

- `allocationRatioPlanned`: Allocation ratio for the active treatment
  versus control.

- `zstatRiskDiffBound`: The critical value on the scale of score test
  statistic for risk difference.

- `pi2star`: The response probability in the control group at which the
  critical value of the test statistic is attained.

## Author

Kaifeng Lu, <kaifenglu@gmail.com>

## Examples

``` r
# Superiority test

getDesignRiskDiffExact(n = 50, pi1 = 0.6, pi2 = 0.25, alpha = 0.025)
#>   alpha attainedAlpha     power  n riskDiffH0 pi1  pi2 allocationRatioPlanned
#> 1 0.025    0.02315239 0.6946112 50          0 0.6 0.25                      1
#>   zstatRiskDiffBound  pi2star
#> 1           2.039508 0.106221


# Non-inferiority test

getDesignRiskDiffExact(beta = 0.2, riskDiffH0 = -0.3,
                       pi1 = 0.9, pi2 = 0.9, alpha = 0.025)
#>   alpha attainedAlpha     power  n riskDiffH0 pi1 pi2 allocationRatioPlanned
#> 1 0.025    0.02433985 0.8672679 47       -0.3 0.9 0.9                      1
#>   zstatRiskDiffBound   pi2star
#> 1           2.057453 0.5125817

```
