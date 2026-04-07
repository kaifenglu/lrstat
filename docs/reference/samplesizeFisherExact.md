# Sample Size for Fisher's Exact Test for Two Proportions

Obtains the sample size given power for Fisher's exact test for two
proportions.

## Usage

``` r
samplesizeFisherExact(
  beta = NA_real_,
  pi1 = NA_real_,
  pi2 = NA_real_,
  allocationRatioPlanned = 1,
  alpha = 0.05,
  max_n_search = 1000L,
  window = 10L
)
```

## Arguments

- beta:

  The type II error.

- pi1:

  The assumed probability for the active treatment group.

- pi2:

  The assumed probability for the control group.

- allocationRatioPlanned:

  Allocation ratio for the active treatment versus control. Defaults to
  1 for equal randomization.

- alpha:

  The two-sided significance level. Defaults to 0.05.

- max_n_search:

  The maximum sample size to search up to. If no sample size up to this
  value satisfies the windowed power criterion, an error is thrown.

- window:

  The number of consecutive sample sizes that must all satisfy the power
  criterion to confirm the found sample size. This is to mitigate
  non-monotonicity of power in sample size for the exact test.

## Value

A data frame with the following variables:

- `alpha`: The two-sided significance level.

- `power`: The power.

- `n`: The sample size.

- `pi1`: The assumed probability for the active treatment group.

- `pi2`: The assumed probability for the control group.

- `allocationRatioPlanned`: Allocation ratio for the active treatment
  versus control.

## Author

Kaifeng Lu, <kaifenglu@gmail.com>

## Examples

``` r
(design1 <- samplesizeFisherExact(
  beta = 0.1, pi1 = 0.25, pi2 = 0.05, alpha = 0.05))
#>   alpha     power   n  pi1  pi2 allocationRatioPlanned
#> 1  0.05 0.9065643 137 0.25 0.05                      1
```
