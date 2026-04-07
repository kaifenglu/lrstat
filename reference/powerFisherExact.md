# Power for Fisher's Exact Test for Two Proportions

Obtains the power given sample size for Fisher's exact test for two
proportions.

## Usage

``` r
powerFisherExact(
  n = NA_integer_,
  pi1 = NA_real_,
  pi2 = NA_real_,
  allocationRatioPlanned = 1,
  alpha = 0.05
)
```

## Arguments

- n:

  The total sample size.

- pi1:

  The assumed probability for the active treatment group.

- pi2:

  The assumed probability for the control group.

- allocationRatioPlanned:

  Allocation ratio for the active treatment versus control. Defaults to
  1 for equal randomization.

- alpha:

  The two-sided significance level. Defaults to 0.05.

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
(design1 <- powerFisherExact(
  n = 136, pi1 = 0.25, pi2 = 0.05, alpha = 0.05))
#>   alpha     power   n  pi1  pi2 allocationRatioPlanned
#> 1  0.05 0.8979826 136 0.25 0.05                      1
```
