# Power and Sample Size for Fisher's Exact Test for Two Proportions

Obtains the power given sample size or obtains the sample size given
power for Fisher's exact test for two proportions.

## Usage

``` r
getDesignFisherExact(
  beta = NA_real_,
  n = NA_real_,
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
(design1 <- getDesignFisherExact(
  beta = 0.2, pi1 = 0.5, pi2 = 0.2, alpha = 0.05))
#>   alpha     power  n pi1 pi2 allocationRatioPlanned
#> 1  0.05 0.8168484 87 0.5 0.2                      1
```
