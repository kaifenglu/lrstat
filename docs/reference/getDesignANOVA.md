# Power and Sample Size for One-Way ANOVA

Obtains the power and sample size for one-way analysis of variance.

## Usage

``` r
getDesignANOVA(
  beta = NA_real_,
  n = NA_real_,
  ngroups = 2,
  means = NA_real_,
  stDev = 1,
  allocationRatioPlanned = NA_real_,
  rounding = TRUE,
  alpha = 0.05
)
```

## Arguments

- beta:

  The type II error.

- n:

  The total sample size.

- ngroups:

  The number of treatment groups.

- means:

  The treatment group means.

- stDev:

  The common standard deviation.

- allocationRatioPlanned:

  Allocation ratio for the treatment groups. It has length `ngroups - 1`
  or `ngroups`. If it is of length `ngroups - 1`, then the last
  treatment group will assume value 1 for allocation ratio.

- rounding:

  Whether to round up sample size. Defaults to 1 for sample size
  rounding.

- alpha:

  The two-sided significance level. Defaults to 0.05.

## Value

An S3 class `designANOVA` object with the following components:

- `power`: The power to reject the null hypothesis that there is no
  difference among the treatment groups.

- `alpha`: The two-sided significance level.

- `n`: The number of subjects.

- `ngroups`: The number of treatment groups.

- `means`: The treatment group means.

- `stDev`: The common standard deviation.

- `effectsize`: The effect size.

- `allocationRatioPlanned`: Allocation ratio for the treatment groups.

- `rounding`: Whether to round up sample size.

## Details

Let \\\\\mu_i: i=1,\ldots,k\\\\ denote the group means, and \\\\r_i:
i=1,\ldots,k\\\\ denote the randomization probabilities to the \\k\\
treatment groups. Let \\\sigma\\ denote the common standard deviation,
and \\n\\ denote the total sample size. Then the \\F\\-statistic \$\$F =
\frac{SSR/(k-1)}{SSE/(n-k)} \sim F\_{k-1, n-k, \lambda}\$\$ where
\$\$\lambda = n \sum\_{i=1}^k r_i (\mu_i - \bar{\mu})^2/\sigma^2\$\$ is
the noncentrality parameter, and \\\bar{\mu} = \sum\_{i=1}^k r_i
\mu_i\\.

## Author

Kaifeng Lu, <kaifenglu@gmail.com>

## Examples

``` r
(design1 <- getDesignANOVA(
  beta = 0.1, ngroups = 4, means = c(1.5, 2.5, 2, 0),
  stDev = 3.5, allocationRatioPlanned = c(2, 2, 2, 1),
  alpha = 0.05))
#>   alpha     power   n ngroups stDev effectsize
#> 1  0.05 0.9007535 279       4   3.5 0.05164515
```
