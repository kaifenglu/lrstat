# Power and Sample Size for One-Way ANOVA Contrast

Obtains the power and sample size for a single contrast in one-way
analysis of variance.

## Usage

``` r
getDesignANOVAContrast(
  beta = NA_real_,
  n = NA_real_,
  ngroups = 2,
  means = NA_real_,
  stDev = 1,
  contrast = NA_real_,
  meanContrastH0 = 0,
  allocationRatioPlanned = NA_real_,
  rounding = TRUE,
  alpha = 0.025
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

- contrast:

  The coefficients for the single contrast.

- meanContrastH0:

  The mean of the contrast under the null hypothesis.

- allocationRatioPlanned:

  Allocation ratio for the treatment groups. It has length `ngroups - 1`
  or `ngroups`. If it is of length `ngroups - 1`, then the last
  treatment group will assume value 1 for allocation ratio.

- rounding:

  Whether to round up sample size. Defaults to 1 for sample size
  rounding.

- alpha:

  The one-sided significance level. Defaults to 0.025.

## Value

An S3 class `designANOVAContrast` object with the following components:

- `power`: The power to reject the null hypothesis for the treatment
  contrast.

- `alpha`: The one-sided significance level.

- `n`: The number of subjects.

- `ngroups`: The number of treatment groups.

- `means`: The treatment group means.

- `stDev`: The common standard deviation.

- `contrast`: The coefficients for the single contrast.

- `meanContrastH0`: The mean of the contrast under the null hypothesis.

- `meanContrast`: The mean of the contrast under the alternative
  hypothesis.

- `effectsize`: The effect size.

- `allocationRatioPlanned`: Allocation ratio for the treatment groups.

- `rounding`: Whether to round up sample size.

## Author

Kaifeng Lu, <kaifenglu@gmail.com>

## Examples

``` r
(design1 <- getDesignANOVAContrast(
  beta = 0.1, ngroups = 4, means = c(1.5, 2.5, 2, 0),
  stDev = 3.5, contrast = c(1, 1, 1, -3),
  allocationRatioPlanned = c(2, 2, 2, 1),
  alpha = 0.025))
#>   alpha     power   n ngroups stDev meanContrastH0 meanContrast effectsize
#> 1 0.025 0.9002752 265       4   3.5              0            6 0.03998334
```
