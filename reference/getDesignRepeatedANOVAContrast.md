# Power and Sample Size for One-Way Repeated Measures ANOVA Contrast

Obtains the power and sample size for a single contrast in one-way
repeated measures analysis of variance.

## Usage

``` r
getDesignRepeatedANOVAContrast(
  beta = NA_real_,
  n = NA_real_,
  ngroups = 2,
  means = NA_real_,
  stDev = 1,
  corr = 0,
  contrast = NA_real_,
  meanContrastH0 = 0,
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

  The total standard deviation.

- corr:

  The correlation among the repeated measures.

- contrast:

  The coefficients for the single contrast.

- meanContrastH0:

  The mean of the contrast under the null hypothesis.

- rounding:

  Whether to round up sample size. Defaults to 1 for sample size
  rounding.

- alpha:

  The one-sided significance level. Defaults to 0.025.

## Value

An S3 class `designRepeatedANOVAContrast` object with the following
components:

- `power`: The power to reject the null hypothesis for the treatment
  contrast.

- `alpha`: The one-sided significance level.

- `n`: The number of subjects.

- `ngroups`: The number of treatment groups.

- `means`: The treatment group means.

- `stDev`: The total standard deviation.

- `corr`: The correlation among the repeated measures.

- `contrast`: The coefficients for the single contrast.

- `meanContrastH0`: The mean of the contrast under the null hypothesis.

- `meanContrast`: The mean of the contrast under the alternative
  hypothesis.

- `effectsize`: The effect size.

- `rounding`: Whether to round up sample size.

## Author

Kaifeng Lu, <kaifenglu@gmail.com>

## Examples

``` r
(design1 <- getDesignRepeatedANOVAContrast(
  beta = 0.1, ngroups = 4, means = c(1.5, 2.5, 2, 0),
  stDev = 5, corr = 0.2, contrast = c(1, 1, 1, -3)/3,
  alpha = 0.025))
#>   alpha     power  n ngroups stDev corr meanContrastH0 meanContrast effectsize
#> 1 0.025 0.9012156 71       4     5  0.2              0            2       0.15
```
