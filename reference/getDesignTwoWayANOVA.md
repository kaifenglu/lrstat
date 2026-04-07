# Power and Sample Size for Two-Way ANOVA

Obtains the power and sample size for two-way analysis of variance.

## Usage

``` r
getDesignTwoWayANOVA(
  beta = NA_real_,
  n = NA_real_,
  nlevelsA = 2,
  nlevelsB = 2,
  means = NA_real_,
  stDev = 1,
  rounding = TRUE,
  alpha = 0.05
)
```

## Arguments

- beta:

  The type II error.

- n:

  The total sample size.

- nlevelsA:

  The number of groups for Factor A.

- nlevelsB:

  The number of levels for Factor B.

- means:

  The matrix of treatment means for Factors A and B combination.

- stDev:

  The common standard deviation.

- rounding:

  Whether to round up sample size. Defaults to 1 for sample size
  rounding.

- alpha:

  The two-sided significance level. Defaults to 0.05.

## Value

An S3 class `designTwoWayANOVA` object with the following components:

- `alpha`: The two-sided significance level.

- `nlevelsA`: The number of levels for Factor A.

- `nlevelsB`: The number of levels for Factor B.

- `means`: The matrix of treatment group means.

- `stDev`: The common standard deviation.

- `effectsizeA`: The effect size for Factor A.

- `effectsizeB`: The effect size for Factor B.

- `effectsizeAB`: The effect size for Factor A and Factor B interaction.

- `rounding`: Whether to round up sample size.

- `powerdf`: The data frame containing the power and sample size
  results. It has the following variables:

  - `n`: The sample size.

  - `powerA`: The power to reject the null hypothesis that there is no
    difference among Factor A levels.

  - `powerB`: The power to reject the null hypothesis that there is no
    difference among Factor B levels.

  - `powerAB`: The power to reject the null hypothesis that there is no
    interaction between Factor A and Factor B.

## Author

Kaifeng Lu, <kaifenglu@gmail.com>

## Examples

``` r
(design1 <- getDesignTwoWayANOVA(
  beta = 0.1, nlevelsA = 2, nlevelsB = 2,
  means = matrix(c(0.5, 4.7, 0.4, 6.9), 2, 2, byrow = TRUE),
  stDev = 2, alpha = 0.05))
#>     n    powerA    powerB   powerAB alpha nlevelsA nlevelsB stDev effectsizeA
#> 1 156 0.9028350 1.0000000 0.9460967  0.05        2        2     2  0.06890625
#> 2  12 0.1267756 0.9807046 0.1424807  0.05        2        2     2  0.06890625
#> 3 132 0.8492519 1.0000000 0.9062945  0.05        2        2     2  0.06890625
#>   effectsizeB effectsizeAB
#> 1    1.788906   0.08265625
#> 2    1.788906   0.08265625
#> 3    1.788906   0.08265625
```
