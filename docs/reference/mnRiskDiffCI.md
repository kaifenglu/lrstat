# Miettinen-Nurminen Score Confidence Interval for Two-Sample Risk Difference

Obtains the Miettinen-Nurminen score confidence interval for two-sample
risk difference possibly with stratification.

## Usage

``` r
mnRiskDiffCI(n1, y1, n2, y2, cilevel = 0.95)
```

## Arguments

- n1:

  The sample size for the active treatment group.

- y1:

  The number of responses for the active treatment group.

- n2:

  The sample size for the control group.

- y2:

  The number of responses for the control group.

- cilevel:

  The confidence interval level.

## Value

A list with two components:

- `data` A data frame containing the input sample size and number of
  responses for each treatment group. It has the following variables:

  - `n1`: The sample size for the active treatment group.

  - `y1`: The number of responses for the active treatment group.

  - `n2`: The sample size for the control group.

  - `y2`: The number of responses for the control group.

- `estimates`: A data frame containing the point estimate and confidence
  interval for risk difference. It has the following variables:

  - `scale`: The scale of treatment effect.

  - `estimate`: The point estimate.

  - `lower`: The lower limit of the confidence interval.

  - `upper`: The upper limit of the confidence interval.

  - `cilevel`: The confidence interval level.

## Details

The Mantel-Haenszel sample size weights are used for stratified samples.

## Author

Kaifeng Lu, <kaifenglu@gmail.com>

## Examples

``` r
mnRiskDiffCI(n1 = c(10, 10), y1 = c(4, 3),
             n2 = c(20, 10), y2 = c(2, 0))
#>             scale estimate      lower     upper cilevel
#> 1 risk difference      0.3 0.08253629 0.5320603    0.95
```
