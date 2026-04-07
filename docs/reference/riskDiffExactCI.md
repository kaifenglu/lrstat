# Exact Unconditional Confidence Interval for Risk Difference

Obtains the exact unconditional confidence interval for risk difference
based on the standardized score statistic.

## Usage

``` r
riskDiffExactCI(
  n1 = NA_integer_,
  y1 = NA_integer_,
  n2 = NA_integer_,
  y2 = NA_integer_,
  cilevel = 0.95
)
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

A data frame containing the following variables:

- `scale`: The scale of treatment effect.

- `estimate`: The point estimate.

- `lower`: The lower limit of the confidence interval.

- `upper`: The upper limit of the confidence interval.

- `cilevel`: The confidence interval level.

## Author

Kaifeng Lu, <kaifenglu@gmail.com>

## Examples

``` r
riskDiffExactCI(n1 = 30, y1 = 2, n2 = 30, y2 = 1, cilevel = 0.95)
#>             scale   estimate     lower    upper cilevel
#> 1 risk difference 0.03333333 -0.111153 0.200071    0.95
```
