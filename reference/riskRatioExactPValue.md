# P-Value for Exact Unconditional Test of Risk Ratio

Obtains the p-value for exact unconditional test of risk ratio.

## Usage

``` r
riskRatioExactPValue(
  n1 = NA_integer_,
  y1 = NA_integer_,
  n2 = NA_integer_,
  y2 = NA_integer_,
  riskRatioH0 = 1,
  directionUpper = TRUE
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

- riskRatioH0:

  The risk ratio under the null hypothesis. Defaults to 1.

- directionUpper:

  Whether larger values represent better responses.

## Value

A data frame containing the following variables:

- `riskRatioH0`: The risk ratio under the null hypothesis.

- `directionUpper`: Whether larger values represent better responses.

- `riskRatio`: The observed risk ratio.

- `zstat`: The observed value of the Z test statistic.

- `pvalue`: The one-sided p-value for the unconditional exact test.

- `pi2star`: The value of pi2 that yields the p-value.

## Author

Kaifeng Lu, <kaifenglu@gmail.com>

## Examples

``` r
riskRatioExactPValue(riskRatioH0 = 1, directionUpper = 1,
                     n1 = 68, y1 = 2, n2 = 65, y2 = 1)
#>   riskRatioH0 directionUpper riskRatio     zstat    pvalue    pi2star
#> 1           1           TRUE  1.911765 0.5445955 0.3543871 0.02046496
```
