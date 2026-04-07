# P-Value for Exact Unconditional Test of Risk Difference

Obtains the p-value for exact unconditional test of risk difference.

## Usage

``` r
riskDiffExactPValue(
  n1 = NA_integer_,
  y1 = NA_integer_,
  n2 = NA_integer_,
  y2 = NA_integer_,
  riskDiffH0 = 0,
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

- riskDiffH0:

  The risk difference under the null hypothesis. Defaults to 0.

- directionUpper:

  Whether larger values represent better responses.

## Value

A data frame containing the following variables:

- `riskDiffH0`: The risk difference under the null hypothesis.

- `directionUpper`: Whether larger values represent better responses.

- `riskDiff`: The observed risk difference.

- `zstat`: The observed value of the Z test statistic.

- `pvalue`: The one-sided p-value for the unconditional exact test.

- `pi2star`: The value of pi2 that yields the p-value.

## Author

Kaifeng Lu, <kaifenglu@gmail.com>

## Examples

``` r
riskDiffExactPValue(riskDiffH0 = 0, directionUpper = 1,
                    n1 = 68, y1 = 2, n2 = 65, y2 = 1)
#>   riskDiffH0 directionUpper   riskDiff     zstat    pvalue    pi2star
#> 1          0           TRUE 0.01402715 0.5445955 0.3543871 0.02046496
```
