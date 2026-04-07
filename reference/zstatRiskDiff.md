# Miettinen-Nurminen Score Test Statistic for Two-Sample Risk difference

Obtains the Miettinen-Nurminen score test statistic for two-sample risk
difference possibly with stratification.

## Usage

``` r
zstatRiskDiff(n1, y1, n2, y2, riskDiffH0 = 0)
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

## Value

The value of the score test statistic.

## Details

The Mantel-Haenszel sample size weights are used for stratified samples.

## Author

Kaifeng Lu, <kaifenglu@gmail.com>

## Examples

``` r
zstatRiskDiff(n1 = c(10, 10), y1 = c(4, 3),
              n2 = c(20, 10), y2 = c(2, 0), riskDiffH0 = 0)
#> [1] 2.627423
```
