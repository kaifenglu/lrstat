# Miettinen-Nurminen Score Test Statistic for Two-Sample Risk Ratio

Obtains the Miettinen-Nurminen score test statistic for two-sample risk
ratio possibly with stratification.

## Usage

``` r
zstatRiskRatio(n1, y1, n2, y2, riskRatioH0 = 1)
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

## Value

The value of the score test statistic.

## Details

The Mantel-Haenszel sample size weights are used for stratified samples.

## Author

Kaifeng Lu, <kaifenglu@gmail.com>

## Examples

``` r
zstatRiskRatio(n1 = c(10, 10), y1 = c(4, 3),
               n2 = c(20, 10), y2 = c(2, 0), riskRatioH0 = 1)
#> [1] 2.627423
```
