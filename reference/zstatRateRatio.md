# Miettinen-Nurminen Score Test Statistic for Two-Sample Rate Ratio

Obtains the Miettinen-Nurminen score test statistic for two-sample rate
ratio possibly with stratification.

## Usage

``` r
zstatRateRatio(t1, y1, t2, y2, rateRatioH0 = 1)
```

## Arguments

- t1:

  The exposure for the active treatment group.

- y1:

  The number of events for the active treatment group.

- t2:

  The exposure for the control group.

- y2:

  The number of events for the control group.

- rateRatioH0:

  The rate ratio under the null hypothesis. Defaults to 1.

## Value

The value of the score test statistic.

## Details

The Mantel-Haenszel weights are used for stratified samples.

## Author

Kaifeng Lu, <kaifenglu@gmail.com>

## Examples

``` r
zstatRateRatio(t1 = c(10, 10), y1 = c(4, 3),
               t2 = c(20, 10), y2 = c(2, 0), rateRatioH0 = 1)
#> [1] 2.424871
```
