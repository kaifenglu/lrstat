# Miettinen-Nurminen Score Confidence Interval for Two-Sample Rate Difference

Obtains the Miettinen-Nurminen score confidence interval for two-sample
rate difference possibly with stratification.

## Usage

``` r
mnRateDiffCI(t1, y1, t2, y2, cilevel = 0.95)
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

- cilevel:

  The confidence interval level.

## Value

A list with two components:

- `data` A data frame containing the input exposure and number of events
  for each treatment group. It has the following variables:

  - `t1`: The exposure for the active treatment group.

  - `y1`: The number of events for the active treatment group.

  - `t2`: The exposure for the control group.

  - `y2`: The number of events for the control group.

- `estimates`: A data frame containing the point estimate and confidence
  interval for rate difference. It has the following variables:

  - `scale`: The scale of treatment effect.

  - `estimate`: The point estimate.

  - `lower`: The lower limit of the confidence interval.

  - `upper`: The upper limit of the confidence interval.

  - `cilevel`: The confidence interval level.

## Details

The Mantel-Haenszel weights are used for stratified samples.

## Author

Kaifeng Lu, <kaifenglu@gmail.com>

## Examples

``` r
mnRateDiffCI(t1 = c(10,10), y1 = c(4,3), t2 = c(20,10), y2 = c(2,0))
#>             scale estimate      lower     upper cilevel
#> 1 rate difference      0.3 0.06511711 0.6870998    0.95
```
