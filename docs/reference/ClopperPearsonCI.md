# Clopper-Pearson Confidence Interval for One-Sample Proportion

Obtains the Clopper-Pearson exact confidence interval for a one-sample
proportion.

## Usage

``` r
ClopperPearsonCI(n, y, cilevel = 0.95)
```

## Arguments

- n:

  The sample size.

- y:

  The number of responses.

- cilevel:

  The confidence interval level.

## Value

A data frame with the following variables:

- `n`: The sample size.

- `y`: The number of responses.

- `phat`: The observed proportion of responses.

- `lower`: The lower limit of the confidence interval.

- `upper`: The upper limit of the confidence interval.

- `cilevel`: The confidence interval level.

## References

Clopper, C. J., & Pearson, E. S. (1934). The use of confidence or
fiducial limits illustrated in the case of the binomial. Biometrika,
26(4), 404-413.

## Author

Kaifeng Lu, <kaifenglu@gmail.com>

## Examples

``` r
ClopperPearsonCI(20, 3)
#>    n y phat      lower     upper cilevel
#> 1 20 3 0.15 0.03207094 0.3789268    0.95
```
