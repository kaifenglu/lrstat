# Power for Binomial One-Sample Exact Test

Obtains the power for binomial one-sample exact test.

## Usage

``` r
powerOnePropExact(
  n = NA_integer_,
  piH0 = NA_real_,
  pi = NA_real_,
  alpha = 0.025
)
```

## Arguments

- n:

  The sample size.

- piH0:

  The response probability under the null hypothesis.

- pi:

  The response probability under the alternative hypothesis.

- alpha:

  The one-sided significance level. Defaults to 0.025.

## Value

A data frame containing the critical value of the number of responses
for rejecting the null hypothesis, the attained type I error, the power
for the exact test, the sample size, the response probabilities under
the null and alternative hypotheses, and the direction of the
alternative.

## Author

Kaifeng Lu, <kaifenglu@gmail.com>

## Examples

``` r
powerOnePropExact(n = 110, piH0 = 0.15, pi = 0.25, alpha = 0.05)
#>   alpha attainedAlpha     power   n piH0   pi  r
#> 1  0.05    0.03547524 0.8096741 110 0.15 0.25 24
```
