# Power for Poisson One-Sample Exact Test

Obtains the power for Poisson one-sample exact test.

## Usage

``` r
powerOneRateExact(
  n = NA_integer_,
  lambdaH0 = NA_real_,
  lambda = NA_real_,
  D = 1,
  alpha = 0.025
)
```

## Arguments

- n:

  The sample size.

- lambdaH0:

  The Poisson rate under the null hypothesis.

- lambda:

  The Poisson rate under the alternative hypothesis.

- D:

  The average exposure per subject.

- alpha:

  The one-sided significance level. Defaults to 0.025.

## Value

A data frame containing the critical value of the number of events for
rejecting the null hypothesis, the attained type I error, the power for
the exact test, the sample size, the Poisson rates under the null and
alternative hypotheses, the average exposure, and the direction of the
alternative.

## Author

Kaifeng Lu, <kaifenglu@gmail.com>

## Examples

``` r
powerOneRateExact(n = 525, lambdaH0 = 0.049, lambda = 0.012,
                  D = 0.5, alpha = 0.025)
#>   alpha attainedAlpha     power   n lambdaH0 lambda   D r
#> 1 0.025    0.01173728 0.9002103 525    0.049  0.012 0.5 5
```
