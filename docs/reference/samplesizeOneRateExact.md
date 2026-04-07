# Sample Size for Poisson One-Sample Exact Test

Obtains the sample size for Poisson one-sample exact test.

## Usage

``` r
samplesizeOneRateExact(
  beta = 0.2,
  lambdaH0 = NA_real_,
  lambda = NA_real_,
  D = 1,
  alpha = 0.025,
  max_n_search = 1000L,
  window = 10L
)
```

## Arguments

- beta:

  The type II error.

- lambdaH0:

  The Poisson rate under the null hypothesis.

- lambda:

  The Poisson rate under the alternative hypothesis.

- D:

  The average exposure per subject.

- alpha:

  The one-sided significance level. Defaults to 0.025.

- max_n_search:

  The maximum sample size to search up to. If no sample size up to this
  value satisfies the windowed power criterion, an error is thrown.

- window:

  The number of consecutive sample sizes that must all satisfy the power
  criterion to confirm the found sample size. This is to mitigate
  non-monotonicity of power in sample size for the exact test.

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
samplesizeOneRateExact(beta = 0.2, lambdaH0 = 0.2, lambda = 0.3,
                       D = 1, alpha = 0.05)
#>   alpha attainedAlpha     power   n lambdaH0 lambda D  r
#> 1  0.05    0.04267447 0.8078048 162      0.2    0.3 1 43
```
