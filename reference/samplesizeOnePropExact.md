# Sample Size for Binomial One-Sample Exact Test

Obtains the sample size for binomial one-sample exact test.

## Usage

``` r
samplesizeOnePropExact(
  beta = 0.2,
  piH0 = NA_real_,
  pi = NA_real_,
  alpha = 0.025,
  max_n_search = 1000L,
  window = 10L
)
```

## Arguments

- beta:

  The type II error.

- piH0:

  The response probability under the null hypothesis.

- pi:

  The response probability under the alternative hypothesis.

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

A data frame containing the critical value of the number of responses
for rejecting the null hypothesis, the attained type I error, the power
for the exact test, the sample size, the response probabilities under
the null and alternative hypotheses, and the direction of the
alternative.

## Author

Kaifeng Lu, <kaifenglu@gmail.com>

## Examples

``` r
samplesizeOnePropExact(beta = 0.2, piH0 = 0.15, pi = 0.25, alpha = 0.025)
#>   alpha attainedAlpha     power   n piH0   pi  r
#> 1 0.025     0.0179192 0.8126761 136 0.15 0.25 30
```
