# Power and Sample Size for One-Sample Poisson Rate Exact Test

Obtains the power given sample size or obtains the sample size given
power for one-sample Poisson rate.

## Usage

``` r
getDesignOneRateExact(
  beta = NA_real_,
  n = NA_real_,
  lambdaH0 = NA_real_,
  lambda = NA_real_,
  D = 1,
  alpha = 0.025
)
```

## Arguments

- beta:

  The type II error.

- n:

  The total sample size.

- lambdaH0:

  The Poisson rate under the null hypothesis.

- lambda:

  The Poisson rate under the alternative hypothesis.

- D:

  The average exposure per subject.

- alpha:

  The one-sided significance level. Defaults to 0.025.

## Value

A data frame containing the following variables:

- `alpha`: The specified significance level.

- `attainedAlpha`: The attained type I error of the exact test.

- `power`: The actual power of the exact test.

- `n`: The sample size.

- `lambdaH0`: The Poisson rate under the null hypothesis.

- `lambda`: The Poisson rate under the alternative hypothesis.

- `D`: The average exposure per subject.

- `r`: The critical value of the number of events for rejecting the null
  hypothesis. Reject H0 if `Y >= r` for upper-tailed test, and reject H0
  if `Y <= r` for lower-tailed test.

## Author

Kaifeng Lu, <kaifenglu@gmail.com>

## Examples

``` r
# Example 1: power calculation
(design1 <- getDesignOneRateExact(
  n = 525, lambdaH0 = 0.049, lambda = 0.012,
  D = 0.5, alpha = 0.025))
#>   alpha attainedAlpha     power   n lambdaH0 lambda   D r
#> 1 0.025    0.01173728 0.9002103 525    0.049  0.012 0.5 5

# Example 2: sample size calculation
(design2 <- getDesignOneRateExact(
  beta = 0.2, lambdaH0 = 0.2, lambda = 0.3,
  D = 1, alpha = 0.05))
#>   alpha attainedAlpha     power   n lambdaH0 lambda D  r
#> 1  0.05    0.04267447 0.8078048 162      0.2    0.3 1 43
```
