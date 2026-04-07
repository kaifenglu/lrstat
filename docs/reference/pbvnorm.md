# Distribution Function of the Standard Bivariate Normal

Computes the cumulative distribution function (CDF) of the standard
bivariate normal distribution with specified lower and upper integration
limits and correlation coefficient.

## Usage

``` r
pbvnorm(lower = NA_real_, upper = NA_real_, rho = 0)
```

## Arguments

- lower:

  A numeric vector of length 2 specifying the lower limits of
  integration.

- upper:

  A numeric vector of length 2 specifying the upper limits of
  integration.

- rho:

  A numeric value specifying the correlation coefficient of the standard
  bivariate normal distribution.

## Value

A numeric value representing the probability that a standard bivariate
normal vector falls within the specified rectangular region.

## Details

This function evaluates the probability \\P(\code{lower\[1\]} \< X \<
\code{upper\[1\]}, \code{lower\[2\]} \< Y \< \code{upper\[2\]})\\ where
\\(X, Y)\\ follows a standard bivariate normal distribution with
correlation `corr`.

## Author

Kaifeng Lu, <kaifenglu@gmail.com>

## Examples

``` r
pbvnorm(c(-1, -1), c(1, 1), 0.5)
#> [1] 0.4979718
```
