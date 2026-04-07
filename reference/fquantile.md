# The Quantiles of a Survival Distribution

Obtains the quantiles of a survival distribution.

## Usage

``` r
fquantile(S, probs, ...)
```

## Arguments

- S:

  The survival function of a univariate survival time.

- probs:

  The numeric vector of probabilities.

- ...:

  Additional arguments to be passed to S.

## Value

A vector of `length(probs)` for the quantiles.

## Author

Kaifeng Lu, <kaifenglu@gmail.com>

## Examples

``` r
fquantile(pweibull, probs = c(0.25, 0.5, 0.75),
          shape = 1.37, scale = 1/0.818, lower.tail = FALSE)
#> [1] 0.4923723 0.9355369 1.5516430
```
