# Exit Probabilities for Multi-Arm Multi-Stage Design

Computes the exit (rejection) probabilities for a multi-arm multi-stage
design.

## Usage

``` r
exitprob_mams(
  M = NA_integer_,
  r = 1,
  theta = NA_real_,
  corr_known = TRUE,
  kMax = NA_integer_,
  b = NA_real_,
  I = NA_real_
)
```

## Arguments

- M:

  Number of active treatment arms.

- r:

  Randomization ratio of each active arm to the common control.

- theta:

  A vector of length \\M\\ representing the true treatment effects for
  each active arm versus the common control.

- corr_known:

  Logical. If `TRUE`, the correlation between Wald statistics is derived
  from the randomization ratio \\r\\ as \\r / (r + 1)\\. If `FALSE`, a
  conservative correlation of 0 is used.

- kMax:

  Number of sequential looks.

- b:

  A vector of critical values (length \\kMax\\).

- I:

  A vector of information levels (length \\kMax\\) for any active arm
  versus the common control.

## Value

A vector `exitProb` of length \\kMax\\ containing the probability of
rejection at each look.

## Details

The function assumes a multivariate normal distribution for the Wald
statistics and all active arms share the same information level.

## Author

Kaifeng Lu, <kaifenglu@gmail.com>

## Examples

``` r
# Setup: 2 active arms vs control and 3 sequential looks.

# Information levels: equal spacing over 3 looks based on max 110 patients
# per arm, SD = 1.0
I <- c(95 / (2 * 1.0^2) * seq(1, 3)/3)

# O'Brien-Fleming critical values
b <- c(3.886563, 2.748215, 2.243908)

# Type I error under the global null hypothesis
p0 <- exitprob_mams(M = 2, theta = c(0, 0), kMax = 3, b = b, I = I)
cumsum(p0)
#> [1] 0.0001007491 0.0058081431 0.0249999769

# Power under alternative: Treatment effects of 0.3 and 0.5
p1 <- exitprob_mams(M = 2, theta = c(0.3, 0.5), kMax = 3, b = b, I = I)
cumsum(p1)
#> [1] 0.03130476 0.55104497 0.90221780
```
