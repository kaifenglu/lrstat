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
  b = NULL,
  a = NULL,
  I = NULL
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

  A vector of efficacy boundaries for the max-Z statistics.

- a:

  A vector of futility boundaries for the max-Z statistics.

- I:

  A vector of information levels for any active arm versus the common
  control.

## Value

A vector `exitProb` of length `kMax` containing the probability of
rejection at each look.

## Details

The function assumes a multivariate normal distribution for the Wald
statistics and all active arms share the same information level.

## References

Ping Gao, Yingqiu Li. Adaptive multiple comparison sequential design
(AMCSD) for clinical trials. Journal of Biopharmaceutical Statistics,
2024, 34(3), 424-440.

## Author

Kaifeng Lu, <kaifenglu@gmail.com>

## Examples

``` r
# Setup: 2 active arms vs control and 3 sequential looks.

# Information levels: equal spacing over 3 looks based on a maximum of
# 95 patients per arm, SD = 1.0
I <- 95 / (2 * 1.0^2) * seq(1, 3)/3

# O'Brien-Fleming critical values
b <- c(3.886562, 2.748214, 2.243907)

# Type I error under the global null hypothesis
p0 <- exitprob_mams(M = 2, theta = c(0, 0), kMax = 3, b = b, I = I)
cumsum(p0$exitProbUpper)
#> [1] 0.0001007465 0.0058081554 0.0250000399

# Power under alternative: Treatment effects of 0.3 and 0.5
p1 <- exitprob_mams(M = 2, theta = c(0.3, 0.5), kMax = 3, b = b, I = I)
cumsum(p1$exitProbUpper)
#> [1] 0.03130483 0.55104538 0.90221799
```
