# Exit Probabilities for Phase 2/3 Seamless Design

Computes the upper and lower exit probabilities for a phase 2/3 seamless
design. In Phase 2, multiple active arms are compared against a common
control arm. If the phase-2 max-Z statistic crosses the efficacy
boundary, the trial stops early for efficacy; if it falls below the
futility boundary, the trial stops early for futility. Otherwise, the
best-performing arm is selected to proceed to Phase 3, where it is
tested against the common control over multiple looks with upper and
optional lower stopping boundaries.

## Usage

``` r
exitprob_seamless(
  M = NA_integer_,
  r = 1,
  theta = NA_real_,
  corr_known = TRUE,
  K = NA_integer_,
  b = NULL,
  a = NULL,
  I = NULL
)
```

## Arguments

- M:

  Number of active treatment arms in Phase 2.

- r:

  Randomization ratio of each active arm to the common control in Phase
  2.

- theta:

  A vector of length \\M\\ representing the true treatment effects for
  each active arm versus the common control.

- corr_known:

  Logical. If `TRUE`, the correlation between Wald statistics in Phase 2
  is derived from the randomization ratio \\r\\ as \\r / (r + 1)\\. If
  `FALSE`, a conservative correlation of 0 is used.

- K:

  Number of sequential looks in Phase 3.

- b:

  A vector of efficacy boundaries (length \\K+1\\). The first element is
  the efficacy boundary for the phase-2 max-Z statistic; the remaining
  \\K\\ elements are efficacy boundaries for the selected arm in Phase
  3.

- a:

  An optional vector of futility boundaries (length \\K+1\\). The first
  element is the futility boundary for the phase-2 max-Z statistic; the
  remaining \\K\\ elements are futility boundaries for the selected arm
  in Phase 3. If omitted, no futility stopping is applied.

- I:

  A vector of information levels (length \\K+1\\) for any active arm
  versus the common control. The first element is for Phase 2; the
  remaining \\K\\ elements are for the looks in Phase 3.

## Value

A list containing the following components:

- `exitProbUpper`: A vector of length \\K + 1\\. The first element is
  the probability of stopping for efficacy in Phase 2; the remaining
  elements are the probabilities of stopping for efficacy at each look
  in Phase 3.

- `exitProbLower`: A vector of length \\K + 1\\. The first element is
  the probability of stopping for futility in Phase 2; the remaining
  elements are the probabilities of stopping for futility at each look
  in Phase 3.

- `exitProbByArmUpper`: A \\(K+1) \times M\\ matrix. The \\(k, m)\\-th
  entry gives the probability of stopping for efficacy at look \\k\\
  given that arm \\m\\ is selected as best.

- `exitProbByArmLower`: A \\(K+1) \times M\\ matrix. The \\(k, m)\\-th
  entry gives the probability of stopping for futility at look \\k\\
  given that arm \\m\\ is selected as best.

- `selectAsBest`: A vector of length \\M\\ containing the probability
  that each active arm is selected to move on to Phase 3.

For backward compatibility, the list also contains:

- `exitProb`: identical to `exitProbUpper`.

- `exitProbByArm`: identical to `exitProbByArmUpper`.

## Details

The function assumes a multivariate normal distribution for the Wald
statistics. The "best" arm is defined as the active arm with the largest
Z-statistic at the end of Phase 2 among designs that continue beyond the
phase-2 analysis.

**Decision Rules:**

- **Phase 2 efficacy stop**: reject if the phase-2 max-Z statistic
  satisfies \\\max_m Z_m(I_0) \ge b_0\\.

- **Phase 2 futility stop**: stop for futility if the phase-2 max-Z
  statistic satisfies \\\max_m Z_m(I_0) \le a_0\\.

- **Continue to Phase 3**: if \\a_0 \< \max_m Z_m(I_0) \< b_0\\, select
  the arm with the largest phase-2 Z-statistic and continue with that
  arm only.

- **Phase 3 efficacy stop**: at look \\k\\, reject if the selected arm's
  Z-statistic exceeds the efficacy boundary and no earlier stop has
  occurred.

- **Phase 3 futility stop**: at look \\k\\, stop for futility if the
  selected arm's Z-statistic is below the futility boundary and no
  earlier stop has occurred.

**Design Assumptions:**

- All active arms share the same information level in Phase 2.

- Exactly one active arm is selected at the end of Phase 2 based on the
  largest observed Z-statistic when the trial continues to Phase 3.

## References

Ping Gao, Yingqiu Li. Adaptive two-stage seamless sequential design for
clinical trials. Journal of Biopharmaceutical Statistics, 2025, 35(4),
565-587.

## Author

Kaifeng Lu, <kaifenglu@gmail.com>

## Examples

``` r
# Setup: 2 active arms vs control in Phase 2; 1 selected arm vs control
# in Phase 3. Phase 3 has 2 sequential looks.

# Information levels: equal spacing over 3 looks based on a maximum of
# 110 patients per arm, SD = 1.0
I <- c(110 / (2 * 1.0^2) * seq(1, 3)/3)

# O'Brien-Fleming efficacy boundaries
b <- c(3.776605, 2.670463, 2.180424)

# No futility stopping
p0 <- exitprob_seamless(M = 2, theta = c(0, 0), K = 2, b = b, I = I)
cumsum(p0$exitProbUpper)
#> [1] 0.0001572756 0.0066431322 0.0250000060

# Add futility stopping
a <- c(0, 0.5, b[3])
p1 <- exitprob_seamless(M = 2, theta = c(0.3, 0.5), K = 2, b = b, a = a, I = I)
cbind(
  cumulativeEfficacy = cumsum(p1$exitProbUpper),
  cumulativeFutility = cumsum(p1$exitProbLower)
)
#>      cumulativeEfficacy cumulativeFutility
#> [1,]         0.05477567        0.007803144
#> [2,]         0.62292767        0.014040545
#> [3,]         0.89800885        0.101991188
```
