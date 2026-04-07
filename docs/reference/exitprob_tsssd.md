# Exit Probabilities for Two-Stage Seamless Sequential Design (TSSSD)

Computes the exit (rejection) probabilities for a two-stage selection
and testing design. In Phase 2, multiple active arms are compared
against a common control arm. The best-performing arm is selected to
proceed to Phase 3, where it is tested against the common control over
multiple looks.

## Usage

``` r
exitprob_tsssd(
  M = NA_integer_,
  r = 1,
  theta = NA_real_,
  corr_known = TRUE,
  K = NA_integer_,
  b = NA_real_,
  I = NA_real_
)
```

## Arguments

- M:

  Number of active treatment arms in Phase 2 (\\M \ge 2\\).

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

  A vector of critical values (length \\K+1\\). The first element is for
  Phase 2; the remaining \\K\\ elements are for the looks in Phase 3.

- I:

  A vector of information levels (length \\K+1\\) for any active arm
  versus the common control. The first element is for Phase 2; the
  remaining \\K\\ elements are for the looks in Phase 3.

## Value

A list containing:

- `exitProb`: A vector of length \\K + 1\\. The first element is the
  probability of rejection in Phase 2; the remaining elements are the
  probabilities of rejection at each look in Phase 3.

- `exitProbByArm`: A \\(K+1) \times M\\ matrix. The \\(k, m)\\-th entry
  represents the probability that the global null is rejected at look
  \\k\\ given that arm \\m\\ was selected as the best arm.

- `selectAsBest`: A vector of length \\M\\ containing the probability
  that each active arm is selected to move on to Phase 3.

## Details

The function assumes a multivariate normal distribution for the Wald
statistics. The "best" arm is defined as the active arm with the largest
score statistic at the end of Phase 2.

**Decision Rules:**

- **Phase 2**: The global null hypothesis is rejected if the Wald
  statistic for the best arm, \\Z(I_0)\\, satisfies \\Z(I_0) \ge b_0\\.

- **Phase 3**: If the trial continues, the hypothesis is rejected at
  look \\k\\ if \\Z(I_k) \ge b_k\\ and all previous looks (including
  Phase 2) failed to reject.

**Design Assumptions:**

- All active arms share the same information level in Phase 2.

- Exactly one active arm is selected at the end of Phase 2 based on the
  largest observed statistic.

## Author

Kaifeng Lu, <kaifenglu@gmail.com>

## Examples

``` r
# Setup: 2 active arms vs control in phase 2; 1 selected arm vs control
# in phase 3. Phase 3 has 2 sequential looks.

# Information levels: equal spacing over 3 looks based on max 110 patients
# per arm, SD = 1.0
I <- c(110 / (2 * 1.0^2) * seq(1, 3)/3)

# O'Brien-Fleming critical values
b <- c(3.776606, 2.670463, 2.180424)

# Type I error under the global null hypothesis
p0 <- exitprob_tsssd(M = 2, theta = c(0, 0), K = 2, b = b, I = I)
cumsum(p0$exitProb)
#> [1] 0.000157278 0.006643135 0.025000008

# Power under alternative: Treatment effects of 0.3 and 0.5
p1 <- exitprob_tsssd(M = 2, theta = c(0.3, 0.5), K = 2, b = b, I = I)
cumsum(p1$exitProb)
#> [1] 0.05477555 0.62309680 0.90160747
```
