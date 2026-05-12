# Efficacy Boundaries for Phase 2/3 Seamless Design

Calculates the efficacy stopping boundaries for a phase 2/3 seamless
design, accounting for the selection of the best arm at the end of Phase
2 and sequential testing in Phase 3.

## Usage

``` r
getBound_seamless(
  M = NA_integer_,
  r = 1,
  corr_known = TRUE,
  k = NA_integer_,
  informationRates = NA_real_,
  alpha = 0.025,
  typeAlphaSpending = "sfOF",
  parameterAlphaSpending = NA_real_,
  userAlphaSpending = NA_real_,
  spendingTime = NA_real_,
  efficacyStopping = NA_integer_
)
```

## Arguments

- M:

  Number of active treatment arms in Phase 2.

- r:

  Randomization ratio of each active arm to the common control in Phase
  2.

- corr_known:

  Logical. If `TRUE`, the correlation between Wald statistics in Phase 2
  is derived from the randomization ratio \\r\\ as \\r / (r + 1)\\. If
  `FALSE`, a conservative correlation of 0 is assumed.

- k:

  The index of the current look in Phase 3.

- informationRates:

  A numeric vector of information rates up to the current look. Values
  must be strictly increasing and \\\le 1\\.

- alpha:

  The significance level. Defaults to 0.025.

- typeAlphaSpending:

  The type of alpha spending. One of the following: `"OF"` for
  O'Brien-Fleming boundaries, `"P"` for Pocock boundaries, `"WT"` for
  Wang & Tsiatis boundaries, `"sfOF"` for O'Brien-Fleming type spending
  function, `"sfP"` for Pocock type spending function, `"sfKD"` for Kim
  & DeMets spending function, `"sfHSD"` for Hwang, Shi & DeCani spending
  function, `"user"` for user defined spending, and `"none"` for no
  early efficacy stopping. Defaults to `"sfOF"`.

- parameterAlphaSpending:

  The parameter value for the alpha spending. Corresponds to \\\Delta\\
  for `"WT"`, \\\rho\\ for `"sfKD"`, and \\\gamma\\ for `"sfHSD"`.

- userAlphaSpending:

  The user defined alpha spending. Cumulative alpha spent up to each
  stage.

- spendingTime:

  A numeric vector of length \\k+1\\ specifying the error spending time
  at each analysis. Values must be strictly increasing and \\\le 1\\. If
  omitted, defaults to `informationRates`.

- efficacyStopping:

  Indicators of whether efficacy stopping is allowed at each stage.
  Defaults to `TRUE` if left unspecified.

## Value

A numeric vector of length \\k + 1\\ containing the critical values (on
the standard normal Z-scale) for each analysis up to the current look.

## Details

The function determines critical values by solving for the boundary that
satisfies the alpha-spending requirement, given the selection of the
"best" arm at the end of Phase 2.

If `typeAlphaSpending` is `"OF"`, `"P"`, `"WT"`, or `"none"`, then
`informationRates`, `efficacyStopping`, and `spendingTime` must be of
full length `kMax`, and `informationRates` and `spendingTime` must end
with 1.

## References

Ping Gao, Yingqiu Li. Adaptive two-stage seamless sequential design for
clinical trials. Journal of Biopharmaceutical Statistics, 2025, 35(4),
565-587.

## Author

Kaifeng Lu, <kaifenglu@gmail.com>

## Examples

``` r
# Determine O'Brien-Fleming boundaries for a seamless design with
# 2 active arms in Phase 2 and 2 looks in Phase 3 (3 looks total).
getBound_seamless(M = 2, k = 2, informationRates = seq(1, 3)/3,
                  alpha = 0.025, typeAlphaSpending = "OF")
#> [1] 3.776605 2.670463 2.180424
```
