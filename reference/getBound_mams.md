# Efficacy Boundaries for Multi-Arm Multi-Stage Design

Calculates the efficacy stopping boundaries for a multiple comparison
sequential design.

## Usage

``` r
getBound_mams(
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

  Number of active treatment arms.

- r:

  Randomization ratio of each active arm to the common control.

- corr_known:

  Logical. If `TRUE`, the correlation between Wald statistics is derived
  from the randomization ratio \\r\\ as \\r / (r + 1)\\. If `FALSE`, a
  conservative correlation of 0 is assumed.

- k:

  The index of the current look.

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

  A numeric vector of length \\k\\ specifying the error spending time at
  each analysis. Values must be strictly increasing and \\\le 1\\. If
  omitted, defaults to `informationRates`.

- efficacyStopping:

  Indicators of whether efficacy stopping is allowed at each stage.
  Defaults to `TRUE` if left unspecified.

## Value

A numeric vector of length \\k\\ containing the critical values (on the
standard normal Z-scale) for each analysis up to the current look.

## Details

The function determines critical values by solving for the boundary that
satisfies the alpha-spending requirement, given the selection of the
"best" arm at the end of Phase 2.

If `typeAlphaSpending` is specified as `"OF"` (O'Brien-Fleming), `"P"`
(Pocock), or `"WT"` (Wang-Tsiatis), the boundaries are calculated
assuming the looks are equally spaced in terms of information.

## Author

Kaifeng Lu, <kaifenglu@gmail.com>

## Examples

``` r
# Determine O'Brien-Fleming boundaries for a TSSSD with
# 2 active arms and 3 looks.
getBound_mams(M = 2, k = 3, informationRates = seq(1, 3)/3,
              alpha = 0.025, typeAlphaSpending = "OF")
#> [1] 3.886562 2.748214 2.243907
```
