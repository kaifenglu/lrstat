# Efficacy Boundaries for Group Sequential Design

Obtains the efficacy stopping boundaries for a group sequential design.

## Usage

``` r
getBound(
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

- k:

  Look number for the current analysis.

- informationRates:

  Information rates up to the current look. Must be increasing and less
  than or equal to 1.

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
  for `"WT"`, \\rho\\ for `"sfKD"`, and \\\gamma\\ for `"sfHSD"`.

- userAlphaSpending:

  The user defined alpha spending. Cumulative alpha spent up to each
  stage.

- spendingTime:

  A vector of length `k` for the error spending time at each analysis.
  Must be increasing and less than or equal to 1. Defaults to missing,
  in which case, it is the same as `informationRates`.

- efficacyStopping:

  Indicators of whether efficacy stopping is allowed at each stage.
  Defaults to `TRUE` if left unspecified.

## Value

A numeric vector of critical values up to the current look.

## Details

If `typeAlphaSpending` is "OF", "P", or "WT", then the boundaries will
be based on equally spaced looks.

## Author

Kaifeng Lu, <kaifenglu@gmail.com>

## Examples

``` r
getBound(k = 2, informationRates = c(0.5,1),
         alpha = 0.025, typeAlphaSpending = "sfOF")
#> [1] 2.962588 1.968596
```
