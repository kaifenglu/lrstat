# Assess Proportional Hazards Assumption Based on Scaled Schoenfeld Residuals

Obtains the scaled Schoenfeld residuals and tests the proportional
hazards assumption using a score test for the interaction between each
covariate and a transformed time variable.

## Usage

``` r
zph_phregr(object, transform = "km")
```

## Arguments

- object:

  The output from the `phregr` call.

- transform:

  A character string indicating how survival times should be transformed
  before the test is performed. Supported values include "identity",
  "log", "rank", and "km" (default).

## Value

A list with the following components:

- `table` A matrix with one row for each parameter and a final row for
  the global test. The columns contain the score test for adding the
  time-dependent term, the degrees of freedom, and the two-sided
  p-value.

- `x` The transformed time values.

- `time` The original (untransformed) event times, with tied event times
  repeated.

- `strata` The stratum index for each event.

- `y` The matrix of scaled Schoenfeld residuals, with one column for
  each parameter and one row for each event. Column names correspond to
  the parameter names.

- `var` An approximate covariance matrix of the scaled Schoenfeld
  residuals, used to construct an approximate standard error band for
  plots.

- `transform` the transformation applied to the time values.

## Details

This corresponds to the `cox.zph` function from the `survival` package
with `terms = FALSE` and `global = TRUE`.

## References

Patricia M. Grambsch and Terry M. Therneau. Proportional hazards tests
and diagnostics based on weighted residuals. Biometrika 1994; 81:515-26.

## Author

Kaifeng Lu, <kaifenglu@gmail.com>

## Examples

``` r
fit <- phregr(data = liver, time = "Time", event = "Status",
              covariates = c("log(Bilirubin)", "log(Protime)",
                             "log(Albumin)", "Age", "Edema"),
              ties = "breslow")

zph <- zph_phregr(fit, transform = "log")

zph$table
#>                      chisq df           p
#> log(Bilirubin)  0.15583626  1 0.693019270
#> log(Protime)   10.45818301  1 0.001221073
#> log(Albumin)    1.79397334  1 0.180442799
#> Age             0.03514403  1 0.851294182
#> Edema           2.86147789  1 0.090724109
#> GLOBAL         13.32206571  5 0.020540524
```
