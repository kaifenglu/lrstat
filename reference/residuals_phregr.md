# Residuals for Proportional Hazards Regression Models

Obtains the martingale, deviance, score, or Schoenfeld residuals for a
proportional hazards regression model.

## Usage

``` r
residuals_phregr(
  object,
  type = c("martingale", "deviance", "score", "schoenfeld", "dfbeta", "dfbetas",
    "scaledsch"),
  collapse = FALSE,
  weighted = (type %in% c("dfbeta", "dfbetas"))
)
```

## Arguments

- object:

  The output from the `phregr` call.

- type:

  The type of residuals desired, with options including `"martingale"`,
  `"deviance"`, `"score"`, `"schoenfeld"`, `"dfbeta"`, `"dfbetas"`, and
  `"scaledsch"`.

- collapse:

  Whether to collapse the residuals by `id`. This is not applicable for
  Schoenfeld type residuals.

- weighted:

  Whether to compute weighted residuals.

## Value

For martingale and deviance residuals, the result is a vector with one
element corresponding to each subject (without `collapse`). For score
residuals, the result is a matrix where each row represents a subject
and each column corresponds to a variable. The row order aligns with the
input data used in the original fit. For Schoenfeld residuals, the
result is a matrix with one row for each event and one column per
variable. These rows are sorted by time within strata, with the
attributes `stratum` and `time` included.

Score residuals represent each individual's contribution to the score
vector. Two commonly used transformations of this are `dfbeta`, which
represents the approximate change in the coefficient vector if the
observation is excluded, and `dfbetas`, which gives the approximate
change in the coefficients scaled by the standard error of the
coefficients.

## Details

For score and Schoenfeld type residuals, the proportional hazards model
must include at least one covariate. The algorithms for `deviance`,
`dfbeta`, `dfbetas`, and `scaledsch` residuals follow the
`residuals.coxph` function in the `survival` package.

## References

Terry M. Therneau, Patricia M. Grambsch, and Thomas M. Fleming.
Martingale based residuals for survival models. Biometrika 1990;
77:147-160.

Patricia M. Grambsch and Terry M. Therneau. Proportional hazards tests
and diagnostics based on weighted residuals. Biometrika 1994; 81:515-26.

## Author

Kaifeng Lu, <kaifenglu@gmail.com>

## Examples

``` r
library(dplyr)

# Example 1 with right-censored data
fit1 <- phregr(data = rawdata %>% filter(iterationNumber == 1) %>%
                 mutate(treat = 1*(treatmentGroup == 1)),
               stratum = "stratum",
               time = "timeUnderObservation", event = "event",
               covariates = "treat")

ressco <- residuals_phregr(fit1, type = "score")
head(ressco)
#> [1]  0.2344213 -0.2232656  0.1964093  0.5382246  0.3367144 -0.2268484

# Example 2 with counting process data
fit2 <- phregr(data = heart %>% mutate(rx = as.numeric(transplant) - 1),
               time = "start", time2 = "stop", event = "event",
               covariates = c("rx", "age"), id = "id", robust = TRUE)

resssch <- residuals_phregr(fit2, type = "scaledsch")
head(resssch)
#>              rx         age
#> [1,] -0.3762107  0.12783993
#> [2,]  0.1024506 -0.03363229
#> [3,] -0.4513918  0.11091644
#> [4,] -0.4282462  0.10487560
#> [5,] -0.8083496  0.14573692
#> [6,]  0.2879278 -0.14038327
```
