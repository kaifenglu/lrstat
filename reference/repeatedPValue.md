# Repeated p-Values for Group Sequential Design

Obtains the repeated p-values for a group sequential design.

## Usage

``` r
repeatedPValue(
  kMax,
  typeAlphaSpending = "sfOF",
  parameterAlphaSpending = NA,
  maxInformation = 1,
  p,
  information,
  spendingTime = NULL,
  nthreads = 0
)
```

## Arguments

- kMax:

  The maximum number of stages.

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

- maxInformation:

  The target maximum information. Defaults to 1, in which case,
  `information` represents `informationRates`.

- p:

  The raw p-values at look 1 to look `k`. It can be a matrix with `k`
  columns for `k <= kMax`.

- information:

  The observed information by look. It can be a matrix with `k` columns.

- spendingTime:

  The error spending time at each analysis, must be increasing and less
  than or equal to 1. Defaults to `NULL`, in which case, it is the same
  as `informationRates` derived from `information` and `maxInformation`.
  It can be a matrix with `k` columns.

- nthreads:

  The number of threads to use in simulations (0 means the default
  RcppParallel behavior).

## Value

The repeated p-values at look 1 to look `k`.

## Author

Kaifeng Lu, <kaifenglu@gmail.com>

## Examples

``` r
# Example 1: informationRates different from spendingTime
repeatedPValue(kMax = 3, typeAlphaSpending = "sfOF",
               maxInformation = 800,
               p = c(0.2, 0.15, 0.1),
               information = c(529, 700, 800),
               spendingTime = c(0.6271186, 0.8305085, 1))
#> [1] 0.3101673 0.2258992 0.1232186

# Example 2: Maurer & Bretz (2013), current look is not the last look
repeatedPValue(kMax = 3, typeAlphaSpending = "sfOF",
               p = matrix(c(0.0062, 0.017,
                            0.009, 0.13,
                            0.0002, 0.0035,
                            0.002, 0.06),
                          nrow=4, ncol=2),
               information = c(1/3, 2/3),
               nthreads = 1)
#>           [,1]        [,2]
#> [1,] 0.1140577 0.002393359
#> [2,] 0.1682137 0.017160592
#> [3,] 0.1315366 0.011648306
#> [4,] 0.3820275 0.128460197
```
