# Correlation Between PFS and OS Given Correlation Between PD and OS

Computes the correlation between PFS and OS given the correlation
between PD and OS.

## Usage

``` r
corr_pfs_os(
  piecewiseSurvivalTime = 0L,
  hazard_pfs = NA_real_,
  hazard_os = NA_real_,
  rho_pd_os = NA_real_
)
```

## Arguments

- piecewiseSurvivalTime:

  A vector that specifies the starting time of piecewise exponential
  survival time intervals. Must start with 0, e.g., `c(0, 6)` breaks the
  time axis into 2 event intervals: \\\[0, 6)\\ and \\\[6, \infty)\\.
  Defaults to 0 for exponential distribution.

- hazard_pfs:

  A scalar or numeric vector specifying the hazard(s) for PFS based on a
  piecewise exponential distribution.

- hazard_os:

  A scalar or numeric vector specifying the hazard(s) for overall
  survival (OS) based on a piecewise exponential distribution.

- rho_pd_os:

  A numeric value specifying the correlation between PD and OS times.

## Value

The estimated correlation between PFS and OS.

## Details

This function first determines the piecewise exponential distribution
for PD such that the implied survival function for PFS time,
\\T\_{\text{pfs}} = \min(T\_{\text{pd}}, T\_{\text{os}})\\, closely
matches the specified piecewise exponential distribution for PFS with
hazard vector \\\lambda\_{\text{pfs}}\\. Then, it calculates the
correlation between PFS and OS times based on the derived piecewise
exponential distribution for PD and the given piecewise exponential
distribution for OS.

## Author

Kaifeng Lu (<kaifenglu@gmail.com>)

## Examples

``` r
u <- c(0, 1, 3, 4)
lambda1 <- c(0.0151, 0.0403, 0.0501, 0.0558)
lambda2 <- 0.0145
rho_pd_os <- 0.5
corr_pfs_os(u, lambda1, lambda2, rho_pd_os)
#> [1] 0.5483736
```
