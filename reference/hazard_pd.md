# Hazard Function for Progressive Disease (PD) Given Correlation Between PD and OS

Computes the hazard function of a piecewise exponential distribution for
progressive disease (PD), such that the resulting hazard function for
progression-free survival (PFS) closely matches a given piecewise hazard
for PFS.

## Usage

``` r
hazard_pd(
  piecewiseSurvivalTime = 0L,
  hazard_pfs = NA_real_,
  hazard_os = NA_real_,
  rho_pd_os = 0.5
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

A list with the following components:

- `piecewiseSurvivalTime`: A vector that specifies the starting time
  points of the intervals for the piecewise exponential distribution for
  PD.

- `hazard_pd`: A numeric vector representing the calculated hazard rates
  for the piecewise exponential distribution of PD.

- `hazard_os`: A numeric vector representing the hazard rates for the
  piecewise exponential distribution of OS at the same time points as
  PD.

- `rho_pd_os`: The correlation between PD and OS times (as input).

## Details

This function determines the hazard vector \\\lambda\_{\text{pd}}\\ for
the piecewise exponential distribution of PD, so that the implied
survival function for PFS time, \\T\_{\text{pfs}} = \min(T\_{\text{pd}},
T\_{\text{os}})\\, closely matches the specified piecewise exponential
distribution for PFS with hazard vector \\\lambda\_{\text{pfs}}\\.

To achieve this, we simulate \\(Z\_{\text{pd}}, Z\_{\text{os}})\\ from a
standard bivariate normal distribution with correlation \\\rho\\. Then,
\\U\_{\text{pd}} = \Phi(Z\_{\text{pd}})\\ and \\U\_{\text{os}} =
\Phi(Z\_{\text{os}})\\ are generated, where \\\Phi\\ denotes the
standard normal CDF.

The times to PD and OS are obtained via the inverse transform method
using quantile functions of the piecewise exponential distribution:
\$\$T\_{\text{pd}} =
\text{qpwexp}(U\_{\text{pd}},u,\lambda\_{\text{pd}})\$\$
\$\$T\_{\text{os}} =
\text{qpwexp}(U\_{\text{os}},u,\lambda\_{\text{os}})\$\$ where
`u = piecewiseSurvivalTime`.

The function solves for \\\lambda\_{\text{pd}}\\ such that the survival
function of \\T\_{\text{pfs}}\\ closely matches that of a piecewise
exponential distribution with hazard \\\lambda\_{\text{pfs}}\\:
\$\$P(\min(T\_{\text{pd}}, T\_{\text{os}}) \> t) =
S\_{\text{pfs}}(t)\$\$ Since \$\$Z\_{\text{pd}} =
\Phi^{-1}(\text{ppwexp}(T\_\text{pd}, u, \lambda\_{\text{pd}}))\$\$ and
\$\$Z\_{\text{os}} = \Phi^{-1}(\text{ppwexp}(T\_\text{os}, u,
\lambda\_{\text{os}}))\$\$ we have \$\$P(\min(T\_{\text{pd}},
T\_{\text{os}}) \> t) = P(Z\_{\text{pd}} \>
\Phi^{-1}(\text{ppwexp}(t,u,\lambda\_{\text{pd}})), Z\_{\text{os}} \>
\Phi^{-1}(\text{ppwexp}(t,u,\lambda\_{\text{os}})))\$\$ while
\$\$S\_{\text{pfs}}(t) = 1 -
\text{ppwexp}(t,u,\lambda\_{\text{pfs}})\$\$

Matching is performed sequentially at the internal cut points \\u_2,
..., u_J\\ and at the point \\u_J + \log(2)/\lambda\_{\text{pfs},J}\\
for the final interval, as well as the percentile points at 10%, 20%,
..., 90%, and 95% to solve for \\\lambda\_{\text{pd},1}, \ldots,
\lambda\_{\text{pd},K}\\, where \\K\\ is the total number of unique cut
points.

## Author

Kaifeng Lu (<kaifenglu@gmail.com>)

## Examples

``` r
u <- c(0, 1, 3, 4)
lambda1 <- c(0.0151, 0.0403, 0.0501, 0.0558)
lambda2 <- 0.0145
rho_pd_os <- 0.5
hazard_pd(u, lambda1, lambda2, rho_pd_os)
#> $piecewiseSurvivalTime
#>  [1]  0.000000  1.000000  3.000000  3.192825  4.000000  5.386085  7.779121
#>  [8] 10.541678 13.809089 16.421992 17.808078 22.963670 30.230070 42.652063
#> 
#> $hazard_pd
#>  [1] 0.0008356782 0.0311408296 0.0427651037 0.0430302691 0.0498781129
#>  [6] 0.0505667622 0.0512705251 0.0518601620 0.0522914905 0.0525291948
#> [11] 0.0528407297 0.0532889452 0.0537593532 0.0541462469
#> 
#> $hazard_os
#>  [1] 0.0145 0.0145 0.0145 0.0145 0.0145 0.0145 0.0145 0.0145 0.0145 0.0145
#> [11] 0.0145 0.0145 0.0145 0.0145
#> 
#> $rho_pd_os
#> [1] 0.5
#> 
```
