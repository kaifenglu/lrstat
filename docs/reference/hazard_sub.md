# Hazard Function for Sub Population

Computes the hazard function of a piecewise exponential distribution for
the biomarker negative sub population, such that the resulting survival
function for the ITT population closely matches a given piecewise
survival function.

## Usage

``` r
hazard_sub(
  piecewiseSurvivalTime = 0L,
  hazard_itt = NA_real_,
  hazard_pos = NA_real_,
  p_pos = NA_real_
)
```

## Arguments

- piecewiseSurvivalTime:

  A vector that specifies the starting time of piecewise exponential
  survival time intervals. Must start with 0, e.g., `c(0, 6)` breaks the
  time axis into 2 event intervals: \\\[0, 6)\\ and \\\[6, \infty)\\.
  Defaults to 0 for exponential distribution.

- hazard_itt:

  A scalar or numeric vector specifying the hazard(s) for the ITT
  population based on a piecewise exponential distribution.

- hazard_pos:

  A scalar or numeric vector specifying the hazard(s) for the biomarker
  positive sub population based on a piecewise exponential distribution.

- p_pos:

  A numeric value specifying the prevalence of the biomarker positive
  sub population.

## Value

A list with the following components:

- `piecewiseSurvivalTime`: A vector that specifies the starting time
  points of the intervals for the piecewise exponential distribution for
  the biomarker negative sub population.

- `hazard_pos`: A numeric vector representing the hazard rates for the
  piecewise exponential distribution of the biomarker positive sub
  population at the same time points as the biomarker negative sub
  population.

- `hazard_neg`: A numeric vector representing the estimated hazard rates
  for the piecewise exponential distribution of the biomarker negative
  sub population.

- `p_pos`: The prevalence of the biomarker positive sub population (as
  input).

## Details

This function determines the hazard vector \\\lambda\_{\text{neg}}\\ for
the piecewise exponential distribution of the biomarker negative sub
population, so that the implied survival function for the ITT population
closely matches the specified piecewise exponential distribution with
hazard vector \\\lambda\_{\text{itt}}\\.

Let \\p\_{\text{pos}}\\ be the prevalence of the biomarker positive sub
population, then the survival function for the ITT population is given
by \$\$S\_{\text{itt}}(t) = p\_{\text{pos}} S\_{\text{pos}}(t) + (1 -
p\_{\text{pos}}) S\_{\text{neg}}(t)\$\$ where \\S\_{\text{pos}}(t)\\ and
\\S\_{\text{neg}}(t)\\ are the survival functions for the biomarker
positive and biomarker negative sub populations, respectively.

Matching is performed sequentially at the internal cutpoints \\u_2, ...,
u_J\\ and at the point \\u_J + \log(2)/\lambda\_{\text{itt},J}\\ for the
final interval, as well as the percentile points at 10%, 20%, ..., 90%,
and 95%, to solve for \\\lambda\_{\text{neg},1}, \ldots,
\lambda\_{\text{neg},K}\\, where \\K\\ is the total number of unique cut
points.

## Author

Kaifeng Lu (<kaifenglu@gmail.com>)

## Examples

``` r
u <- c(0, 1, 3, 4)
lambda_itt <- c(0.0151, 0.0403, 0.0501, 0.0558)
lambda_pos <- c(0.0115, 0.0302, 0.0351, 0.0404)
p_pos <- 0.3
hazard_sub(u, lambda_itt, lambda_pos, p_pos)
#> $piecewiseSurvivalTime
#>  [1]  0.000000  1.000000  3.000000  3.192825  4.000000  5.386085  7.779121
#>  [8] 10.541678 13.809089 16.421992 17.808078 22.963670 30.230070 42.652063
#> 
#> $hazard_pos
#>  [1] 0.0115 0.0302 0.0351 0.0351 0.0404 0.0404 0.0404 0.0404 0.0404 0.0404
#> [11] 0.0404 0.0404 0.0404 0.0404
#> 
#> $hazard_neg
#>  [1] 0.01664683 0.04471458 0.05676603 0.05683882 0.06288919 0.06319981
#>  [7] 0.06365275 0.06423148 0.06485005 0.06530670 0.06613675 0.06801176
#> [13] 0.07224688 0.08162772
#> 
#> $p_pos
#> [1] 0.3
#> 
```
