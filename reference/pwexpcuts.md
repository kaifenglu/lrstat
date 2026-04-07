# Piecewise Exponential Approximation to a Survival Distribution

Obtains the piecewise exponential distribution that approximates a
survival distribution.

## Usage

``` r
pwexpcuts(S, ..., tol = 1e-04)
```

## Arguments

- S:

  The survival function of a univariate survival time.

- ...:

  Additional arguments to be passed to S.

- tol:

  The tolerance for convergence of the profile log-likelihood. Defaults
  to 0.0001.

## Value

A list with three components:

- `piecewiseSurvivalTime`: A vector that specifies the starting time of
  piecewise exponential survival time intervals. Must start with 0,
  e.g., c(0, 6) breaks the time axis into 2 event intervals: \[0, 6) and
  \[6, Inf).

- `lambda`: A vector of hazard rates for the event. One for each
  analysis time interval.

- `loglik`: The sequence of the asymptotic limit of the piecewise
  exponential log-likelihood for an increasing number of change points.

## Details

This function computes the piecewise exponential approximation to a
survival distribution. The piecewise exponential model divides the time
axis into \\J\\ intervals defined by the change points, where each
interval \\\[t_j, t\_{j+1})\\ has a constant hazard rate \\\lambda_j\\.
The time intervals are specified as: \$\$\[t_1, t_2), \[t_2, t_3),
\ldots, \[t\_{J}, t\_{J+1})\$\$ where \\t_1 = 0\\, \\t\_{J+1} =
\infty\\, and \\t_j = \tau\_{j-1}\\ for \\j = 2, \ldots, J\\. The
function starts with \\J = 2\\ (1 change point) and gradually increases
\\J\\ by adding one change point at a time until the maximized profile
log-likelihood for change points stabilizes, i.e., the relative increase
in the maximum of the profile log-likelihood function is less than
`tol`. If the relative change in the hazard rate is also less than
`tol`, the function stops and returns the results.

## Author

Kaifeng Lu, <kaifenglu@gmail.com>

## Examples

``` r
# Example 1: Piecewise exponential
pwexpcuts(ptpwexp, piecewiseSurvivalTime = c(0, 3.4, 5.5),
          lambda = c(0.0168, 0.0833, 0.0431), lowerBound = 0,
          lower.tail = FALSE)
#> $piecewiseSurvivalTime
#> [1] 0.000000 3.399971 5.500012
#> 
#> $lambda
#> [1] 0.01680000 0.08329879 0.04310000
#> 
#> $loglik
#> [1] -4.120403 -4.096664
#> 

# Example 2: Weibull
pwexpcuts(pweibull, shape = 1.37, scale = 1/0.818, lower.tail = FALSE)
#> $piecewiseSurvivalTime
#> [1] 0.0000000 0.5266367 1.2814623 2.0762636 2.8645205 3.6338592
#> 
#> $lambda
#> [1] 0.5842842 0.9756606 1.2392703 1.4363304 1.5939072 1.7628962
#> 
#> $loglik
#> [1] -1.066978 -1.058643 -1.057256 -1.057019 -1.056978
#> 
```
