# Profile Log-Likelihood Function for Change Points in Piecewise Exponential Approximation

Obtains the profile log-likelihood function for change points in the
piecewise exponential approximation to a survival function.

## Usage

``` r
pwexploglik(tau, S, ...)
```

## Arguments

- tau:

  The numeric vector of change points.

- S:

  The survival function of a univariate survival time.

- ...:

  Additional arguments to be passed to S.

## Value

A list with the following three components:

- `piecewiseSurvivalTime`: A vector that specifies the starting time of
  piecewise exponential survival time intervals.

- `lambda`: A vector of hazard rates for the event. One for each
  analysis time interval.

- `loglik`: The value of the profile log-likelihood.

## Details

This function computes the profile log-likelihood for change points in a
piecewise exponential survival model.

Let \\S(t)\\ denote the survival function of a univariate survival time,
and \\\tau\\ be a vector of \\J-1\\ change points. The piecewise
exponential survival model divides the time axis into \\J\\ intervals
defined by the change points \\\tau\\, where each interval \\\[t_j,
t\_{j+1})\\ has a constant hazard rate \\\lambda_j\\. The time intervals
are specified as: \$\$\[t_1, t_2), \[t_2, t_3), \ldots, \[t\_{J},
t\_{J+1})\$\$ where \\t_1 = 0\\, \\t\_{J+1} = \infty\\, and \\t_j =
\tau\_{j-1}\\ for \\j = 2, \ldots, J\\.

For each subject, the expected number of events occurring in the
\\j\\-th interval is \$\$d_j = E\\I(t_j \< Y \leq t\_{j+1})\\ = S(t_j) -
S(t\_{j+1})\$\$ The expected exposure in the \\j\\-th interval is:
\$\$e_j = E\\(Y-t_j)I(t_j \< Y \leq t\_{j+1}) + (t\_{j+1} - t_j)I(Y \>
t\_{j+1})\\\$\$ which can be shown to be equivalent to \$\$e_j =
\int\_{t_j}^{t\_{j+1}} S(t) dt\$\$

The log-likelihood for the piecewise exponential model is:
\$\$\ell(\tau,\lambda) = \sum\_{j=1}^J \\d_j \log(\lambda_j) - e_j
\lambda_j\\\$\$ The profile log-likelihood for \\\tau\\ is obtained by
maximizing \\\ell(\tau,\lambda)\\ with respect to \\\lambda\\ for fixed
\\\tau\\. The maximum likelihood estimate of the hazard rate in the
\\j\\-th interval is \$\$\lambda_j = \frac{d_j}{e_j}\$\$ Substituting
back, the profile log-likelihood is \$\$\ell(\tau) = \sum\_{j=1}^J d_j
\log(d_j/e_j) - 1\$\$ where we use the fact that \\\sum\_{j=1}^J d_j =
1\\.

## Author

Kaifeng Lu, <kaifenglu@gmail.com>

## Examples

``` r
pwexploglik(tau = c(0.5, 1.2, 2.8), pweibull,
            shape = 1.37, scale = 1/0.818, lower.tail = FALSE)
#> $piecewiseSurvivalTime
#> [1] 0.0 0.5 1.2 2.8
#> 
#> $lambda
#> [1] 0.5741586 0.9558528 1.2688072 1.6284728
#> 
#> $loglik
#> [1] -1.05696
#> 
```
