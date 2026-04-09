# Equicoordinate Quantiles of the Multivariate Normal Distribution

Computes the equicoordinate quantile \\q\\ such that \\P(X_1 \le q, X_2
\le q, \ldots, X_k \le q) = p\\ for a multivariate normal random vector
\\X\\.

## Usage

``` r
qmvnormr(
  p,
  mean = NULL,
  sigma,
  pivot = FALSE,
  fast = TRUE,
  n0 = 1024,
  n_max = 16384,
  R = 8,
  abseps = 1e-04,
  releps = 0,
  seed = 0,
  parallel = TRUE,
  nthreads = 0
)
```

## Arguments

- p:

  The probability level (cumulative probability).

- mean:

  The mean vector. If `NULL` (default), a zero vector of appropriate
  length is used.

- sigma:

  The covariance (or correlation) matrix of the distribution.

- pivot:

  Logical; if `TRUE`, applies an initial pivoting step to reorder the
  integration variables for improved efficiency.

- fast:

  Logical; if `TRUE`, uses a fast approximation of the univariate normal
  CDF and quantile functions.

- n0:

  Initial number of samples per replication for the Monte Carlo
  integration.

- n_max:

  Maximum number of samples allowed per replication.

- R:

  Number of independent replications used to estimate the error.

- abseps:

  Absolute error tolerance for the probability calculation.

- releps:

  Relative error tolerance for the probability calculation.

- seed:

  Random seed for reproducibility. If 0, a seed is generated from the
  computer clock.

- parallel:

  Logical; if `TRUE`, computations are performed in parallel.

- nthreads:

  Number of threads for parallel execution. If 0, the default
  RcppParallel behavior is used.

## Value

A numeric value representing the calculated equicoordinate quantile.

## Details

This function finds the value \\q\\ using a root-finding algorithm
applied to the `pmvnormr` function. It solves for the value where the
multivariate normal cumulative distribution function equals the target
probability \\p\\.

## Author

Kaifeng Lu, <kaifenglu@gmail.com>

## Examples

``` r
n <- 5
mean <- rep(0, n)
sigma <- matrix(0.5, n, n)
diag(sigma) <- 1
qmvnormr(0.5, mean = mean, sigma = sigma)
#> [1] 0.8150347
```
