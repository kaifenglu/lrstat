# Multivariate Normal Probabilities over a Hyperrectangle

Computes the probability that a multivariate normal random vector falls
within a rectangular region defined by lower and upper bounds.

## Usage

``` r
pmvnormr(
  lower,
  upper,
  mean = NULL,
  sigma,
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

- lower:

  A numeric vector of lower integration limits.

- upper:

  A numeric vector of upper integration limits.

- mean:

  The mean vector. If `NULL` (default), a zero vector of appropriate
  length is used.

- sigma:

  The covariance (or correlation) matrix of the distribution.

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

The estimated probability with the following attributes:

- `method`: `"exact"` for analytic methods or `"qmc"` for the Monte
  Carlo approach.

- `error`: The estimated error (half-width of 95% confidence interval).

- `nsamples`: The total number of samples used across all replications.

## Details

The function automatically selects the most efficient computation method
based on the input:

- Analytic Methods: Used for univariate cases or multivariate
  distributions with a compound symmetry correlation structure and
  non-negative correlations.

- Monte Carlo Estimation: Used for all other cases. The algorithm
  employs a randomized quasi-Monte Carlo (QMC) approach using
  generalized Halton sequences.

To improve efficiency and accuracy, the QMC approach incorporates:

- Sequential Conditioning: Mimics the standardization and transformation
  approach used in
  [`mvtnorm::lpmvnorm`](https://rdrr.io/pkg/mvtnorm/man/lpmvnorm.html),
  reducing the \\J\\-dimensional integral to a \\(J-1)\\-dimensional
  problem over a hypercube.

- Initial Pivoting: Reorders the integration variables to minimize the
  variance of the integrand.

- Adaptive Sampling: The number of samples per replication increases
  dynamically until the estimated error falls below `abseps` or
  `releps`, or until `n_max` is reached.

The standard error is derived from \\R\\ independent replications. For
high-dimensional problems, computations can be accelerated by setting
`parallel = TRUE`, which distributes the replications across multiple
CPU threads via `nthreads`.

## Author

Kaifeng Lu, <kaifenglu@gmail.com>

## Examples

``` r
# Example 1: Compound symmetry covariance structure and analytic method
n <- 5
mean <- rep(0, n)
lower <- rep(-1, n)
upper <- rep(3, n)
sigma <- matrix(0.5, n, n)
diag(sigma) <- 1
pmvnormr(lower, upper, mean, sigma)
#> [1] 0.5800477
#> attr(,"method")
#> [1] "analytic"
#> attr(,"error")
#> [1] 0
#> attr(,"nsamples")
#> [1] 1

# Example 2: General covariance structure and Monte Carlo method
n <- 5
mean <- rep(0, n)
lower <- rep(-1, n)
upper <- rep(3, n)
sigma <- matrix(c(1, 0.5, 0.3, 0.2, 0.1,
                  0.5, 1, 0.4, 0.3, 0.2,
                  0.3, 0.4, 1, 0.5, 0.3,
                  0.2, 0.3, 0.5, 1, 0.4,
                  0.1, 0.2, 0.3, 0.4, 1), nrow = n)
pmvnormr(lower, upper, mean, sigma)
#> [1] 0.5259119
#> attr(,"method")
#> [1] "qmc"
#> attr(,"error")
#> [1] 7.900235e-05
#> attr(,"nsamples")
#> [1] 8192
```
