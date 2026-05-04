#' @title Multivariate Normal Probabilities over a Hyperrectangle
#' @description Computes the probability that a multivariate normal random
#' vector falls within a rectangular region defined by lower and upper bounds.
#'
#' @param lower A numeric vector of lower integration limits.
#' @param upper A numeric vector of upper integration limits.
#' @param mean  The mean vector. If \code{NULL} (default), a zero vector of
#'   appropriate length is used.
#' @param sigma The covariance (or correlation) matrix of the distribution.
#' @param n0 Initial number of samples per replication for the Monte Carlo
#'   integration.
#' @param n_max Maximum number of samples allowed per replication.
#' @param R Number of independent replications used to estimate the error.
#' @param abseps Absolute error tolerance for the probability calculation.
#' @param releps Relative error tolerance for the probability calculation.
#' @param seed Random seed for reproducibility. If 0, a seed is generated
#'   from the computer clock.
#' @param parallel Logical; if \code{TRUE}, computations are performed in parallel.
#' @param nthreads Number of threads for parallel execution. If 0, the
#'   default RcppParallel behavior is used.
#'
#' @details
#' The function automatically selects the most efficient computation method
#' based on the input:
#'
#' * Analytic Methods: Used for univariate cases or multivariate distributions
#'   with a compound symmetry correlation structure and non-negative
#'   correlations.
#'
#' * Monte Carlo Estimation: Used for all other cases. The algorithm employs
#'   a randomized quasi-Monte Carlo (QMC) approach using generalized Halton
#'   sequences.
#'
#' To improve efficiency and accuracy, the QMC approach incorporates:
#'
#' * Sequential Conditioning: Mimics the standardization and transformation
#'   approach used in \code{mvtnorm::lpmvnorm}, reducing the \eqn{J}-dimensional
#'   integral to a \eqn{(J-1)}-dimensional problem over a hypercube.
#'
#' * Adaptive Sampling: The number of samples per replication increases
#'   dynamically until the estimated error falls below \code{abseps} or
#'   \code{releps}, or until \code{n_max} is reached.
#'
#' The standard error is derived from \eqn{R} independent replications.
#' For high-dimensional problems, computations can be accelerated by
#' setting \code{parallel = TRUE}, which distributes the replications across
#' multiple CPU threads via \code{nthreads}.
#'
#' @return The estimated probability with the following attributes:
#'
#' * \code{method}: \code{"exact"} for analytic methods or \code{"qmc"} for
#'   the Monte Carlo approach.
#'
#' * \code{error}: The estimated error (half-width of 95% confidence interval).
#'
#' * \code{nsamples}: The total number of samples used across all replications.
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @examples
#'
#' # Example 1: Compound symmetry covariance structure and analytic method
#' n <- 5
#' mean <- rep(0, n)
#' lower <- rep(-1, n)
#' upper <- rep(3, n)
#' sigma <- matrix(0.5, n, n)
#' diag(sigma) <- 1
#' pmvnormr(lower, upper, mean, sigma)
#'
#' # Example 2: General covariance structure and Monte Carlo method
#' n <- 5
#' mean <- rep(0, n)
#' lower <- rep(-1, n)
#' upper <- rep(3, n)
#' sigma <- matrix(c(1, 0.5, 0.3, 0.2, 0.1,
#'                   0.5, 1, 0.4, 0.3, 0.2,
#'                   0.3, 0.4, 1, 0.5, 0.3,
#'                   0.2, 0.3, 0.5, 1, 0.4,
#'                   0.1, 0.2, 0.3, 0.4, 1), nrow = n)
#' pmvnormr(lower, upper, mean, sigma, seed = 314159)
#'
#' @export
pmvnormr <- function(lower = NULL, upper = NULL, mean = NULL, sigma,
                     n0 = 1024, n_max = 16384, R = 8, abseps = 1e-4,
                     releps = 0.0, seed = 0, parallel = TRUE, nthreads = 0) {
  if (!is.matrix(sigma) && length(sigma) == 1)
    sigma <- matrix(sigma, nrow = 1, ncol = 1)
  if (is.null(dim(sigma)) || length(dim(sigma)) != 2L)
    stop("sigma must be a matrix")
  if (nrow(sigma) != ncol(sigma)) stop("sigma must be square")
  K <- nrow(sigma)
  if (is.null(mean)) mean <- rep(0, K)
  if (is.null(lower)) lower <- rep(-Inf, K)
  if (is.null(upper)) upper <- rep(Inf, K)

  # Respect user-requested number of threads (best effort)
  if (nthreads > 0) {
    n_physical_cores <- parallel::detectCores(logical = FALSE)
    RcppParallel::setThreadOptions(min(nthreads, n_physical_cores))
  }

  out <- pmvnormRcpp(lower = lower, upper = upper, mean = mean, sigma = sigma,
                     n0 = n0, n_max = n_max, R = R, abseps = abseps,
                     releps = releps, seed = seed, parallel = parallel)
  prob <- out$prob
  attr(prob, "method") <- out$method
  attr(prob, "error") <- out$error
  attr(prob, "nsamples") <- out$nsamples
  prob
}


#' @title Equicoordinate Quantiles of the Multivariate Normal Distribution
#' @description Computes the equicoordinate quantile \eqn{q} such that
#' \eqn{P(X_1 \le q, X_2 \le q, \ldots, X_k \le q) = p} for a multivariate
#' normal random vector \eqn{X}.
#'
#' @param p The probability level (cumulative probability).
#' @param mean  The mean vector. If \code{NULL} (default), a zero vector of
#'   appropriate length is used.
#' @param sigma The covariance (or correlation) matrix of the distribution.
#' @param n0 Initial number of samples per replication for the Monte Carlo
#'   integration.
#' @param n_max Maximum number of samples allowed per replication.
#' @param R Number of independent replications used to estimate the error.
#' @param abseps Absolute error tolerance for the probability calculation.
#' @param releps Relative error tolerance for the probability calculation.
#' @param seed Random seed for reproducibility. If 0, a seed is generated
#'   from the computer clock.
#' @param parallel Logical; if \code{TRUE}, computations are performed in parallel.
#' @param nthreads Number of threads for parallel execution. If 0, the
#'   default RcppParallel behavior is used.
#'
#' @details
#' This function finds the value \eqn{q} using a root-finding algorithm
#' applied to the \code{pmvnormr} function. It solves for the value where
#' the multivariate normal cumulative distribution function equals the
#' target probability \eqn{p}.
#'
#' @return A numeric value representing the calculated equicoordinate quantile.
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @examples
#'
#' n <- 5
#' mean <- rep(0, n)
#' sigma <- matrix(0.5, n, n)
#' diag(sigma) <- 1
#' qmvnormr(0.5, mean = mean, sigma = sigma)
#'
#' @export
qmvnormr <- function(p, mean = NULL, sigma,
                     n0 = 1024, n_max = 16384, R = 8, abseps = 1e-4,
                     releps = 0.0, seed = 0, parallel = TRUE, nthreads = 0) {
  if (!is.matrix(sigma) && length(sigma) == 1)
    sigma <- matrix(sigma, nrow = 1, ncol = 1)
  if (is.null(dim(sigma)) || length(dim(sigma)) != 2L)
    stop("sigma must be a matrix")
  if (nrow(sigma) != ncol(sigma)) stop("sigma must be square")
  if (is.null(mean)) mean <- rep(0, nrow(sigma))

  # Respect user-requested number of threads (best effort)
  if (nthreads > 0) {
    n_physical_cores <- parallel::detectCores(logical = FALSE)
    RcppParallel::setThreadOptions(min(nthreads, n_physical_cores))
  }

  qmvnormRcpp(p = p, mean = mean, sigma = sigma,
              n0 = n0, n_max = n_max, R = R, abseps = abseps,
              releps = releps, seed = seed, parallel = parallel)
}
