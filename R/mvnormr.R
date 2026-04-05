#' @title Probability of Multivariate Normal Distribution over a Hyperrectangle
#' @description Computes the probability that a multivariate normal random vector
#' falls within a specified hyperrectangle defined by lower and upper bounds.
#'
#' @param lower the lower bounds of the hyperrectangle.
#' @param upper the upper bounds of the hyperrectangle.
#' @param mean  the mean vector (optional). If NULL, defaults to
#'   \code{rep(0, nrow(sigma))}.
#' @param sigma the covariance matrix of the multivariate normal distribution.
#' @param fast whether to use the fast version of pnorm and qnorm.
#' @param n0 the initial number of samples for the Monte Carlo estimation
#'   per replication.
#' @param n_max the maximum number of samples for the Monte Carlo estimation
#'   per replication.
#' @param R the number of replications for the Monte Carlo estimation.
#' @param abseps the absolute error tolerance for the Monte Carlo estimation.
#' @param releps the relative error tolerance for the Monte Carlo estimation.
#' @param seed the random seed for reproducibility. If 0 (default), seed
#'   will be generated randomly.
#' @param parallel whether to use parallel computation for the Monte Carlo
#'   estimation.
#' @param nthreads the number of threads to use in simulations (0 means
#'   the default RcppParallel behavior).
#'
#' @return list(prob=..., error=..., n=...)
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @examples
#' n <- 5
#' mean <- rep(0, 5)
#' lower <- rep(-1, 5)
#' upper <- rep(3, 5)
#' corr <- diag(5)
#' corr[lower.tri(corr)] <- 0.5
#' corr[upper.tri(corr)] <- 0.5
#' prob <- pmvnormr(lower, upper, mean, corr)
#' print(prob)
#'
#' @export
pmvnormr <- function(lower, upper, mean = NULL, sigma, fast = TRUE,
                     n0 = 128, n_max = 16384, R = 8, abseps = 1e-4,
                     releps = 0.0, seed = 0, parallel = TRUE, nthreads = 0) {
  if (!is.matrix(sigma) && length(sigma) == 1)
    sigma = matrix(sigma, nrow = 1, ncol = 1)
  if (is.null(dim(sigma)) || length(dim(sigma)) != 2L)
    stop("sigma must be a matrix")
  if (nrow(sigma) != ncol(sigma)) stop("sigma must be square")
  if (is.null(mean)) mean <- rep(0, nrow(sigma))

  # Respect user-requested number of threads (best effort)
  if (nthreads > 0) {
    n_physical_cores <- parallel::detectCores(logical = FALSE)
    RcppParallel::setThreadOptions(min(nthreads, n_physical_cores))
  }

  out <- pmvnormRcpp(lower = lower, upper = upper, mean = mean, sigma = sigma,
                     fast = fast, n0 = n0, n_max = n_max, R = R, abseps = abseps,
                     releps = releps, seed = seed, parallel = parallel)
  prob <- out$prob
  attr(prob, "error") <- out$error
  attr(prob, "n") <- out$n
  prob
}


#' @title Quantile of the Multivariate Normal Distribution
#' @description Computes the equicoordinate quantile function of the
#' multivariate normal distribution.
#'
#' @param p the probability.
#' @param mean  the mean vector (optional). If NULL, defaults to
#'   \code{rep(0, nrow(sigma))}.
#' @param sigma the covariance matrix of the multivariate normal distribution.
#' @param fast whether to use the fast version of pnorm and qnorm.
#' @param n0 the initial number of samples for the Monte Carlo estimation
#'   per replication.
#' @param n_max the maximum number of samples for the Monte Carlo estimation
#'   per replication.
#' @param R the number of replications for the Monte Carlo estimation.
#' @param abseps the absolute error tolerance for the Monte Carlo estimation.
#' @param releps the relative error tolerance for the Monte Carlo estimation.
#' @param seed the random seed for reproducibility. If 0 (default), seed
#'   will be generated randomly.
#' @param parallel whether to use parallel computation for the Monte Carlo
#'   estimation.
#' @param nthreads the number of threads to use in simulations (0 means
#'   the default RcppParallel behavior).
#'
#' @return The calculated quantile.
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @examples
#' n <- 5
#' mean <- rep(0, 5)
#' corr <- diag(5)
#' corr[lower.tri(corr)] <- 0.5
#' corr[upper.tri(corr)] <- 0.5
#' qmvnormr(0.5, mean, corr)
#'
#' @export
qmvnormr <- function(p, mean = NULL, sigma, fast = TRUE,
                     n0 = 128, n_max = 16384, R = 8, abseps = 1e-4,
                     releps = 0.0, seed = 0, parallel = TRUE, nthreads = 0) {
  if (!is.matrix(sigma) && length(sigma) == 1)
    sigma = matrix(sigma, nrow = 1, ncol = 1)
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
              fast = fast, n0 = n0, n_max = n_max, R = R, abseps = abseps,
              releps = releps, seed = seed, parallel = parallel)
}
