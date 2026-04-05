#' @title Probability of Multivariate Normal Distribution over a Hyperrectangle
#' @description Computes the probability that a multivariate normal random vector
#' falls within a specified hyperrectangle defined by lower and upper bounds.
#'
#' @param lower the lower bounds of the hyperrectangle.
#' @param upper the upper bounds of the hyperrectangle.
#' @param mean  the mean vector (optional). If NULL, defaults to rep(0, nrow(sigma)).
#' @param sigma the covariance matrix of the multivariate normal distribution.
#' @param ... forwarded to pmvnormRcpp (fast, n0, n_max, R, abseps, releps, seed)
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
pmvnormr <- function(lower, upper, mean = NULL, sigma, ...) {
  if (!is.matrix(sigma) && length(sigma) == 1) sigma = matrix(sigma, nrow = 1, ncol = 1)
  if (is.null(dim(sigma)) || length(dim(sigma)) != 2L) stop("sigma must be a matrix")
  if (nrow(sigma) != ncol(sigma)) stop("sigma must be square")
  if (is.null(mean)) mean <- rep(0, nrow(sigma))

  out <- pmvnormRcpp(lower = lower, upper = upper, mean = mean, sigma = sigma, ...)
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
#' @param mean  the mean vector (optional). If NULL, defaults to rep(0, nrow(sigma)).
#' @param sigma the covariance matrix of the multivariate normal distribution.
#' @param ... forwarded to qmvnormRcpp (fast, n0, n_max, R, abseps, releps, seed)
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
qmvnormr <- function(p, mean = NULL, sigma, ...) {
  if (!is.matrix(sigma) && length(sigma) == 1) sigma = matrix(sigma, nrow = 1, ncol = 1)
  if (is.null(dim(sigma)) || length(dim(sigma)) != 2L) stop("sigma must be a matrix")
  if (nrow(sigma) != ncol(sigma)) stop("sigma must be square")
  if (is.null(mean)) mean <- rep(0, nrow(sigma))

  qmvnormRcpp(p = p, mean = mean, sigma = sigma, ...)
}
