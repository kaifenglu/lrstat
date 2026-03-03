


#' @title Density Function of Truncated Piecewise Exponential
#' Distribution
#' @description Obtains the density of a truncated piecewise exponential
#' distribution.
#'
#' @param q The vector of quantiles.
#' @inheritParams param_piecewiseSurvivalTime
#' @inheritParams param_lambda
#' @param lowerBound The left truncation time point for the survival time.
#'   Defaults to 0 for no truncation.
#' @param log.d Logical; if TRUE, densities d are given as log(d).

#' @return The density d such that
#' d = lambda(q) * P(X > q | X > lowerBound).
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @examples
#' dtpwexp(q = c(8, 18), piecewiseSurvivalTime = c(0, 6, 9, 15),
#'         lambda = c(0.025, 0.04, 0.015, 0.007))
#'
#' @export
dtpwexp <- function(q, piecewiseSurvivalTime = 0, lambda = 0.0578,
                    lowerBound = 0, log.d = FALSE) {

  if (any(q < 0)) {
    stop("q must be nonnegative")
  }

  if (piecewiseSurvivalTime[1] != 0) {
    stop("piecewiseSurvivalTime must start with 0")
  }

  if (length(piecewiseSurvivalTime) > 1 &&
      any(diff(piecewiseSurvivalTime) <= 0)) {
    stop("piecewiseSurvivalTime should be increasing")
  }

  if (length(lambda) != length(piecewiseSurvivalTime)) {
    stop("lambda and piecewiseSurvivalTime must have the same length")
  }

  if (any(lambda < 0)) {
    stop("lambda must be nonnegative")
  }

  if (lowerBound < 0) {
    stop("lowerBound must be nonnegative")
  }

  dtpwexpcpp(q = q, piecewiseSurvivalTime = piecewiseSurvivalTime,
             lambda = lambda, lowerBound = lowerBound, logd = log.d)
}


#' @title Distribution Function of Truncated Piecewise Exponential
#' Distribution
#' @description Obtains the probability of a truncated piecewise exponential
#' distribution.
#'
#' @param q The vector of quantiles.
#' @inheritParams param_piecewiseSurvivalTime
#' @inheritParams param_lambda
#' @param lowerBound The left truncation time point for the survival time.
#'   Defaults to 0 for no truncation.
#' @param lower.tail Logical; if TRUE (default), probabilities are P(X <= x),
#'   otherwise, P(X > x).
#' @param log.p Logical; if TRUE, probabilities p are given as log(p).

#' @return The probability p such that
#' P(X > q | X > lowerBound) = 1 - p.
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @examples
#' ptpwexp(q = c(8, 18), piecewiseSurvivalTime = c(0, 6, 9, 15),
#'         lambda = c(0.025, 0.04, 0.015, 0.007))
#'
#' @export
ptpwexp <- function(q, piecewiseSurvivalTime = 0, lambda = 0.0578,
                    lowerBound = 0, lower.tail = TRUE, log.p = FALSE) {

  if (any(q < 0)) {
    stop("q must be nonnegative")
  }

  if (piecewiseSurvivalTime[1] != 0) {
    stop("piecewiseSurvivalTime must start with 0")
  }

  if (length(piecewiseSurvivalTime) > 1 &&
      any(diff(piecewiseSurvivalTime) <= 0)) {
    stop("piecewiseSurvivalTime should be increasing")
  }

  if (length(lambda) != length(piecewiseSurvivalTime)) {
    stop("lambda and piecewiseSurvivalTime must have the same length")
  }

  if (any(lambda < 0)) {
    stop("lambda must be nonnegative")
  }

  if (lowerBound < 0) {
    stop("lowerBound must be nonnegative")
  }

  ptpwexpcpp(q = q, piecewiseSurvivalTime = piecewiseSurvivalTime,
             lambda = lambda, lowerBound = lowerBound,
             lowertail = lower.tail, logp = log.p)
}


#' @title Quantile Function of Truncated Piecewise Exponential Distribution
#' @description Obtains the quantile of a truncated piecewise exponential
#' distribution.
#'
#' @param p The vector of probabilities.
#' @inheritParams param_piecewiseSurvivalTime
#' @inheritParams param_lambda
#' @param lowerBound The left truncation time point for the survival time.
#'   Defaults to 0 for no truncation.
#' @param lower.tail Logical; if TRUE (default), probabilities are P(X <= x),
#'   otherwise, P(X > x).
#' @param log.p Logical; if TRUE, probabilities p are given as log(p).
#'
#' @return The quantile q such that
#' P(X > q | X > lowerBound) = 1 - p.
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @examples
#' qtpwexp(p = c(0.205, 0.317), piecewiseSurvivalTime = c(0, 6, 9, 15),
#'         lambda = c(0.025, 0.04, 0.015, 0.007))
#'
#' @export
qtpwexp <- function(p, piecewiseSurvivalTime = 0, lambda = 0.0578,
                    lowerBound = 0, lower.tail = TRUE, log.p = FALSE) {
  if (any(p < 0 | p > 1)) {
    stop("p must lie between 0 and 1")
  }

  if (piecewiseSurvivalTime[1] != 0) {
    stop("piecewiseSurvivalTime must start with 0")
  }

  if (length(piecewiseSurvivalTime) > 1 &&
      any(diff(piecewiseSurvivalTime) <= 0)) {
    stop("piecewiseSurvivalTime should be increasing")
  }

  if (length(lambda) != length(piecewiseSurvivalTime)) {
    stop("lambda and piecewiseSurvivalTime must have the same length")
  }

  if (any(lambda < 0)) {
    stop("lambda must be nonnegative")
  }

  if (lowerBound < 0) {
    stop("lowerBound must be nonnegative")
  }

  qtpwexpcpp(p = p, piecewiseSurvivalTime = piecewiseSurvivalTime,
             lambda = lambda, lowerBound = lowerBound,
             lowertail = lower.tail, logp = log.p)
}



#' @title Random Number Generation Function of Truncated Piecewise
#' Exponential Distribution
#' @description Obtains random samples from a truncated piecewise
#' exponential distribution.
#'
#' @param n The number of observations.
#' @inheritParams param_piecewiseSurvivalTime
#' @inheritParams param_lambda
#' @param lowerBound The left truncation time point for the survival time.
#'   Defaults to 0 for no truncation.
#'
#' @return The random numbers from truncated piecewise exponential
#' distribution.
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @examples
#' rtpwexp(n = 10, piecewiseSurvivalTime = c(0, 6, 9, 15),
#'         lambda = c(0.025, 0.04, 0.015, 0.007))
#'
#' @export
rtpwexp <- function(n, piecewiseSurvivalTime = 0, lambda = 0.0578,
                    lowerBound = 0) {
  if (n <= 0 || n != round(n)) {
    stop("n must be a positive integer")
  }

  if (piecewiseSurvivalTime[1] != 0) {
    stop("piecewiseSurvivalTime must start with 0")
  }

  if (length(piecewiseSurvivalTime) > 1 &&
      any(diff(piecewiseSurvivalTime) <= 0)) {
    stop("piecewiseSurvivalTime should be increasing")
  }

  if (length(lambda) != length(piecewiseSurvivalTime)) {
    stop("lambda and piecewiseSurvivalTime must have the same length")
  }

  if (any(lambda < 0)) {
    stop("lambda must be nonnegative")
  }

  if (lowerBound < 0) {
    stop("lowerBound must be nonnegative")
  }

  qtpwexpcpp(p = runif(n), piecewiseSurvivalTime = piecewiseSurvivalTime,
             lambda = lambda, lowerBound = lowerBound,
             lowertail = TRUE, logp = FALSE)
}

