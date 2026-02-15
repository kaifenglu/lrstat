#' @title Error Spending
#' @description Obtains the error spent at given spending times
#' for the specified error spending function.
#'
#' @param t A vector of spending times, typically equal to information
#'   fractions.
#' @param error The total error to spend.
#' @param sf The spending function. One of the following: "sfOF" for
#'   O'Brien-Fleming type spending function, "sfP" for Pocock type spending
#'   function, "sfKD" for Kim & DeMets spending function, and "sfHSD" for
#'   Hwang, Shi & DeCani spending function. Defaults to "sfOF".
#' @param sfpar The parameter for the spending function. Corresponds to
#'   rho for "sfKD" and gamma for "sfHSD".
#'
#' @details
#' This function implements a variety of error spending functions commonly
#' used in group sequential designs, assuming one-sided hypothesis testing.
#'
#' **O'Brien-Fleming-Type Spending Function**
#'
#' This spending function allocates very little alpha early on and more alpha
#' later in the trial. It is defined as:
#' \deqn{
#' \alpha(t) = 2 - 2\Phi\left(\frac{z_{\alpha/2}}{\sqrt{t}}\right),
#' }
#' where \eqn{\Phi} is the standard normal cumulative distribution function,
#' \eqn{z_{\alpha/2}} is the critical value from the standard normal
#' distribution, and \eqn{t \in [0, 1]} denotes the information fraction.
#'
#' **Pocock-Type Spending Function**
#'
#' This function spends alpha more evenly throughout the study:
#' \deqn{
#' \alpha(t) = \alpha \log(1 + (e - 1)t),
#' }
#' where \eqn{e} is Euler's number (approximately 2.718).
#'
#' **Kim and DeMets Power-Type Spending Function**
#'
#' This family of spending functions is defined as:
#' \deqn{
#' \alpha(t) = \alpha t^{\rho}, \quad \rho > 0.
#' }
#' - When \eqn{\rho = 1}, the function mimics Pocock-type boundaries.
#' - When \eqn{\rho = 3}, it approximates O’Brien-Fleming-type boundaries.
#'
#' **Hwang, Shih, and DeCani Spending Function**
#'
#' This flexible family of functions is given by:
#' \deqn{
#' \alpha(t) =
#' \begin{cases}
#' \alpha \frac{1 - e^{-\gamma t}}{1 - e^{-\gamma}}, & \text{if }
#' \gamma \ne 0 \\ \alpha t, & \text{if } \gamma = 0.
#' \end{cases}
#' }
#' - When \eqn{\gamma = -4}, the spending function resembles
#'   O’Brien-Fleming boundaries.
#' - When \eqn{\gamma = 1}, it resembles Pocock boundaries.
#'
#' @return A vector of errors spent up to the interim look.
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @examples
#'
#' errorSpent(t = 0.5, error = 0.025, sf = "sfOF")
#'
#' errorSpent(t = c(0.5, 0.75, 1), error = 0.025, sf = "sfHSD", sfpar = -4)
#'
#' @export
errorSpent <- function(t, error = 0.025, sf = "sfOF", sfpar = NA) {
  sapply(t, errorSpentcpp, error = error, sf = sf, sfpar = sfpar)
}


#' @title Stagewise Exit Probabilities
#' @description Obtains the stagewise exit probabilities for both efficacy
#' and futility stopping.
#'
#' @param b Upper boundaries on the z-test statistic scale.
#' @param a Lower boundaries on the z-test statistic scale. Defaults to
#'   \code{c(rep(-6.0, kMax-1), b[kMax])} if left unspecified, where
#'   \code{kMax = length(b)}.
#' @param theta Stagewise parameter of interest, e.g., \code{-U/V} for
#'   weighted log-rank test, where \code{U} is the mean and \code{V} is
#'   the variance of the weighted log-rank test score statistic at each
#'   stage. For proportional hazards and conventional log-rank test, use the
#'   scalar input, \code{theta = -log(HR)}. Defaults to 0 corresponding to
#'   the null hypothesis.
#' @param I Stagewise cumulative information, e.g., \code{V}, the variance
#'   of the weighted log-rank test score statistic at each stage. For
#'   conventional log-rank test, information can be approximated by
#'   \code{phi*(1-phi)*D}, where \code{phi} is the probability of being
#'   allocated to the active arm, and \code{D} is the total number of events
#'   at each stage. Defaults to \code{seq(1, kMax)} if left unspecified.
#'
#' @return A list of stagewise exit probabilities:
#'
#' * \code{exitProbUpper}: The vector of efficacy stopping probabilities
#'
#' * \code{exitProbLower}: The vector of futility stopping probabilities.
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @examples
#' exitprob(b = c(3.471, 2.454, 2.004), theta = -log(0.6),
#'          I = c(50, 100, 150)/4)
#'
#' exitprob(b = c(2.963, 2.359, 2.014),
#'          a = c(-0.264, 0.599, 2.014),
#'          theta = c(0.141, 0.204, 0.289),
#'          I = c(81, 121, 160))
#'
#' @export
exitprob <- function(b, a = NA, theta = 0, I = NA) {
  exitprobcpp(b = b, a = a, theta = theta, I = I)
}


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


#' @title Mean and Variance of Truncated Piecewise Exponential Distribution
#' @description Obtains the mean and variance from a truncated piecewise
#' exponential distribution.
#'
#' @inheritParams param_piecewiseSurvivalTime
#' @inheritParams param_lambda
#' @param lowerBound The left truncation time point for the survival time.
#'   Defaults to 0 for no truncation.
#'
#' @return A list with two components, one for the mean, and the other for
#' the variance of the truncated piecewise exponential distribution.
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @examples
#' mtpwexp(piecewiseSurvivalTime = c(0, 6, 9, 15),
#'         lambda = c(0.025, 0.04, 0.015, 0.007))
#'
#' @export
mtpwexp <- function(piecewiseSurvivalTime = 0, lambda = 0.0578,
                    lowerBound = 0) {
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

  mtpwexpcpp(piecewiseSurvivalTime = piecewiseSurvivalTime,
             lambda = lambda, lowerBound = lowerBound)
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
             lambda = lambda, lowerBound = lowerBound)
}


#' @title Distribution Function of the Standard Bivariate Normal
#' @description Computes the cumulative distribution function (CDF) of
#' the standard bivariate normal distribution with specified lower and
#' upper integration limits and correlation coefficient.
#'
#' @param lower A numeric vector of length 2 specifying the lower limits
#'   of integration.
#' @param upper A numeric vector of length 2 specifying the upper limits
#'   of integration.
#' @param rho A numeric value specifying the correlation coefficient of
#'   the standard bivariate normal distribution.
#'
#' @details This function evaluates the probability
#' \eqn{P(\code{lower[1]} < X < \code{upper[1]},
#' \code{lower[2]} < Y < \code{upper[2]})} where
#' \eqn{(X, Y)} follows a standard bivariate normal
#' distribution with correlation \code{corr}.
#'
#' @return A numeric value representing the probability that a standard
#' bivariate normal vector falls within the specified rectangular region.
#'
#' @examples
#' pbvnorm(c(-1, -1), c(1, 1), 0.5)
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @export
pbvnorm <- function(lower = c(-Inf, Inf),
                    upper = c(Inf, Inf), rho = 0) {

  if (length(lower) != 1 && length(lower) != 2) {
    stop("lower must be a numerical vector of length 1 or 2")
  } else if (length(lower) == 1) {
    lower = rep(lower, 2)
  }

  if (length(upper) != 1 && length(upper) != 2) {
    stop("upper must be a numerical vector of length 1 or 2")
  } else if (length(upper) == 1) {
    upper = rep(upper, 2)
  }

  if (any(lower >= upper)) {
    stop("lower must be less than upper")
  }

  if (rho >= 1 || rho <= -1) {
    stop("corr must lie between -1 and 1")
  }

  pbvnormcpp(lower, upper, rho)
}


#' @title Hazard Function for Progressive Disease (PD) Given Correlation
#' Between PD and OS
#'
#' @description
#' Computes the hazard function of a piecewise exponential
#' distribution for progressive disease (PD), such that the
#' resulting hazard function for progression-free survival (PFS)
#' closely matches a given piecewise hazard for PFS.
#'
#' @inheritParams param_piecewiseSurvivalTime
#' @param hazard_pfs A scalar or numeric vector specifying the
#'   hazard(s) for PFS based on a piecewise exponential distribution.
#' @param hazard_os A scalar or numeric vector specifying the
#'   hazard(s) for overall survival (OS) based on a piecewise
#'   exponential distribution.
#' @param rho_pd_os A numeric value specifying the correlation
#'   between PD and OS times.
#'
#' @details
#' This function determines the hazard vector \eqn{\lambda_{\text{pd}}}
#' for the piecewise exponential distribution of PD, so that the
#' implied survival function for PFS time,
#' \eqn{T_{\text{pfs}} = \min(T_{\text{pd}}, T_{\text{os}})}, closely
#' matches the specified piecewise exponential distribution for PFS
#' with hazard vector \eqn{\lambda_{\text{pfs}}}.
#'
#' To achieve this, we simulate
#' \eqn{(Z_{\text{pd}}, Z_{\text{os}})} from
#' a standard bivariate normal distribution with correlation
#' \eqn{\rho}. Then, \eqn{U_{\text{pd}} = \Phi(Z_{\text{pd}})}
#' and \eqn{U_{\text{os}} = \Phi(Z_{\text{os}})} are generated, where
#' \eqn{\Phi} denotes the standard normal CDF.
#'
#' The times to PD and OS are obtained via the inverse transform
#' method using quantile functions of the piecewise exponential distribution:
#' \deqn{T_{\text{pd}} = \text{qpwexp}(U_{\text{pd}},u,\lambda_{\text{pd}})}
#' \deqn{T_{\text{os}} = \text{qpwexp}(U_{\text{os}},u,\lambda_{\text{os}})}
#' where \code{u = piecewiseSurvivalTime}.
#'
#' The function solves for \eqn{\lambda_{\text{pd}}} such that
#' the survival function of \eqn{T_{\text{pfs}}} closely matches that
#' of a piecewise exponential distribution with hazard
#' \eqn{\lambda_{\text{pfs}}}:
#' \deqn{P(\min(T_{\text{pd}}, T_{\text{os}}) > t) = S_{\text{pfs}}(t)}
#' Since \deqn{Z_{\text{pd}} =
#'   \Phi^{-1}(\text{ppwexp}(T_\text{pd}, u, \lambda_{\text{pd}}))} and
#' \deqn{Z_{\text{os}} =
#'   \Phi^{-1}(\text{ppwexp}(T_\text{os}, u, \lambda_{\text{os}}))}
#' we have
#' \deqn{P(\min(T_{\text{pd}}, T_{\text{os}}) > t) =
#' P(Z_{\text{pd}} >
#'     \Phi^{-1}(\text{ppwexp}(t,u,\lambda_{\text{pd}})),
#'   Z_{\text{os}} >
#'     \Phi^{-1}(\text{ppwexp}(t,u,\lambda_{\text{os}})))}
#' while
#' \deqn{S_{\text{pfs}}(t) = 1 - \text{ppwexp}(t,u,\lambda_{\text{pfs}})}
#'
#' Matching is performed sequentially at the internal cut points
#' \eqn{u_2, ..., u_J} and at the point
#' \eqn{u_J + \log(2)/\lambda_{\text{pfs},J}} for the final interval,
#' as well as the percentile points at 10%, 20%, ..., 90%, and 95%
#' to solve for \eqn{\lambda_{\text{pd},1}, \ldots,
#' \lambda_{\text{pd},K}}, where \eqn{K} is the total number of
#' unique cut points.
#'
#' @return A list with the following components:
#'
#' * \code{piecewiseSurvivalTime}: A vector that specifies the starting time
#'   points of the intervals for the piecewise exponential distribution
#'   for PD.
#'
#' * \code{hazard_pd}: A numeric vector representing the calculated hazard
#'   rates for the piecewise exponential distribution of PD.
#'
#' * \code{hazard_os}: A numeric vector representing the hazard rates for
#'   the piecewise exponential distribution of OS at the same time points
#'   as PD.
#'
#' * \code{rho_pd_os}: The correlation between PD and OS times (as input).
#'
#' @author
#' Kaifeng Lu (\email{kaifenglu@gmail.com})
#'
#' @examples
#' u <- c(0, 1, 3, 4)
#' lambda1 <- c(0.0151, 0.0403, 0.0501, 0.0558)
#' lambda2 <- 0.0145
#' rho_pd_os <- 0.5
#' hazard_pd(u, lambda1, lambda2, rho_pd_os)
#'
#' @export
hazard_pd <- function(piecewiseSurvivalTime = 0,
                      hazard_pfs = 0.0578,
                      hazard_os = 0.02,
                      rho_pd_os = 0.5) {

  if (piecewiseSurvivalTime[1] != 0) {
    stop("piecewiseSurvivalTime must start with 0")
  }

  if (length(piecewiseSurvivalTime) > 1 &&
      any(diff(piecewiseSurvivalTime) <= 0)) {
    stop("piecewiseSurvivalTime should be increasing")
  }

  k = length(piecewiseSurvivalTime)
  if (!(length(hazard_pfs) %in% c(1, k))) {
    stop(paste("hazard_pfs must be a scalar or have the same length",
               "as piecewiseSurvivalTime"))
  } else if (length(hazard_pfs) == 1) {
    hazard_pfs = rep(hazard_pfs, k)
  }

  if (!(length(hazard_os) %in% c(1, k))) {
    stop(paste("hazard_os must be a scalar or have the same length",
               "as piecewiseSurvivalTime"))
  } else if (length(hazard_os) == 1) {
    hazard_os = rep(hazard_os, k)
  }

  if (any(hazard_os <= 0)) {
    stop("hazard_os must be positive")
  }

  if (any(hazard_pfs <= hazard_os)) {
    stop("hazard_pfs must be greater than hazard_os")
  }

  if (rho_pd_os <= -1 || rho_pd_os >= 1) {
    stop("corr_pd_os must lie between -1 and 1")
  }

  hazard_pdcpp(piecewiseSurvivalTime, hazard_pfs, hazard_os, rho_pd_os)
}


#' @title Correlation Between PFS and OS Given Correlation Between PD and OS
#'
#' @description
#' Computes the correlation between PFS and OS given the correlation
#' between PD and OS.
#'
#' @inheritParams param_piecewiseSurvivalTime
#' @param hazard_pfs A scalar or numeric vector specifying the
#'   hazard(s) for PFS based on a piecewise exponential distribution.
#' @param hazard_os A scalar or numeric vector specifying the
#'   hazard(s) for overall survival (OS) based on a piecewise
#'   exponential distribution.
#' @param rho_pd_os A numeric value specifying the correlation
#'   between PD and OS times.
#'
#' @details
#' This function first determines the piecewise exponential distribution
#' for PD such that the implied survival function for PFS time,
#' \eqn{T_{\text{pfs}} = \min(T_{\text{pd}}, T_{\text{os}})}, closely
#' matches the specified piecewise exponential distribution for PFS
#' with hazard vector \eqn{\lambda_{\text{pfs}}}. Then, it calculates
#' the correlation between PFS and OS times based on the derived
#' piecewise exponential distribution for PD and the given piecewise
#' exponential distribution for OS.
#'
#' @return The estimated correlation between PFS and OS.
#'
#' @author
#' Kaifeng Lu (\email{kaifenglu@gmail.com})
#'
#' @examples
#' u <- c(0, 1, 3, 4)
#' lambda1 <- c(0.0151, 0.0403, 0.0501, 0.0558)
#' lambda2 <- 0.0145
#' rho_pd_os <- 0.5
#' corr_pfs_os(u, lambda1, lambda2, rho_pd_os)
#'
#' @export
corr_pfs_os <- function(piecewiseSurvivalTime, hazard_pfs,
                        hazard_os, rho_pd_os) {

  if (piecewiseSurvivalTime[1] != 0) {
    stop("piecewiseSurvivalTime must start with 0")
  }

  if (length(piecewiseSurvivalTime) > 1 &&
      any(diff(piecewiseSurvivalTime) <= 0)) {
    stop("piecewiseSurvivalTime should be increasing")
  }

  k = length(piecewiseSurvivalTime)
  if (!(length(hazard_pfs) %in% c(1, k))) {
    stop(paste("hazard_pfs must be a scalar or have the same length",
               "as piecewiseSurvivalTime"))
  } else if (length(hazard_pfs) == 1) {
    hazard_pfs = rep(hazard_pfs, k)
  }

  if (!(length(hazard_os) %in% c(1, k))) {
    stop(paste("hazard_os must be a scalar or have the same length",
               "as piecewiseSurvivalTime"))
  } else if (length(hazard_os) == 1) {
    hazard_os = rep(hazard_os, k)
  }

  if (any(hazard_os <= 0)) {
    stop("hazard_os must be positive")
  }

  if (any(hazard_pfs <= hazard_os)) {
    stop("hazard_pfs must be greater than hazard_os")
  }

  if (rho_pd_os <= -1 || rho_pd_os >= 1) {
    stop("corr_pd_os must lie between -1 and 1")
  }

  corr_pfs_oscpp(piecewiseSurvivalTime, hazard_pfs, hazard_os, rho_pd_os);
}


#' @title Hazard Function for Sub Population
#'
#' @description
#' Computes the hazard function of a piecewise exponential
#' distribution for the biomarker negative sub population, such that the
#' resulting survival function for the ITT population
#' closely matches a given piecewise survival function.
#'
#' @inheritParams param_piecewiseSurvivalTime
#' @param hazard_itt A scalar or numeric vector specifying the
#'   hazard(s) for the ITT population based on a piecewise exponential
#'   distribution.
#' @param hazard_pos A scalar or numeric vector specifying the
#'   hazard(s) for the biomarker positive sub population
#'   based on a piecewise exponential distribution.
#' @param p_pos A numeric value specifying the prevalence of the
#'   biomarker positive sub population.
#'
#' @details
#' This function determines the hazard vector \eqn{\lambda_{\text{neg}}}
#' for the piecewise exponential distribution of the biomarker negative
#' sub population, so that the implied survival function for the ITT
#' population closely matches the specified piecewise exponential
#' distribution with hazard vector \eqn{\lambda_{\text{itt}}}.
#'
#' Let \eqn{p_{\text{pos}}} be the
#' prevalence of the biomarker positive sub population,
#' then the survival function for the ITT population is given by
#' \deqn{S_{\text{itt}}(t) = p_{\text{pos}} S_{\text{pos}}(t) +
#' (1 - p_{\text{pos}}) S_{\text{neg}}(t)}
#' where \eqn{S_{\text{pos}}(t)} and \eqn{S_{\text{neg}}(t)} are
#' the survival functions for the biomarker positive and
#' biomarker negative sub populations, respectively.
#'
#' Matching is performed sequentially at the internal cutpoints
#' \eqn{u_2, ..., u_J} and at the point
#' \eqn{u_J + \log(2)/\lambda_{\text{itt},J}} for the final interval,
#' as well as the percentile points at 10%, 20%, ..., 90%, and 95%,
#' to solve for \eqn{\lambda_{\text{neg},1}, \ldots,
#' \lambda_{\text{neg},K}}, where \eqn{K} is the total number of
#' unique cut points.
#'
#' @return A list with the following components:
#'
#' * \code{piecewiseSurvivalTime}: A vector that specifies the starting time
#'   points of the intervals for the piecewise exponential distribution
#'   for the biomarker negative sub population.
#'
#' * \code{hazard_pos}: A numeric vector representing the hazard rates for
#'   the piecewise exponential distribution of the biomarker positive
#'   sub population at the same time points as the biomarker negative
#'   sub population.
#'
#' * \code{hazard_neg}: A numeric vector representing the estimated hazard
#'   rates for the piecewise exponential distribution of the biomarker
#'   negative sub population.
#'
#' * \code{p_pos}: The prevalence of the biomarker positive sub population
#'  (as input).
#'
#' @author
#' Kaifeng Lu (\email{kaifenglu@gmail.com})
#'
#' @examples
#' u <- c(0, 1, 3, 4)
#' lambda_itt <- c(0.0151, 0.0403, 0.0501, 0.0558)
#' lambda_pos <- c(0.0115, 0.0302, 0.0351, 0.0404)
#' p_pos <- 0.3
#' hazard_sub(u, lambda_itt, lambda_pos, p_pos)
#'
#' @export
hazard_sub <- function(piecewiseSurvivalTime = 0,
                       hazard_itt = 0.0578,
                       hazard_pos = 0.02,
                       p_pos = 0.5) {

  if (piecewiseSurvivalTime[1] != 0) {
    stop("piecewiseSurvivalTime must start with 0")
  }

  if (length(piecewiseSurvivalTime) > 1 &&
      any(diff(piecewiseSurvivalTime) <= 0)) {
    stop("piecewiseSurvivalTime should be increasing")
  }

  k = length(piecewiseSurvivalTime)
  if (!(length(hazard_itt) %in% c(1, k))) {
    stop(paste("hazard_itt must be a scalar or have the same length",
               "as piecewiseSurvivalTime"))
  } else if (length(hazard_itt) == 1) {
    hazard_itt = rep(hazard_itt, k)
  }

  if (!(length(hazard_pos) %in% c(1, k))) {
    stop(paste("hazard_pos must be a scalar or have the same length",
               "as piecewiseSurvivalTime"))
  } else if (length(hazard_pos) == 1) {
    hazard_pos = rep(hazard_pos, k)
  }

  if (any(hazard_pos <= 0)) {
    stop("hazard_pos must be positive")
  }

  if (any(hazard_itt <= hazard_pos)) {
    stop("hazard_itt must be greater than hazard_pos")
  }

  if (p_pos <= 0 || p_pos >= 1) {
    stop("p_pos must lie between 0 and 1")
  }

  hazard_subcpp(piecewiseSurvivalTime, hazard_itt, hazard_pos, p_pos)
}

