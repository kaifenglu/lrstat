#' @name lrstat-package
#' @aliases lrstat-package
#' @docType package
#'
#' @title Power and Sample Size Calculation for Non-Proportional Hazards
#'
#' @description Performs power and sample size calculation for non-proportional
#' hazards model using the Fleming-Harrington family of weighted log-rank tests.
#'
#' @details For proportional hazards, the power is determined by the
#' total number of events and the constant hazard ratio along with
#' information rates and spending functions. For
#' non-proportional hazards, the hazard ratio varies over time and
#' the calendar time plays a key role in determining the mean and
#' variance of the log-rank test score statistic. It requires an
#' iterative algorithm to find the calendar time at which the
#' targeted number of events will be reached for each interim
#' analysis. The lrstat package uses the analytic method 
#' in Lu (2021) to find the mean and variance of the weighted
#' log-rank test score statistic at each interim analysis. In
#' addition, the package approximates the variance and covariance
#' matrix of the sequentially calculated log-rank test statistics
#' under the alternative hypothesis with that under the null hypothesis
#' to take advantage of the independent increments structure in
#' Tsiatis (1982) applicable for the Fleming-Harrington family
#' of weighted log-rank tests.
#'
#' The most useful functions in the package are \code{lrstat}, \code{lrpower},
#' \code{lrsamplesize}, and \code{lrsim}, which calculate the mean and variance
#' of log-rank test score statistic at a sequence of given calendar
#' times, the power of the log-rank test, the sample size in terms
#' of accrual duration and follow-up duration, and the log-rank
#' test simulation, respectively. The \code{accrual} function calculates
#' the number of patients accrued at given calendar times. The
#' \code{caltime} function finds the calendar times to reach the targeted
#' number of events. The \code{exitprob} function calculates the stagewise
#' exit probabilities for specified boundaries with a varying mean
#' parameter over time based on an adaptation of the recursive
#' integration algorithm described in Chapter 19 of Jennison and
#' Turnbull (2000).
#'
#' The development of the lrstat package is heavily influenced by
#' the rpact package. We find their function arguments to be
#' self-explanatory. We have used the same names whenever appropriate
#' so that users familiar with the rpact package can learn the
#' lrstat package quickly. However, there are notable differences:
#' \itemize{
#'   \item lrstat uses direct approximation, while rpact uses the Schoenfeld
#'     method for log-rank test power and sample size calculation.
#'   \item lrstat uses \code{accrualDuration} to explicitly set the end of
#'     accrual period, while rpact incorporates the end of accrual period
#'     in \code{accrualTime}.
#'   \item lrstat considers the trial a failure at the last stage if
#'     the log-rank test cannot reject the null hypothesis up to this
#'     stage and cannot stop for futility at an earlier stage.
#'   \item the \code{lrsim} function uses the variance of the log-rank
#'     test score statistic as the information.
#' }
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @references
#' Anastasios A. Tsiatis. Repeated significance testing for a general class
#' of statistics used in censored survival analysis. J Am Stat Assoc.
#' 1982;77:855-861. 
#'
#' Christopher Jennison, Bruce W. Turnbull. Group Sequential Methods with
#' Applications to Clinical Trials. Chapman & Hall/CRC: Boca Raton, 2000,
#' ISBN:0849303168
#'
#' Kaifeng Lu. Sample size calculation for logrank test and prediction of
#' number of events over time. Pharm Stat. 2021;20:229-244.
#' 
#'
#' @seealso rpact, gsDesign
#'
#' @examples
#' lrpower(kMax = 2, informationRates = c(0.8, 1),
#'         criticalValues = c(2.250, 2.025), accrualIntensity = 20,
#'         piecewiseSurvivalTime = c(0, 6),
#'         lambda1 = c(0.0533, 0.0309), lambda2 = c(0.0533, 0.0533),
#'         gamma1 = 0.00427, gamma2 = 0.00427,
#'         accrualDuration = 22, followupTime = 18)
#'
#' @useDynLib lrstat, .registration = TRUE
#' @importFrom Rcpp evalCpp
#'
NULL

