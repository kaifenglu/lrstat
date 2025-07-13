#' @title Group Sequential Design for One-Sample Mean
#' @description Obtains the power given sample size or obtains the sample
#' size given power for a group sequential design for one-sample mean.
#'
#' @param beta The type II error.
#' @param n The total sample size.
#' @param meanH0 The mean under the null hypothesis.
#'   Defaults to 0.
#' @param mean The mean under the alternative hypothesis.
#' @param stDev The standard deviation.
#' @param normalApproximation The type of computation of the p-values.
#'   If \code{TRUE}, the variance is assumed to be known, otherwise
#'   the calculations are performed with the t distribution. The exact
#'   calculation using the t distribution is only implemented for the
#'   fixed design.
#' @param rounding Whether to round up sample size. Defaults to 1 for
#'   sample size rounding.
#' @inheritParams param_kMax
#' @param informationRates The information rates. Fixed prior to the trial.
#'   Defaults to \code{(1:kMax) / kMax} if left unspecified.
#' @inheritParams param_efficacyStopping
#' @inheritParams param_futilityStopping
#' @inheritParams param_criticalValues
#' @inheritParams param_alpha
#' @inheritParams param_typeAlphaSpending
#' @inheritParams param_parameterAlphaSpending
#' @inheritParams param_userAlphaSpending
#' @inheritParams param_futilityBounds
#' @inheritParams param_typeBetaSpending
#' @inheritParams param_parameterBetaSpending
#' @inheritParams param_userBetaSpending
#' @param spendingTime A vector of length \code{kMax} for the error spending
#'   time at each analysis. Defaults to missing, in which case, it is the
#'   same as \code{informationRates}.
#'
#' @return An S3 class \code{designOneMean} object with three components:
#'
#' * \code{overallResults}: A data frame containing the following variables:
#'
#'     - \code{overallReject}: The overall rejection probability.
#'
#'     - \code{alpha}: The overall significance level.
#'
#'     - \code{attainedAlpha}: The attained significance level, which is
#'       different from the overall significance level in the presence of
#'       futility stopping.
#'
#'     - \code{kMax}: The number of stages.
#'
#'     - \code{theta}: The parameter value.
#'
#'     - \code{information}: The maximum information.
#'
#'     - \code{expectedInformationH1}: The expected information under H1.
#'
#'     - \code{expectedInformationH0}: The expected information under H0.
#'
#'     - \code{drift}: The drift parameter, equal to
#'       \code{theta*sqrt(information)}.
#'
#'     - \code{inflationFactor}: The inflation factor (relative to the
#'       fixed design).
#'
#'     - \code{numberOfSubjects}: The maximum number of subjects.
#'
#'     - \code{expectedNumberOfSubjectsH1}: The expected number of subjects
#'       under H1.
#'
#'     - \code{expectedNumberOfSubjectsH0}: The expected number of subjects
#'       under H0.
#'
#'     - \code{meanH0}: The mean under the null hypothesis.
#'
#'     - \code{mean}: The mean under the alternative hypothesis.
#'
#'     - \code{stDev}: The standard deviation.
#'
#' * \code{byStageResults}: A data frame containing the following variables:
#'
#'     - \code{informationRates}: The information rates.
#'
#'     - \code{efficacyBounds}: The efficacy boundaries on the Z-scale.
#'
#'     - \code{futilityBounds}: The futility boundaries on the Z-scale.
#'
#'     - \code{rejectPerStage}: The probability for efficacy stopping.
#'
#'     - \code{futilityPerStage}: The probability for futility stopping.
#'
#'     - \code{cumulativeRejection}: The cumulative probability for efficacy
#'       stopping.
#'
#'     - \code{cumulativeFutility}: The cumulative probability for futility
#'       stopping.
#'
#'     - \code{cumulativeAlphaSpent}: The cumulative alpha spent.
#'
#'     - \code{efficacyP}: The efficacy boundaries on the p-value scale.
#'
#'     - \code{futilityP}: The futility boundaries on the p-value scale.
#'
#'     - \code{information}: The cumulative information.
#'
#'     - \code{efficacyStopping}: Whether to allow efficacy stopping.
#'
#'     - \code{futilityStopping}: Whether to allow futility stopping.
#'
#'     - \code{rejectPerStageH0}: The probability for efficacy stopping
#'       under H0.
#'
#'     - \code{futilityPerStageH0}: The probability for futility stopping
#'       under H0.
#'
#'     - \code{cumulativeRejectionH0}: The cumulative probability for
#'       efficacy stopping under H0.
#'
#'     - \code{cumulativeFutilityH0}: The cumulative probability for futility
#'       stopping under H0.
#'
#'     - \code{efficacyMean}: The efficacy boundaries on the mean scale.
#'
#'     - \code{futilityMean}: The futility boundaries on the mean scale.
#'
#'     - \code{numberOfSubjects}: The number of subjects.
#'
#' * \code{settings}: A list containing the following input parameters:
#'
#'     - \code{typeAlphaSpending}: The type of alpha spending.
#'
#'     - \code{parameterAlphaSpending}: The parameter value for alpha
#'       spending.
#'
#'     - \code{userAlphaSpending}: The user defined alpha spending.
#'
#'     - \code{typeBetaSpending}: The type of beta spending.
#'
#'     - \code{parameterBetaSpending}: The parameter value for beta spending.
#'
#'     - \code{userBetaSpending}: The user defined beta spending.
#'
#'     - \code{spendingTime}: The error spending time at each analysis.
#'
#'     - \code{normalApproximation}: The type of computation of the p-values.
#'       If \code{TRUE}, the variance is assumed to be known, otherwise
#'       the calculations are performed with the t distribution.
#'
#'     - \code{rounding}: Whether to round up sample size.
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @examples
#'
#' # Example 1: group sequential trial power calculation
#' (design1 <- getDesignOneMean(
#'   beta = 0.1, n = NA, meanH0 = 7, mean = 6, stDev = 2.5,
#'   kMax = 5, alpha = 0.025, typeAlphaSpending = "sfOF",
#'   typeBetaSpending = "sfP"))
#'
#' # Example 2: sample size calculation for one-sample t-test
#' (design2 <- getDesignOneMean(
#'   beta = 0.1, n = NA, meanH0 = 7, mean = 6, stDev = 2.5,
#'   normalApproximation = FALSE, alpha = 0.025))
#'
#' @export
getDesignOneMean <- function(
    beta = NA_real_,
    n = NA_real_,
    meanH0 = 0,
    mean = 0.5,
    stDev = 1,
    normalApproximation = TRUE,
    rounding = TRUE,
    kMax = 1L,
    informationRates = NA_real_,
    efficacyStopping = NA_integer_,
    futilityStopping = NA_integer_,
    criticalValues = NA_real_,
    alpha = 0.025,
    typeAlphaSpending = "sfOF",
    parameterAlphaSpending = NA_real_,
    userAlphaSpending = NA_real_,
    futilityBounds = NA_real_,
    typeBetaSpending = "none",
    parameterBetaSpending = NA_real_,
    userBetaSpending = NA_real_,
    spendingTime = NA_real_) {

  if (is.na(beta) && is.na(n)) {
    stop("beta and n cannot be both missing")
  }

  if (!is.na(beta) && !is.na(n)) {
    stop("Only one of beta and n should be provided")
  }

  if (!is.na(beta) && (beta >= 1-alpha || beta < 0.0001)) {
    stop("beta must lie in [0.0001, 1-alpha)")
  }

  if (!is.na(n) && n <= 0) {
    stop("n must be positive")
  }

  if (stDev <= 0) {
    stop("stDev must be positive")
  }

  if (alpha < 0.00001 || alpha >= 1) {
    stop("alpha must lie in [0.00001, 1)")
  }

  if (any(is.na(informationRates))) {
    informationRates = (1:kMax)/kMax
  }

  directionUpper = mean > meanH0

  theta = ifelse(directionUpper, mean - meanH0, meanH0 - mean)

  # variance for one sampling unit
  v1 = stDev^2

  if (is.na(beta)) { # power calculation
    if (rounding) {
      n = ceiling(n - 1.0e-12)
      informationRates = round(n*informationRates)/n
    }

    des = getDesign(
      beta = NA, IMax = n/v1, theta,
      kMax, informationRates,
      efficacyStopping, futilityStopping,
      criticalValues, alpha, typeAlphaSpending,
      parameterAlphaSpending, userAlphaSpending,
      futilityBounds, typeBetaSpending,
      parameterBetaSpending, userBetaSpending,
      spendingTime, 1)

    if (kMax == 1 && !normalApproximation) { # t-test for fixed design
      if (!any(is.na(criticalValues))) {
        alpha = 1 - pnorm(criticalValues)
      }

      b = qt(1-alpha, n-1)
      ncp = theta*sqrt(n/v1)
      power = pt(b, n-1, ncp, lower.tail = FALSE)

      delta = b/sqrt(n/v1)

      des$overallResults$overallReject = power
      des$byStageResults$rejectPerStage = power
      des$byStageResults$futilityPerStage = 1 - power
      des$byStageResults$cumulativeRejection = power
      des$byStageResults$cumulativeFutility = 1 - power
      des$byStageResults$efficacyBounds = b
      des$byStageResults$futilityBounds = b

      if (directionUpper) {
        des$byStageResults$efficacyMean = delta + meanH0
        des$byStageResults$futilityMean = delta + meanH0
      } else {
        des$byStageResults$efficacyMean = -delta + meanH0
        des$byStageResults$futilityMean = -delta + meanH0
      }
    } else {
      if (directionUpper) {
        des$byStageResults$efficacyMean =
          des$byStageResults$efficacyTheta + meanH0
        des$byStageResults$futilityMean =
          des$byStageResults$futilityTheta + meanH0
      } else {
        des$byStageResults$efficacyMean =
          -des$byStageResults$efficacyTheta + meanH0
        des$byStageResults$futilityMean =
          -des$byStageResults$futilityTheta + meanH0
      }
    }
  } else { # sample size calculation
    if (kMax == 1 && !normalApproximation) { # t-test for fixed design
      if (!any(is.na(criticalValues))) {
        alpha = 1 - pnorm(criticalValues)
      }

      n0 = (qnorm(1-alpha) + qnorm(1-beta))^2*v1/theta^2

      n = uniroot(function(n) {
        b = qt(1-alpha, n-1)
        ncp = theta*sqrt(n/v1)
        pt(b, n-1, ncp, lower.tail = FALSE) - (1 - beta)
      }, c(n0, 2*n0))$root

      if (rounding) n = ceiling(n - 1.0e-12)

      des = getDesign(
        beta = NA, IMax = n/v1, theta,
        kMax, informationRates,
        efficacyStopping, futilityStopping,
        criticalValues, alpha, typeAlphaSpending,
        parameterAlphaSpending, userAlphaSpending,
        futilityBounds, typeBetaSpending,
        parameterBetaSpending, userBetaSpending,
        spendingTime, 1)

      b = qt(1-alpha, n-1)
      ncp = theta*sqrt(n/v1)
      power = pt(b, n-1, ncp, lower.tail = FALSE)

      delta = b/sqrt(n/v1)

      des$overallResults$overallReject = power
      des$byStageResults$rejectPerStage = power
      des$byStageResults$futilityPerStage = 1 - power
      des$byStageResults$cumulativeRejection = power
      des$byStageResults$cumulativeFutility = 1 - power
      des$byStageResults$efficacyBounds = b
      des$byStageResults$futilityBounds = b
      if (directionUpper) {
        des$byStageResults$efficacyMean = delta + meanH0
        des$byStageResults$futilityMean = delta + meanH0
      } else {
        des$byStageResults$efficacyMean = -delta + meanH0
        des$byStageResults$futilityMean = -delta + meanH0
      }
    } else {
      des = getDesign(
        beta, IMax = NA, theta,
        kMax, informationRates,
        efficacyStopping, futilityStopping,
        criticalValues, alpha, typeAlphaSpending,
        parameterAlphaSpending, userAlphaSpending,
        futilityBounds, typeBetaSpending,
        parameterBetaSpending, userBetaSpending,
        spendingTime, 1)

      n = des$overallResults$information*v1

      if (rounding) {
        n = ceiling(n - 1.0e-12)
        informationRates = des$byStageResults$informationRates
        informationRates = round(n*informationRates)/n

        des = getDesign(
          beta = NA, IMax = n/v1, theta,
          kMax, informationRates,
          efficacyStopping, futilityStopping,
          criticalValues, alpha, typeAlphaSpending,
          parameterAlphaSpending, userAlphaSpending,
          futilityBounds, typeBetaSpending,
          parameterBetaSpending, userBetaSpending,
          spendingTime, 1)
      }

      if (directionUpper) {
        des$byStageResults$efficacyMean =
          des$byStageResults$efficacyTheta + meanH0
        des$byStageResults$futilityMean =
          des$byStageResults$futilityTheta + meanH0
      } else {
        des$byStageResults$efficacyMean =
          -des$byStageResults$efficacyTheta + meanH0
        des$byStageResults$futilityMean =
          -des$byStageResults$futilityTheta + meanH0
      }
    }
  }

  des$overallResults$theta = theta
  des$overallResults$numberOfSubjects = n
  des$overallResults$expectedNumberOfSubjectsH1 =
    des$overallResults$expectedInformationH1*v1
  des$overallResults$expectedNumberOfSubjectsH0 =
    des$overallResults$expectedInformationH0*v1
  des$overallResults$meanH0 = meanH0
  des$overallResults$mean = mean
  des$overallResults$stDev = stDev

  des$byStageResults$efficacyTheta = NULL
  des$byStageResults$futilityTheta = NULL
  des$byStageResults$numberOfSubjects =
    des$byStageResults$informationRates*n

  des$settings$normalApproximation = normalApproximation
  des$settings$rounding = rounding
  des$settings$varianceRatio = NULL

  attr(des, "class") = "designOneMean"

  des
}

#' @title Group Sequential Design for Paired Mean Difference
#' @description Obtains the power given sample size or obtains the sample
#' size given power for a group sequential design for paired mean
#' difference.
#'
#' @param beta The type II error.
#' @param n The total sample size.
#' @param pairedDiffH0 The paired difference under the null hypothesis.
#'   Defaults to 0.
#' @param pairedDiff The paired difference under the alternative hypothesis.
#' @param stDev The standard deviation for paired difference.
#' @param normalApproximation The type of computation of the p-values.
#'   If \code{TRUE}, the variance is assumed to be known, otherwise
#'   the calculations are performed with the t distribution. The exact
#'   calculation using the t distribution is only implemented for the
#'   fixed design.
#' @param rounding Whether to round up sample size. Defaults to 1 for
#' sample size rounding.
#' @inheritParams param_kMax
#' @param informationRates The information rates. Fixed prior to the trial.
#'   Defaults to \code{(1:kMax) / kMax} if left unspecified.
#' @inheritParams param_efficacyStopping
#' @inheritParams param_futilityStopping
#' @inheritParams param_criticalValues
#' @inheritParams param_alpha
#' @inheritParams param_typeAlphaSpending
#' @inheritParams param_parameterAlphaSpending
#' @inheritParams param_userAlphaSpending
#' @inheritParams param_futilityBounds
#' @inheritParams param_typeBetaSpending
#' @inheritParams param_parameterBetaSpending
#' @inheritParams param_userBetaSpending
#' @param spendingTime A vector of length \code{kMax} for the error spending
#'   time at each analysis. Defaults to missing, in which case, it is the
#'   same as \code{informationRates}.
#'
#' @return An S3 class \code{designPairedMeanDiff} object with three
#' components:
#'
#' * \code{overallResults}: A data frame containing the following variables:
#'
#'     - \code{overallReject}: The overall rejection probability.
#'
#'     - \code{alpha}: The overall significance level.
#'
#'     - \code{attainedAlpha}: The attained significance level, which is
#'       different from the overall significance level in the presence of
#'       futility stopping.
#'
#'     - \code{kMax}: The number of stages.
#'
#'     - \code{theta}: The parameter value.
#'
#'     - \code{information}: The maximum information.
#'
#'     - \code{expectedInformationH1}: The expected information under H1.
#'
#'     - \code{expectedInformationH0}: The expected information under H0.
#'
#'     - \code{drift}: The drift parameter, equal to
#'       \code{theta*sqrt(information)}.
#'
#'     - \code{inflationFactor}: The inflation factor (relative to the
#'       fixed design).
#'
#'     - \code{numberOfSubjects}: The maximum number of subjects.
#'
#'     - \code{expectedNumberOfSubjectsH1}: The expected number of subjects
#'       under H1.
#'
#'     - \code{expectedNumberOfSubjectsH0}: The expected number of subjects
#'       under H0.
#'
#'     - \code{pairedDiffH0}: The paired difference under the null
#'       hypothesis.
#'
#'     - \code{pairedDiff}: The paired difference under the alternative
#'       hypothesis.
#'
#'     - \code{stDev}: The standard deviation for paired difference.
#'
#' * \code{byStageResults}: A data frame containing the following variables:
#'
#'     - \code{informationRates}: The information rates.
#'
#'     - \code{efficacyBounds}: The efficacy boundaries on the Z-scale.
#'
#'     - \code{futilityBounds}: The futility boundaries on the Z-scale.
#'
#'     - \code{rejectPerStage}: The probability for efficacy stopping.
#'
#'     - \code{futilityPerStage}: The probability for futility stopping.
#'
#'     - \code{cumulativeRejection}: The cumulative probability for efficacy
#'       stopping.
#'
#'     - \code{cumulativeFutility}: The cumulative probability for futility
#'       stopping.
#'
#'     - \code{cumulativeAlphaSpent}: The cumulative alpha spent.
#'
#'     - \code{efficacyP}: The efficacy boundaries on the p-value scale.
#'
#'     - \code{futilityP}: The futility boundaries on the p-value scale.
#'
#'     - \code{information}: The cumulative information.
#'
#'     - \code{efficacyStopping}: Whether to allow efficacy stopping.
#'
#'     - \code{futilityStopping}: Whether to allow futility stopping.
#'
#'     - \code{rejectPerStageH0}: The probability for efficacy stopping
#'       under H0.
#'
#'     - \code{futilityPerStageH0}: The probability for futility stopping
#'       under H0.
#'
#'     - \code{cumulativeRejectionH0}: The cumulative probability for
#'       efficacy stopping under H0.
#'
#'     - \code{cumulativeFutilityH0}: The cumulative probability for futility
#'       stopping under H0.
#'
#'     - \code{efficacyPairedDiff}: The efficacy boundaries on the paired
#'       difference scale.
#'
#'     - \code{futilityPairedDiff}: The futility boundaries on the paired
#'       difference scale.
#'
#'     - \code{numberOfSubjects}: The number of subjects.
#'
#' * \code{settings}: A list containing the following input parameters:
#'
#'     - \code{typeAlphaSpending}: The type of alpha spending.
#'
#'     - \code{parameterAlphaSpending}: The parameter value for alpha
#'       spending.
#'
#'     - \code{userAlphaSpending}: The user defined alpha spending.
#'
#'     - \code{typeBetaSpending}: The type of beta spending.
#'
#'     - \code{parameterBetaSpending}: The parameter value for beta spending.
#'
#'     - \code{userBetaSpending}: The user defined beta spending.
#'
#'     - \code{spendingTime}: The error spending time at each analysis.
#'
#'     - \code{normalApproximation}: The type of computation of the p-values.
#'       If \code{TRUE}, the variance is assumed to be known, otherwise
#'       the calculations are performed with the t distribution.
#'
#'     - \code{rounding}: Whether to round up sample size.
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @examples
#'
#' # Example 1: group sequential trial power calculation
#' (design1 <- getDesignPairedMeanDiff(
#'   beta = 0.1, n = NA, pairedDiffH0 = 0, pairedDiff = -2, stDev = 5,
#'   kMax = 5, alpha = 0.05, typeAlphaSpending = "sfOF"))
#'
#' # Example 2: sample size calculation for one-sample t-test
#' (design2 <- getDesignPairedMeanDiff(
#'   beta = 0.1, n = NA, pairedDiffH0 = 0, pairedDiff = -2, stDev = 5,
#'   normalApproximation = FALSE, alpha = 0.025))
#'
#' @export
getDesignPairedMeanDiff <- function(
    beta = NA_real_,
    n = NA_real_,
    pairedDiffH0 = 0,
    pairedDiff = 0.5,
    stDev = 1,
    normalApproximation = TRUE,
    rounding = TRUE,
    kMax = 1L,
    informationRates = NA_real_,
    efficacyStopping = NA_integer_,
    futilityStopping = NA_integer_,
    criticalValues = NA_real_,
    alpha = 0.025,
    typeAlphaSpending = "sfOF",
    parameterAlphaSpending = NA_real_,
    userAlphaSpending = NA_real_,
    futilityBounds = NA_real_,
    typeBetaSpending = "none",
    parameterBetaSpending = NA_real_,
    userBetaSpending = NA_real_,
    spendingTime = NA_real_) {

  des = getDesignOneMean(
    beta, n, pairedDiffH0, pairedDiff, stDev,
    normalApproximation,
    rounding, kMax, informationRates,
    efficacyStopping, futilityStopping,
    criticalValues, alpha, typeAlphaSpending,
    parameterAlphaSpending, userAlphaSpending,
    futilityBounds, typeBetaSpending,
    parameterBetaSpending, userBetaSpending,
    spendingTime)

  nov = names(des$overallResults)
  names(des$overallResults)[nov == "meanH0"] <- "pairedDiffH0"
  names(des$overallResults)[nov == "mean"] <- "pairedDiff"

  nby = names(des$byStageResults)
  names(des$byStageResults)[nby == "efficacyMean"] <- "efficacyPairedDiff"
  names(des$byStageResults)[nby == "futilityMean"] <- "futilityPairedDiff"

  attr(des, "class") = "designPairedMeanDiff"

  des
}


#' @title Group Sequential Design for Paired Mean Ratio
#' @description Obtains the power given sample size or obtains the sample
#' size given power for a group sequential design for paired mean ratio.
#'
#' @param beta The type II error.
#' @param n The total sample size.
#' @param pairedRatioH0 The paired ratio under the null hypothesis.
#' @param pairedRatio The paired ratio under the alternative
#'   hypothesis.
#' @param CV The coefficient of variation for paired ratio.
#' @param normalApproximation The type of computation of the p-values.
#'   If \code{TRUE}, the variance is assumed to be known, otherwise
#'   the calculations are performed with the t distribution. The exact
#'   calculation using the t distribution is only implemented for the
#'   fixed design.
#' @param rounding Whether to round up sample size. Defaults to 1 for
#'   sample size rounding.
#' @inheritParams param_kMax
#' @param informationRates The information rates. Fixed prior to the trial.
#'   Defaults to \code{(1:kMax) / kMax} if left unspecified.
#' @inheritParams param_efficacyStopping
#' @inheritParams param_futilityStopping
#' @inheritParams param_criticalValues
#' @inheritParams param_alpha
#' @inheritParams param_typeAlphaSpending
#' @inheritParams param_parameterAlphaSpending
#' @inheritParams param_userAlphaSpending
#' @inheritParams param_futilityBounds
#' @inheritParams param_typeBetaSpending
#' @inheritParams param_parameterBetaSpending
#' @inheritParams param_userBetaSpending
#' @param spendingTime A vector of length \code{kMax} for the error spending
#'   time at each analysis. Defaults to missing, in which case, it is the
#'   same as \code{informationRates}.
#'
#' @return An S3 class \code{designPairedMeanRatio} object with three
#' components:
#'
#' * \code{overallResults}: A data frame containing the following variables:
#'
#'     - \code{overallReject}: The overall rejection probability.
#'
#'     - \code{alpha}: The overall significance level.
#'
#'     - \code{attainedAlpha}: The attained significance level, which is
#'       different from the overall significance level in the presence of
#'       futility stopping.
#'
#'     - \code{kMax}: The number of stages.
#'
#'     - \code{theta}: The parameter value.
#'
#'     - \code{information}: The maximum information.
#'
#'     - \code{expectedInformationH1}: The expected information under H1.
#'
#'     - \code{expectedInformationH0}: The expected information under H0.
#'
#'     - \code{drift}: The drift parameter, equal to
#'       \code{theta*sqrt(information)}.
#'
#'     - \code{inflationFactor}: The inflation factor (relative to the
#'       fixed design).
#'
#'     - \code{numberOfSubjects}: The maximum number of subjects.
#'
#'     - \code{expectedNumberOfSubjectsH1}: The expected number of subjects
#'       under H1.
#'
#'     - \code{expectedNumberOfSubjectsH0}: The expected number of subjects
#'       under H0.
#'
#'     - \code{pairedRatioH0}: The paired ratio under the null hypothesis.
#'
#'     - \code{pairedRatio}: The paired ratio under the alternative
#'       hypothesis.
#'
#'     - \code{CV}: The coefficient of variation for paired ratio.
#'
#' * \code{byStageResults}: A data frame containing the following variables:
#'
#'     - \code{informationRates}: The information rates.
#'
#'     - \code{efficacyBounds}: The efficacy boundaries on the Z-scale.
#'
#'     - \code{futilityBounds}: The futility boundaries on the Z-scale.
#'
#'     - \code{rejectPerStage}: The probability for efficacy stopping.
#'
#'     - \code{futilityPerStage}: The probability for futility stopping.
#'
#'     - \code{cumulativeRejection}: The cumulative probability for efficacy
#'       stopping.
#'
#'     - \code{cumulativeFutility}: The cumulative probability for futility
#'       stopping.
#'
#'     - \code{cumulativeAlphaSpent}: The cumulative alpha spent.
#'
#'     - \code{efficacyP}: The efficacy boundaries on the p-value scale.
#'
#'     - \code{futilityP}: The futility boundaries on the p-value scale.
#'
#'     - \code{information}: The cumulative information.
#'
#'     - \code{efficacyStopping}: Whether to allow efficacy stopping.
#'
#'     - \code{futilityStopping}: Whether to allow futility stopping.
#'
#'     - \code{rejectPerStageH0}: The probability for efficacy stopping
#'       under H0.
#'
#'     - \code{futilityPerStageH0}: The probability for futility stopping
#'       under H0.
#'
#'     - \code{cumulativeRejectionH0}: The cumulative probability for
#'       efficacy stopping under H0.
#'
#'     - \code{cumulativeFutilityH0}: The cumulative probability for futility
#'       stopping under H0.
#'
#'     - \code{numberOfSubjects}: The number of subjects.
#'
#'     - \code{efficacyPairedRatio}: The efficacy boundaries on the paired
#'       ratio scale.
#'
#'     - \code{futilityPairedRatio}: The futility boundaries on the paired
#'       ratio scale.
#'
#' * \code{settings}: A list containing the following input parameters:
#'
#'     - \code{typeAlphaSpending}: The type of alpha spending.
#'
#'     - \code{parameterAlphaSpending}: The parameter value for alpha
#'       spending.
#'
#'     - \code{userAlphaSpending}: The user defined alpha spending.
#'
#'     - \code{typeBetaSpending}: The type of beta spending.
#'
#'     - \code{parameterBetaSpending}: The parameter value for beta spending.
#'
#'     - \code{userBetaSpending}: The user defined beta spending.
#'
#'     - \code{spendingTime}: The error spending time at each analysis.
#'
#'     - \code{normalApproximation}: The type of computation of the p-values.
#'       If \code{TRUE}, the variance is assumed to be known, otherwise
#'       the calculations are performed with the t distribution.
#'
#'     - \code{rounding}: Whether to round up sample size.
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @examples
#'
#' # Example 1: group sequential trial power calculation
#' (design1 <- getDesignPairedMeanRatio(
#'   beta = 0.1, n = NA, pairedRatio = 1.2, CV = 0.35,
#'   kMax = 5, alpha = 0.05, typeAlphaSpending = "sfOF"))
#'
#' # Example 2: sample size calculation for one-sample t-test
#' (design2 <- getDesignPairedMeanRatio(
#'   beta = 0.1, n = NA, pairedRatio = 1.2, CV = 0.35,
#'   normalApproximation = FALSE, alpha = 0.05))
#'
#' @export
getDesignPairedMeanRatio <- function(
    beta = NA_real_,
    n = NA_real_,
    pairedRatioH0 = 1,
    pairedRatio = 1.2,
    CV = 1,
    normalApproximation = TRUE,
    rounding = TRUE,
    kMax = 1L,
    informationRates = NA_real_,
    efficacyStopping = NA_integer_,
    futilityStopping = NA_integer_,
    criticalValues = NA_real_,
    alpha = 0.025,
    typeAlphaSpending = "sfOF",
    parameterAlphaSpending = NA_real_,
    userAlphaSpending = NA_real_,
    futilityBounds = NA_real_,
    typeBetaSpending = "none",
    parameterBetaSpending = NA_real_,
    userBetaSpending = NA_real_,
    spendingTime = NA_real_) {

  if (pairedRatioH0 <= 0) {
    stop("pairedRatioH0 must be positive")
  }

  if (pairedRatio <= 0) {
    stop("pairedRatio must be positive")
  }

  if (CV <= 0) {
    stop("CV must be positive")
  }

  if (alpha < 0.00001 || alpha >= 1) {
    stop("alpha must lie in [0.00001, 1)")
  }

  des = getDesignPairedMeanDiff(
    beta, n, pairedDiffH0 = log(pairedRatioH0),
    pairedDiff = log(pairedRatio),
    stDev = sqrt(log(1 + CV^2)),
    normalApproximation,
    rounding, kMax, informationRates,
    efficacyStopping, futilityStopping,
    criticalValues, alpha, typeAlphaSpending,
    parameterAlphaSpending, userAlphaSpending,
    futilityBounds, typeBetaSpending,
    parameterBetaSpending, userBetaSpending,
    spendingTime)

  des$overallResults$pairedRatioH0 = pairedRatioH0
  des$overallResults$pairedRatio = pairedRatio
  des$overallResults$CV = CV
  des$overallResults$pairedDiffH0 = NULL
  des$overallResults$pairedDiff = NULL
  des$overallResults$stDev = NULL

  des$byStageResults$efficacyPairedRatio =
    exp(des$byStageResults$efficacyPairedDiff)
  des$byStageResults$futilityPairedRatio =
    exp(des$byStageResults$futilityPairedDiff)
  des$byStageResults$efficacyPairedDiff = NULL
  des$byStageResults$futilityPairedDiff = NULL

  attr(des, "class") = "designPairedMeanRatio"

  des
}


#' @title Group Sequential Design for Two-Sample Mean Difference
#' @description Obtains the power given sample size or obtains the sample
#' size given power for a group sequential design for two-sample mean
#' difference.
#'
#' @param beta The type II error.
#' @param n The total sample size.
#' @param meanDiffH0 The mean difference under the null hypothesis.
#'   Defaults to 0.
#' @param meanDiff The mean difference under the alternative hypothesis.
#' @param stDev The standard deviation.
#' @param allocationRatioPlanned Allocation ratio for the active treatment
#'   versus control. Defaults to 1 for equal randomization.
#' @param normalApproximation The type of computation of the p-values.
#'   If \code{TRUE}, the variance is assumed to be known, otherwise
#'   the calculations are performed with the t distribution. The exact
#'   calculation using the t distribution is only implemented for the
#'   fixed design.
#' @param rounding Whether to round up sample size. Defaults to 1 for
#'   sample size rounding.
#' @inheritParams param_kMax
#' @param informationRates The information rates. Fixed prior to the trial.
#'   Defaults to \code{(1:kMax) / kMax} if left unspecified.
#' @inheritParams param_efficacyStopping
#' @inheritParams param_futilityStopping
#' @inheritParams param_criticalValues
#' @inheritParams param_alpha
#' @inheritParams param_typeAlphaSpending
#' @inheritParams param_parameterAlphaSpending
#' @inheritParams param_userAlphaSpending
#' @inheritParams param_futilityBounds
#' @inheritParams param_typeBetaSpending
#' @inheritParams param_parameterBetaSpending
#' @inheritParams param_userBetaSpending
#' @param spendingTime A vector of length \code{kMax} for the error spending
#'   time at each analysis. Defaults to missing, in which case, it is the
#'   same as \code{informationRates}.
#'
#' @return An S3 class \code{designMeanDiff} object with three components:
#'
#' * \code{overallResults}: A data frame containing the following variables:
#'
#'     - \code{overallReject}: The overall rejection probability.
#'
#'     - \code{alpha}: The overall significance level.
#'
#'     - \code{attainedAlpha}: The attained significance level, which is
#'       different from the overall significance level in the presence of
#'       futility stopping.
#'
#'     - \code{kMax}: The number of stages.
#'
#'     - \code{theta}: The parameter value.
#'
#'     - \code{information}: The maximum information.
#'
#'     - \code{expectedInformationH1}: The expected information under H1.
#'
#'     - \code{expectedInformationH0}: The expected information under H0.
#'
#'     - \code{drift}: The drift parameter, equal to
#'       \code{theta*sqrt(information)}.
#'
#'     - \code{inflationFactor}: The inflation factor (relative to the
#'       fixed design).
#'
#'     - \code{numberOfSubjects}: The maximum number of subjects.
#'
#'     - \code{expectedNumberOfSubjectsH1}: The expected number of subjects
#'       under H1.
#'
#'     - \code{expectedNumberOfSubjectsH0}: The expected number of subjects
#'       under H0.
#'
#'     - \code{meanDiffH0}: The mean difference under the null hypothesis.
#'
#'     - \code{meanDiff}: The mean difference under the alternative
#'       hypothesis.
#'
#'     - \code{stDev}: The standard deviation.
#'
#' * \code{byStageResults}: A data frame containing the following variables:
#'
#'     - \code{informationRates}: The information rates.
#'
#'     - \code{efficacyBounds}: The efficacy boundaries on the Z-scale.
#'
#'     - \code{futilityBounds}: The futility boundaries on the Z-scale.
#'
#'     - \code{rejectPerStage}: The probability for efficacy stopping.
#'
#'     - \code{futilityPerStage}: The probability for futility stopping.
#'
#'     - \code{cumulativeRejection}: The cumulative probability for efficacy
#'       stopping.
#'
#'     - \code{cumulativeFutility}: The cumulative probability for futility
#'       stopping.
#'
#'     - \code{cumulativeAlphaSpent}: The cumulative alpha spent.
#'
#'     - \code{efficacyP}: The efficacy boundaries on the p-value scale.
#'
#'     - \code{futilityP}: The futility boundaries on the p-value scale.
#'
#'     - \code{information}: The cumulative information.
#'
#'     - \code{efficacyStopping}: Whether to allow efficacy stopping.
#'
#'     - \code{futilityStopping}: Whether to allow futility stopping.
#'
#'     - \code{rejectPerStageH0}: The probability for efficacy stopping
#'       under H0.
#'
#'     - \code{futilityPerStageH0}: The probability for futility stopping
#'       under H0.
#'
#'     - \code{cumulativeRejectionH0}: The cumulative probability for
#'       efficacy stopping under H0.
#'
#'     - \code{cumulativeFutilityH0}: The cumulative probability for futility
#'       stopping under H0.
#'
#'     - \code{efficacyMeanDiff}: The efficacy boundaries on the mean
#'       difference scale.
#'
#'     - \code{futilityMeanDiff}: The futility boundaries on the mean
#'       difference scale.
#'
#'     - \code{numberOfSubjects}: The number of subjects.
#'
#' * \code{settings}: A list containing the following input parameters:
#'
#'     - \code{typeAlphaSpending}: The type of alpha spending.
#'
#'     - \code{parameterAlphaSpending}: The parameter value for alpha
#'       spending.
#'
#'     - \code{userAlphaSpending}: The user defined alpha spending.
#'
#'     - \code{typeBetaSpending}: The type of beta spending.
#'
#'     - \code{parameterBetaSpending}: The parameter value for beta spending.
#'
#'     - \code{userBetaSpending}: The user defined beta spending.
#'
#'     - \code{spendingTime}: The error spending time at each analysis.
#'
#'     - \code{allocationRatioPlanned}: Allocation ratio for the active
#'       treatment versus control.
#'
#'     - \code{normalApproximation}: The type of computation of the p-values.
#'       If \code{TRUE}, the variance is assumed to be known, otherwise
#'       the calculations are performed with the t distribution.
#'
#'     - \code{rounding}: Whether to round up sample size.
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @examples
#'
#' # Example 1: group sequential trial power calculation
#' (design1 <- getDesignMeanDiff(
#'   beta = NA, n = 456, meanDiff = 9, stDev = 32,
#'   kMax = 5, alpha = 0.025, typeAlphaSpending = "sfOF",
#'   typeBetaSpending = "sfP"))
#'
#' # Example 2: sample size calculation for two-sample t-test
#' (design2 <- getDesignMeanDiff(
#'   beta = 0.1, n = NA, meanDiff = 0.3, stDev = 1,
#'   normalApproximation = FALSE, alpha = 0.025))
#'
#' @export
getDesignMeanDiff <- function(
    beta = NA_real_,
    n = NA_real_,
    meanDiffH0 = 0,
    meanDiff = 0.5,
    stDev = 1,
    allocationRatioPlanned = 1,
    normalApproximation = TRUE,
    rounding = TRUE,
    kMax = 1L,
    informationRates = NA_real_,
    efficacyStopping = NA_integer_,
    futilityStopping = NA_integer_,
    criticalValues = NA_real_,
    alpha = 0.025,
    typeAlphaSpending = "sfOF",
    parameterAlphaSpending = NA_real_,
    userAlphaSpending = NA_real_,
    futilityBounds = NA_real_,
    typeBetaSpending = "none",
    parameterBetaSpending = NA_real_,
    userBetaSpending = NA_real_,
    spendingTime = NA_real_) {

  if (is.na(beta) && is.na(n)) {
    stop("beta and n cannot be both missing")
  }

  if (!is.na(beta) && !is.na(n)) {
    stop("Only one of beta and n should be provided")
  }

  if (!is.na(beta) && (beta >= 1-alpha || beta < 0.0001)) {
    stop("beta must lie in [0.0001, 1-alpha)")
  }

  if (!is.na(n) && n <= 0) {
    stop("n must be positive")
  }

  if (stDev <= 0) {
    stop("stDev must be positive")
  }

  if (allocationRatioPlanned <= 0) {
    stop("allocationRatioPlanned must be positive")
  }

  if (alpha < 0.00001 || alpha >= 1) {
    stop("alpha must lie in [0.00001, 1)")
  }

  if (any(is.na(informationRates))) {
    informationRates = (1:kMax)/kMax
  }


  r = allocationRatioPlanned/(1 + allocationRatioPlanned)

  directionUpper = meanDiff > meanDiffH0

  theta = ifelse(directionUpper, meanDiff - meanDiffH0,
                 meanDiffH0 - meanDiff)

  # variance for one sampling unit
  v1 = stDev^2/(r*(1-r))

  if (is.na(beta)) { # power calculation
    if (rounding) {
      n = ceiling(n - 1.0e-12)
      informationRates = round(n*informationRates)/n
    }

    des = getDesign(
      beta = NA, IMax = n/v1, theta,
      kMax, informationRates,
      efficacyStopping, futilityStopping,
      criticalValues, alpha, typeAlphaSpending,
      parameterAlphaSpending, userAlphaSpending,
      futilityBounds, typeBetaSpending,
      parameterBetaSpending, userBetaSpending,
      spendingTime, 1)

    if (kMax == 1 && !normalApproximation) { # t-test for fixed design
      if (!any(is.na(criticalValues))) {
        alpha = 1 - pnorm(criticalValues)
      }

      b = qt(1-alpha, n-2)
      ncp = theta*sqrt(n/v1)
      power = pt(b, n-2, ncp, lower.tail = FALSE)

      delta = b/sqrt(n/v1)

      des$overallResults$overallReject = power
      des$byStageResults$rejectPerStage = power
      des$byStageResults$futilityPerStage = 1 - power
      des$byStageResults$cumulativeRejection = power
      des$byStageResults$cumulativeFutility = 1 - power
      des$byStageResults$efficacyBounds = b
      des$byStageResults$futilityBounds = b

      if (directionUpper) {
        des$byStageResults$efficacyMeanDiff = delta + meanDiffH0
        des$byStageResults$futilityMeanDiff = delta + meanDiffH0
      } else {
        des$byStageResults$efficacyMeanDiff = -delta + meanDiffH0
        des$byStageResults$futilityMeanDiff = -delta + meanDiffH0
      }
    } else {
      if (directionUpper) {
        des$byStageResults$efficacyMeanDiff =
          des$byStageResults$efficacyTheta + meanDiffH0
        des$byStageResults$futilityMeanDiff =
          des$byStageResults$futilityTheta + meanDiffH0
      } else {
        des$byStageResults$efficacyMeanDiff =
          -des$byStageResults$efficacyTheta + meanDiffH0
        des$byStageResults$futilityMeanDiff =
          -des$byStageResults$futilityTheta + meanDiffH0
      }
    }
  } else { # sample size calculation
    if (kMax == 1 && !normalApproximation) { # t-test for fixed design
      if (!any(is.na(criticalValues))) {
        alpha = 1 - pnorm(criticalValues)
      }

      n0 = (qnorm(1-alpha) + qnorm(1-beta))^2*v1/theta^2

      n = uniroot(function(n) {
        b = qt(1-alpha, n-2)
        ncp = theta*sqrt(n/v1)
        pt(b, n-2, ncp, lower.tail = FALSE) - (1 - beta)
      }, c(n0, 2*n0))$root

      if (rounding) n = ceiling(n - 1.0e-12)

      des = getDesign(
        beta = NA, IMax = n/v1, theta,
        kMax, informationRates,
        efficacyStopping, futilityStopping,
        criticalValues, alpha, typeAlphaSpending,
        parameterAlphaSpending, userAlphaSpending,
        futilityBounds, typeBetaSpending,
        parameterBetaSpending, userBetaSpending,
        spendingTime, 1)

      b = qt(1-alpha, n-2)
      ncp = theta*sqrt(n/v1)
      power = pt(b, n-2, ncp, lower.tail = FALSE)

      delta = b/sqrt(n/v1)

      des$overallResults$overallReject = power
      des$byStageResults$rejectPerStage = power
      des$byStageResults$futilityPerStage = 1 - power
      des$byStageResults$cumulativeRejection = power
      des$byStageResults$cumulativeFutility = 1 - power
      des$byStageResults$efficacyBounds = b
      des$byStageResults$futilityBounds = b
      if (directionUpper) {
        des$byStageResults$efficacyMeanDiff = delta + meanDiffH0
        des$byStageResults$futilityMeanDiff = delta + meanDiffH0
      } else {
        des$byStageResults$efficacyMeanDiff = -delta + meanDiffH0
        des$byStageResults$futilityMeanDiff = -delta + meanDiffH0
      }
    } else {
      des = getDesign(
        beta, IMax = NA, theta,
        kMax, informationRates,
        efficacyStopping, futilityStopping,
        criticalValues, alpha, typeAlphaSpending,
        parameterAlphaSpending, userAlphaSpending,
        futilityBounds, typeBetaSpending,
        parameterBetaSpending, userBetaSpending,
        spendingTime, 1)

      n = des$overallResults$information*v1

      if (rounding) {
        n = ceiling(n - 1.0e-12)
        informationRates = des$byStageResults$informationRates
        informationRates = round(n*informationRates)/n

        des = getDesign(
          beta = NA, IMax = n/v1, theta,
          kMax, informationRates,
          efficacyStopping, futilityStopping,
          criticalValues, alpha, typeAlphaSpending,
          parameterAlphaSpending, userAlphaSpending,
          futilityBounds, typeBetaSpending,
          parameterBetaSpending, userBetaSpending,
          spendingTime, 1)
      }

      if (directionUpper) {
        des$byStageResults$efficacyMeanDiff =
          des$byStageResults$efficacyTheta + meanDiffH0
        des$byStageResults$futilityMeanDiff =
          des$byStageResults$futilityTheta + meanDiffH0
      } else {
        des$byStageResults$efficacyMeanDiff =
          -des$byStageResults$efficacyTheta + meanDiffH0
        des$byStageResults$futilityMeanDiff =
          -des$byStageResults$futilityTheta + meanDiffH0
      }
    }
  }

  des$overallResults$theta = theta
  des$overallResults$numberOfSubjects = n
  des$overallResults$expectedNumberOfSubjectsH1 =
    des$overallResults$expectedInformationH1*v1
  des$overallResults$expectedNumberOfSubjectsH0 =
    des$overallResults$expectedInformationH0*v1
  des$overallResults$meanDiffH0 = meanDiffH0
  des$overallResults$meanDiff = meanDiff
  des$overallResults$stDev = stDev

  des$byStageResults$efficacyTheta = NULL
  des$byStageResults$futilityTheta = NULL
  des$byStageResults$numberOfSubjects =
    des$byStageResults$informationRates*n

  des$settings$allocationRatioPlanned = allocationRatioPlanned
  des$settings$normalApproximation = normalApproximation
  des$settings$rounding = rounding
  des$settings$varianceRatio = NULL

  attr(des, "class") = "designMeanDiff"

  des
}


#' @title Group Sequential Design for Two-Sample Mean Ratio
#' @description Obtains the power given sample size or obtains the sample
#' size given power for a group sequential design for two-sample mean
#' ratio.
#'
#' @param beta The type II error.
#' @param n The total sample size.
#' @param meanRatioH0 The mean ratio under the null hypothesis.
#'   Defaults to 1.
#' @param meanRatio The mean ratio under the alternative hypothesis.
#' @param CV The coefficient of variation. The standard deviation on the
#'   log scale is equal to \code{sqrt(log(1 + CV^2))}.
#' @param allocationRatioPlanned Allocation ratio for the active treatment
#'   versus control. Defaults to 1 for equal randomization.
#' @param normalApproximation The type of computation of the p-values.
#'   If \code{TRUE}, the variance is assumed to be known, otherwise
#'   the calculations are performed with the t distribution. The exact
#'   calculation using the t distribution is only implemented for the
#'   fixed design.
#' @param rounding Whether to round up sample size. Defaults to 1 for
#'   sample size rounding.
#' @inheritParams param_kMax
#' @param informationRates The information rates. Fixed prior to the trial.
#'   Defaults to \code{(1:kMax) / kMax} if left unspecified.
#' @inheritParams param_efficacyStopping
#' @inheritParams param_futilityStopping
#' @inheritParams param_criticalValues
#' @inheritParams param_alpha
#' @inheritParams param_typeAlphaSpending
#' @inheritParams param_parameterAlphaSpending
#' @inheritParams param_userAlphaSpending
#' @inheritParams param_futilityBounds
#' @inheritParams param_typeBetaSpending
#' @inheritParams param_parameterBetaSpending
#' @inheritParams param_userBetaSpending
#' @param spendingTime A vector of length \code{kMax} for the error spending
#'   time at each analysis. Defaults to missing, in which case, it is the
#'   same as \code{informationRates}.
#'
#' @return An S3 class \code{designMeanRatio} object with three components:
#'
#' * \code{overallResults}: A data frame containing the following variables:
#'
#'     - \code{overallReject}: The overall rejection probability.
#'
#'     - \code{alpha}: The overall significance level.
#'
#'     - \code{attainedAlpha}: The attained significance level, which is
#'       different from the overall significance level in the presence of
#'       futility stopping.
#'
#'     - \code{kMax}: The number of stages.
#'
#'     - \code{theta}: The parameter value.
#'
#'     - \code{information}: The maximum information.
#'
#'     - \code{expectedInformationH1}: The expected information under H1.
#'
#'     - \code{expectedInformationH0}: The expected information under H0.
#'
#'     - \code{drift}: The drift parameter, equal to
#'       \code{theta*sqrt(information)}.
#'
#'     - \code{inflationFactor}: The inflation factor (relative to the
#'       fixed design).
#'
#'     - \code{numberOfSubjects}: The maximum number of subjects.
#'
#'     - \code{expectedNumberOfSubjectsH1}: The expected number of subjects
#'       under H1.
#'
#'     - \code{expectedNumberOfSubjectsH0}: The expected number of subjects
#'       under H0.
#'
#'     - \code{meanRatioH0}: The mean ratio under the null hypothesis.
#'
#'     - \code{meanRatio}: The mean ratio under the alternative hypothesis.
#'
#'     - \code{CV}: The coefficient of variation.
#'
#' * \code{byStageResults}: A data frame containing the following variables:
#'
#'     - \code{informationRates}: The information rates.
#'
#'     - \code{efficacyBounds}: The efficacy boundaries on the Z-scale.
#'
#'     - \code{futilityBounds}: The futility boundaries on the Z-scale.
#'
#'     - \code{rejectPerStage}: The probability for efficacy stopping.
#'
#'     - \code{futilityPerStage}: The probability for futility stopping.
#'
#'     - \code{cumulativeRejection}: The cumulative probability for efficacy
#'       stopping.
#'
#'     - \code{cumulativeFutility}: The cumulative probability for futility
#'       stopping.
#'
#'     - \code{cumulativeAlphaSpent}: The cumulative alpha spent.
#'
#'     - \code{efficacyP}: The efficacy boundaries on the p-value scale.
#'
#'     - \code{futilityP}: The futility boundaries on the p-value scale.
#'
#'     - \code{information}: The cumulative information.
#'
#'     - \code{efficacyStopping}: Whether to allow efficacy stopping.
#'
#'     - \code{futilityStopping}: Whether to allow futility stopping.
#'
#'     - \code{rejectPerStageH0}: The probability for efficacy stopping
#'       under H0.
#'
#'     - \code{futilityPerStageH0}: The probability for futility stopping
#'       under H0.
#'
#'     - \code{cumulativeRejectionH0}: The cumulative probability for
#'       efficacy stopping under H0.
#'
#'     - \code{cumulativeFutilityH0}: The cumulative probability for futility
#'       stopping under H0.
#'
#'     - \code{numberOfSubjects}: The number of subjects.
#'
#'     - \code{efficacyMeanRatio}: The efficacy boundaries on the mean
#'       ratio scale.
#'
#'     - \code{futilityMeanRatio}: The futility boundaries on the mean
#'       ratio scale.
#'
#' * \code{settings}: A list containing the following input parameters:
#'
#'     - \code{typeAlphaSpending}: The type of alpha spending.
#'
#'     - \code{parameterAlphaSpending}: The parameter value for alpha
#'       spending.
#'
#'     - \code{userAlphaSpending}: The user defined alpha spending.
#'
#'     - \code{typeBetaSpending}: The type of beta spending.
#'
#'     - \code{parameterBetaSpending}: The parameter value for beta spending.
#'
#'     - \code{userBetaSpending}: The user defined beta spending.
#'
#'     - \code{spendingTime}: The error spending time at each analysis.
#'
#'     - \code{allocationRatioPlanned}: Allocation ratio for the active
#'       treatment versus control.
#'
#'     - \code{normalApproximation}: The type of computation of the p-values.
#'       If \code{TRUE}, the variance is assumed to be known, otherwise
#'       the calculations are performed with the t distribution.
#'
#'     - \code{rounding}: Whether to round up sample size.
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @examples
#'
#' (design1 <- getDesignMeanRatio(
#'   beta = 0.1, n = NA, meanRatio = 1.25, CV = 0.25,
#'   alpha = 0.05, normalApproximation = FALSE))
#'
#' @export
getDesignMeanRatio <- function(
    beta = NA_real_,
    n = NA_real_,
    meanRatioH0 = 1,
    meanRatio = 1.25,
    CV = 1,
    allocationRatioPlanned = 1,
    normalApproximation = TRUE,
    rounding = TRUE,
    kMax = 1L,
    informationRates = NA_real_,
    efficacyStopping = NA_integer_,
    futilityStopping = NA_integer_,
    criticalValues = NA_real_,
    alpha = 0.025,
    typeAlphaSpending = "sfOF",
    parameterAlphaSpending = NA_real_,
    userAlphaSpending = NA_real_,
    futilityBounds = NA_real_,
    typeBetaSpending = "none",
    parameterBetaSpending = NA_real_,
    userBetaSpending = NA_real_,
    spendingTime = NA_real_) {

  if (meanRatioH0 <= 0) {
    stop("meanRatioH0 must be positive")
  }

  if (meanRatio <= 0) {
    stop("meanRatio must be positive")
  }

  if (CV <= 0) {
    stop("CV must be positive")
  }

  if (allocationRatioPlanned <= 0) {
    stop("allocationRatioPlanned must be positive")
  }

  if (alpha < 0.00001 || alpha >= 1) {
    stop("alpha must lie in [0.00001, 1)")
  }


  des = getDesignMeanDiff(
    beta, n, meanDiffH0 = log(meanRatioH0),
    meanDiff = log(meanRatio),
    stDev = sqrt(log(1 + CV^2)),
    allocationRatioPlanned,
    normalApproximation, rounding,
    kMax, informationRates,
    efficacyStopping, futilityStopping,
    criticalValues, alpha, typeAlphaSpending,
    parameterAlphaSpending, userAlphaSpending,
    futilityBounds, typeBetaSpending,
    parameterBetaSpending, userBetaSpending,
    spendingTime)

  des$overallResults$meanRatioH0 = meanRatioH0
  des$overallResults$meanRatio = meanRatio
  des$overallResults$CV = CV

  des$overallResults$meanDiffH0 = NULL
  des$overallResults$meanDiff = NULL
  des$overallResults$stDev = NULL

  des$byStageResults$efficacyMeanRatio =
    exp(des$byStageResults$efficacyMeanDiff)
  des$byStageResults$futilityMeanRatio =
    exp(des$byStageResults$futilityMeanDiff)

  des$byStageResults$efficacyMeanDiff = NULL
  des$byStageResults$futilityMeanDiff = NULL

  attr(des, "class") = "designMeanRatio"

  des
}


#' @title Group Sequential Design for Mean Difference in 2x2 Crossover
#' @description Obtains the power given sample size or obtains the sample
#' size given power for a group sequential design for two-sample mean
#' difference in 2x2 crossover.
#'
#' @param beta The type II error.
#' @param n The total sample size.
#' @param meanDiffH0 The mean difference under the null hypothesis.
#'   Defaults to 0.
#' @param meanDiff The mean difference under the alternative hypothesis.
#' @param stDev The standard deviation for within-subject random error.
#' @param allocationRatioPlanned Allocation ratio for sequence A/B
#'   versus sequence B/A. Defaults to 1 for equal randomization.
#' @param normalApproximation The type of computation of the p-values.
#'   If \code{TRUE}, the variance is assumed to be known, otherwise
#'   the calculations are performed with the t distribution. The exact
#'   calculation using the t distribution is only implemented for the
#'   fixed design.
#' @param rounding Whether to round up sample size. Defaults to 1 for
#'   sample size rounding.
#' @inheritParams param_kMax
#' @param informationRates The information rates. Fixed prior to the trial.
#'   Defaults to \code{(1:kMax) / kMax} if left unspecified.
#' @inheritParams param_efficacyStopping
#' @inheritParams param_futilityStopping
#' @inheritParams param_criticalValues
#' @inheritParams param_alpha
#' @inheritParams param_typeAlphaSpending
#' @inheritParams param_parameterAlphaSpending
#' @inheritParams param_userAlphaSpending
#' @inheritParams param_futilityBounds
#' @inheritParams param_typeBetaSpending
#' @inheritParams param_parameterBetaSpending
#' @inheritParams param_userBetaSpending
#' @param spendingTime A vector of length \code{kMax} for the error spending
#'   time at each analysis. Defaults to missing, in which case, it is the
#'   same as \code{informationRates}.
#'
#' @return An S3 class \code{designMeanDiffXO} object with three components:
#'
#' * \code{overallResults}: A data frame containing the following variables:
#'
#'     - \code{overallReject}: The overall rejection probability.
#'
#'     - \code{alpha}: The overall significance level.
#'
#'     - \code{attainedAlpha}: The attained significance level, which is
#'       different from the overall significance level in the presence of
#'       futility stopping.
#'
#'     - \code{kMax}: The number of stages.
#'
#'     - \code{theta}: The parameter value.
#'
#'     - \code{information}: The maximum information.
#'
#'     - \code{expectedInformationH1}: The expected information under H1.
#'
#'     - \code{expectedInformationH0}: The expected information under H0.
#'
#'     - \code{drift}: The drift parameter, equal to
#'       \code{theta*sqrt(information)}.
#'
#'     - \code{inflationFactor}: The inflation factor (relative to the
#'       fixed design).
#'
#'     - \code{numberOfSubjects}: The maximum number of subjects.
#'
#'     - \code{expectedNumberOfSubjectsH1}: The expected number of subjects
#'       under H1.
#'
#'     - \code{expectedNumberOfSubjectsH0}: The expected number of subjects
#'       under H0.
#'
#'     - \code{meanDiffH0}: The mean difference under the null hypothesis.
#'
#'     - \code{meanDiff}: The mean difference under the alternative
#'       hypothesis.
#'
#'     - \code{stDev}: The standard deviation for within-subject random
#'       error.
#'
#' * \code{byStageResults}: A data frame containing the following variables:
#'
#'     - \code{informationRates}: The information rates.
#'
#'     - \code{efficacyBounds}: The efficacy boundaries on the Z-scale.
#'
#'     - \code{futilityBounds}: The futility boundaries on the Z-scale.
#'
#'     - \code{rejectPerStage}: The probability for efficacy stopping.
#'
#'     - \code{futilityPerStage}: The probability for futility stopping.
#'
#'     - \code{cumulativeRejection}: The cumulative probability for efficacy
#'       stopping.
#'
#'     - \code{cumulativeFutility}: The cumulative probability for futility
#'       stopping.
#'
#'     - \code{cumulativeAlphaSpent}: The cumulative alpha spent.
#'
#'     - \code{efficacyP}: The efficacy boundaries on the p-value scale.
#'
#'     - \code{futilityP}: The futility boundaries on the p-value scale.
#'
#'     - \code{information}: The cumulative information.
#'
#'     - \code{efficacyStopping}: Whether to allow efficacy stopping.
#'
#'     - \code{futilityStopping}: Whether to allow futility stopping.
#'
#'     - \code{rejectPerStageH0}: The probability for efficacy stopping
#'       under H0.
#'
#'     - \code{futilityPerStageH0}: The probability for futility stopping
#'       under H0.
#'
#'     - \code{cumulativeRejectionH0}: The cumulative probability for
#'       efficacy stopping under H0.
#'
#'     - \code{cumulativeFutilityH0}: The cumulative probability for futility
#'       stopping under H0.
#'
#'     - \code{efficacyMeanDiff}: The efficacy boundaries on the mean
#'       difference scale.
#'
#'     - \code{futilityMeanDiff}: The futility boundaries on the mean
#'       difference scale.
#'
#'     - \code{numberOfSubjects}: The number of subjects.
#'
#' * \code{settings}: A list containing the following input parameters:
#'
#'     - \code{typeAlphaSpending}: The type of alpha spending.
#'
#'     - \code{parameterAlphaSpending}: The parameter value for alpha
#'       spending.
#'
#'     - \code{userAlphaSpending}: The user defined alpha spending.
#'
#'     - \code{typeBetaSpending}: The type of beta spending.
#'
#'     - \code{parameterBetaSpending}: The parameter value for beta spending.
#'
#'     - \code{userBetaSpending}: The user defined beta spending.
#'
#'     - \code{spendingTime}: The error spending time at each analysis.
#'
#'     - \code{allocationRatioPlanned}: Allocation ratio for sequence A/B
#'       versus sequence B/A.
#'
#'     - \code{normalApproximation}: The type of computation of the p-values.
#'       If \code{TRUE}, the variance is assumed to be known, otherwise
#'       the calculations are performed with the t distribution.
#'
#'     - \code{rounding}: Whether to round up sample size.
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @examples
#'
#' (design1 <- getDesignMeanDiffXO(
#'   beta = 0.2, n = NA, meanDiff = 75, stDev = 150,
#'   normalApproximation = FALSE, alpha = 0.05))
#'
#' @export
getDesignMeanDiffXO <- function(
    beta = NA_real_,
    n = NA_real_,
    meanDiffH0 = 0,
    meanDiff = 0.5,
    stDev = 1,
    allocationRatioPlanned = 1,
    normalApproximation = TRUE,
    rounding = TRUE,
    kMax = 1L,
    informationRates = NA_real_,
    efficacyStopping = NA_integer_,
    futilityStopping = NA_integer_,
    criticalValues = NA_real_,
    alpha = 0.025,
    typeAlphaSpending = "sfOF",
    parameterAlphaSpending = NA_real_,
    userAlphaSpending = NA_real_,
    futilityBounds = NA_real_,
    typeBetaSpending = "none",
    parameterBetaSpending = NA_real_,
    userBetaSpending = NA_real_,
    spendingTime = NA_real_) {

  if (is.na(beta) && is.na(n)) {
    stop("beta and n cannot be both missing")
  }

  if (!is.na(beta) && !is.na(n)) {
    stop("Only one of beta and n should be provided")
  }

  if (!is.na(beta) && (beta >= 1-alpha || beta < 0.0001)) {
    stop("beta must lie in [0.0001, 1-alpha)")
  }

  if (!is.na(n) && n <= 0) {
    stop("n must be positive")
  }

  if (stDev <= 0) {
    stop("stDev must be positive")
  }

  if (allocationRatioPlanned <= 0) {
    stop("allocationRatioPlanned must be positive")
  }

  if (alpha < 0.00001 || alpha >= 1) {
    stop("alpha must lie in [0.00001, 1)")
  }

  if (any(is.na(informationRates))) {
    informationRates = (1:kMax)/kMax
  }


  r = allocationRatioPlanned/(1 + allocationRatioPlanned)

  directionUpper = meanDiff > meanDiffH0

  theta = ifelse(directionUpper, meanDiff - meanDiffH0,
                 meanDiffH0 - meanDiff)

  # variance for one sampling unit
  v1 = stDev^2/(2*r*(1-r))

  if (is.na(beta)) { # power calculation
    if (rounding) {
      n = ceiling(n - 1.0e-12)
      informationRates = round(n*informationRates)/n
    }

    des = getDesign(
      beta = NA, IMax = n/v1, theta,
      kMax, informationRates,
      efficacyStopping, futilityStopping,
      criticalValues, alpha, typeAlphaSpending,
      parameterAlphaSpending, userAlphaSpending,
      futilityBounds, typeBetaSpending,
      parameterBetaSpending, userBetaSpending,
      spendingTime, 1)

    if (kMax == 1 && !normalApproximation) { # t-test for fixed design
      if (!any(is.na(criticalValues))) {
        alpha = 1 - pnorm(criticalValues)
      }

      b = qt(1-alpha, n-2)
      ncp = theta*sqrt(n/v1)
      power = pt(b, n-2, ncp, lower.tail = FALSE)

      delta = b/sqrt(n/v1)

      des$overallResults$overallReject = power
      des$byStageResults$rejectPerStage = power
      des$byStageResults$futilityPerStage = 1 - power
      des$byStageResults$cumulativeRejection = power
      des$byStageResults$cumulativeFutility = 1 - power
      des$byStageResults$efficacyBounds = b
      des$byStageResults$futilityBounds = b

      if (directionUpper) {
        des$byStageResults$efficacyMeanDiff = delta + meanDiffH0
        des$byStageResults$futilityMeanDiff = delta + meanDiffH0
      } else {
        des$byStageResults$efficacyMeanDiff = -delta + meanDiffH0
        des$byStageResults$futilityMeanDiff = -delta + meanDiffH0
      }
    } else {
      if (directionUpper) {
        des$byStageResults$efficacyMeanDiff =
          des$byStageResults$efficacyTheta + meanDiffH0
        des$byStageResults$futilityMeanDiff =
          des$byStageResults$futilityTheta + meanDiffH0
      } else {
        des$byStageResults$efficacyMeanDiff =
          -des$byStageResults$efficacyTheta + meanDiffH0
        des$byStageResults$futilityMeanDiff =
          -des$byStageResults$futilityTheta + meanDiffH0
      }
    }
  } else { # sample size calculation
    if (kMax == 1 && !normalApproximation) { # t-test for fixed design
      if (!any(is.na(criticalValues))) {
        alpha = 1 - pnorm(criticalValues)
      }

      n0 = (qnorm(1-alpha) + qnorm(1-beta))^2*v1/theta^2

      n = uniroot(function(n) {
        b = qt(1-alpha, n-2)
        ncp = theta*sqrt(n/v1)
        pt(b, n-2, ncp, lower.tail = FALSE) - (1 - beta)
      }, c(n0, 2*n0))$root

      if (rounding) n = ceiling(n - 1.0e-12)

      des = getDesign(
        beta = NA, IMax = n/v1, theta,
        kMax, informationRates,
        efficacyStopping, futilityStopping,
        criticalValues, alpha, typeAlphaSpending,
        parameterAlphaSpending, userAlphaSpending,
        futilityBounds, typeBetaSpending,
        parameterBetaSpending, userBetaSpending,
        spendingTime, 1)

      b = qt(1-alpha, n-2)
      ncp = theta*sqrt(n/v1)
      power = pt(b, n-2, ncp, lower.tail = FALSE)

      delta = b/sqrt(n/v1)

      des$overallResults$overallReject = power
      des$byStageResults$rejectPerStage = power
      des$byStageResults$futilityPerStage = 1 - power
      des$byStageResults$cumulativeRejection = power
      des$byStageResults$cumulativeFutility = 1 - power
      des$byStageResults$efficacyBounds = b
      des$byStageResults$futilityBounds = b
      if (directionUpper) {
        des$byStageResults$efficacyMeanDiff = delta + meanDiffH0
        des$byStageResults$futilityMeanDiff = delta + meanDiffH0
      } else {
        des$byStageResults$efficacyMeanDiff = -delta + meanDiffH0
        des$byStageResults$futilityMeanDiff = -delta + meanDiffH0
      }
    } else {
      des = getDesign(
        beta, IMax = NA, theta,
        kMax, informationRates,
        efficacyStopping, futilityStopping,
        criticalValues, alpha, typeAlphaSpending,
        parameterAlphaSpending, userAlphaSpending,
        futilityBounds, typeBetaSpending,
        parameterBetaSpending, userBetaSpending,
        spendingTime, 1)

      n = des$overallResults$information*v1

      if (rounding) {
        n = ceiling(n - 1.0e-12)
        informationRates = des$byStageResults$informationRates
        informationRates = round(n*informationRates)/n

        des = getDesign(
          beta = NA, IMax = n/v1, theta,
          kMax, informationRates,
          efficacyStopping, futilityStopping,
          criticalValues, alpha, typeAlphaSpending,
          parameterAlphaSpending, userAlphaSpending,
          futilityBounds, typeBetaSpending,
          parameterBetaSpending, userBetaSpending,
          spendingTime, 1)
      }

      if (directionUpper) {
        des$byStageResults$efficacyMeanDiff =
          des$byStageResults$efficacyTheta + meanDiffH0
        des$byStageResults$futilityMeanDiff =
          des$byStageResults$futilityTheta + meanDiffH0
      } else {
        des$byStageResults$efficacyMeanDiff =
          -des$byStageResults$efficacyTheta + meanDiffH0
        des$byStageResults$futilityMeanDiff =
          -des$byStageResults$futilityTheta + meanDiffH0
      }
    }
  }

  des$overallResults$theta = theta
  des$overallResults$numberOfSubjects = n
  des$overallResults$expectedNumberOfSubjectsH1 =
    des$overallResults$expectedInformationH1*v1
  des$overallResults$expectedNumberOfSubjectsH0 =
    des$overallResults$expectedInformationH0*v1
  des$overallResults$meanDiffH0 = meanDiffH0
  des$overallResults$meanDiff = meanDiff
  des$overallResults$stDev = stDev

  des$byStageResults$efficacyTheta = NULL
  des$byStageResults$futilityTheta = NULL
  des$byStageResults$numberOfSubjects =
    des$byStageResults$informationRates*n

  des$settings$allocationRatioPlanned = allocationRatioPlanned
  des$settings$normalApproximation = normalApproximation
  des$settings$rounding = rounding
  des$settings$varianceRatio = NULL

  attr(des, "class") = "designMeanDiffXO"

  des
}


#' @title Group Sequential Design for Mean Ratio in 2x2 Crossover
#' @description Obtains the power given sample size or obtains the sample
#' size given power for a group sequential design for two-sample mean
#' ratio in 2x2 crossover.
#'
#' @param beta The type II error.
#' @param n The total sample size.
#' @param meanRatioH0 The mean ratio under the null hypothesis.
#'   Defaults to 1.
#' @param meanRatio The mean ratio under the alternative hypothesis.
#' @param CV The coefficient of variation. The standard deviation on the
#'   log scale is equal to \code{sqrt(log(1 + CV^2))}.
#' @param allocationRatioPlanned Allocation ratio for sequence A/B
#'   versus sequence B/A. Defaults to 1 for equal randomization.
#' @param normalApproximation The type of computation of the p-values.
#'   If \code{TRUE}, the variance is assumed to be known, otherwise
#'   the calculations are performed with the t distribution. The exact
#'   calculation using the t distribution is only implemented for the
#'   fixed design.
#' @param rounding Whether to round up sample size. Defaults to 1 for
#'   sample size rounding.
#' @inheritParams param_kMax
#' @param informationRates The information rates. Fixed prior to the trial.
#'   Defaults to \code{(1:kMax) / kMax} if left unspecified.
#' @inheritParams param_efficacyStopping
#' @inheritParams param_futilityStopping
#' @inheritParams param_criticalValues
#' @inheritParams param_alpha
#' @inheritParams param_typeAlphaSpending
#' @inheritParams param_parameterAlphaSpending
#' @inheritParams param_userAlphaSpending
#' @inheritParams param_futilityBounds
#' @inheritParams param_typeBetaSpending
#' @inheritParams param_parameterBetaSpending
#' @inheritParams param_userBetaSpending
#' @param spendingTime A vector of length \code{kMax} for the error spending
#'   time at each analysis. Defaults to missing, in which case, it is the
#'   same as \code{informationRates}.
#'
#' @return An S3 class \code{designMeanRatioXO} object with three components:
#'
#' * \code{overallResults}: A data frame containing the following variables:
#'
#'     - \code{overallReject}: The overall rejection probability.
#'
#'     - \code{alpha}: The overall significance level.
#'
#'     - \code{attainedAlpha}: The attained significance level, which is
#'       different from the overall significance level in the presence of
#'       futility stopping.
#'
#'     - \code{kMax}: The number of stages.
#'
#'     - \code{theta}: The parameter value.
#'
#'     - \code{information}: The maximum information.
#'
#'     - \code{expectedInformationH1}: The expected information under H1.
#'
#'     - \code{expectedInformationH0}: The expected information under H0.
#'
#'     - \code{drift}: The drift parameter, equal to
#'       \code{theta*sqrt(information)}.
#'
#'     - \code{inflationFactor}: The inflation factor (relative to the
#'       fixed design).
#'
#'     - \code{numberOfSubjects}: The maximum number of subjects.
#'
#'     - \code{expectedNumberOfSubjectsH1}: The expected number of subjects
#'       under H1.
#'
#'     - \code{expectedNumberOfSubjectsH0}: The expected number of subjects
#'       under H0.
#'
#'     - \code{meanRatioH0}: The mean ratio under the null hypothesis.
#'
#'     - \code{meanRatio}: The mean ratio under the alternative hypothesis.
#'
#'     - \code{CV}: The coefficient of variation.
#'
#' * \code{byStageResults}: A data frame containing the following variables:
#'
#'     - \code{informationRates}: The information rates.
#'
#'     - \code{efficacyBounds}: The efficacy boundaries on the Z-scale.
#'
#'     - \code{futilityBounds}: The futility boundaries on the Z-scale.
#'
#'     - \code{rejectPerStage}: The probability for efficacy stopping.
#'
#'     - \code{futilityPerStage}: The probability for futility stopping.
#'
#'     - \code{cumulativeRejection}: The cumulative probability for efficacy
#'       stopping.
#'
#'     - \code{cumulativeFutility}: The cumulative probability for futility
#'       stopping.
#'
#'     - \code{cumulativeAlphaSpent}: The cumulative alpha spent.
#'
#'     - \code{efficacyMeanRatio}: The efficacy boundaries on the mean
#'       ratio scale.
#'
#'     - \code{futilityMeanRatio}: The futility boundaries on the mean
#'       ratio scale.
#'
#'     - \code{efficacyP}: The efficacy boundaries on the p-value scale.
#'
#'     - \code{futilityP}: The futility boundaries on the p-value scale.
#'
#'     - \code{information}: The cumulative information.
#'
#'     - \code{efficacyStopping}: Whether to allow efficacy stopping.
#'
#'     - \code{futilityStopping}: Whether to allow futility stopping.
#'
#'     - \code{rejectPerStageH0}: The probability for efficacy stopping
#'       under H0.
#'
#'     - \code{futilityPerStageH0}: The probability for futility stopping
#'       under H0.
#'
#'     - \code{cumulativeRejectionH0}: The cumulative probability for
#'       efficacy stopping under H0.
#'
#'     - \code{cumulativeFutilityH0}: The cumulative probability for futility
#'       stopping under H0.
#'
#'     - \code{numberOfSubjects}: The number of subjects.
#'
#' * \code{settings}: A list containing the following input parameters:
#'
#'     - \code{typeAlphaSpending}: The type of alpha spending.
#'
#'     - \code{parameterAlphaSpending}: The parameter value for alpha
#'       spending.
#'
#'     - \code{userAlphaSpending}: The user defined alpha spending.
#'
#'     - \code{typeBetaSpending}: The type of beta spending.
#'
#'     - \code{parameterBetaSpending}: The parameter value for beta spending.
#'
#'     - \code{userBetaSpending}: The user defined beta spending.
#'
#'     - \code{spendingTime}: The error spending time at each analysis.
#'
#'     - \code{allocationRatioPlanned}: Allocation ratio for sequence A/B
#'       versus sequence B/A.
#'
#'     - \code{normalApproximation}: The type of computation of the p-values.
#'       If \code{TRUE}, the variance is assumed to be known, otherwise
#'       the calculations are performed with the t distribution.
#'
#'     - \code{rounding}: Whether to round up sample size.
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @examples
#'
#' (design1 <- getDesignMeanRatioXO(
#'   beta = 0.1, n = NA, meanRatio = 1.25, CV = 0.25,
#'   alpha = 0.05, normalApproximation = FALSE))
#'
#' @export
getDesignMeanRatioXO <- function(
    beta = NA_real_,
    n = NA_real_,
    meanRatioH0 = 1,
    meanRatio = 1.25,
    CV = 1,
    allocationRatioPlanned = 1,
    normalApproximation = TRUE,
    rounding = TRUE,
    kMax = 1L,
    informationRates = NA_real_,
    efficacyStopping = NA_integer_,
    futilityStopping = NA_integer_,
    criticalValues = NA_real_,
    alpha = 0.025,
    typeAlphaSpending = "sfOF",
    parameterAlphaSpending = NA_real_,
    userAlphaSpending = NA_real_,
    futilityBounds = NA_real_,
    typeBetaSpending = "none",
    parameterBetaSpending = NA_real_,
    userBetaSpending = NA_real_,
    spendingTime = NA_real_) {

  if (meanRatioH0 <= 0) {
    stop("meanRatioH0 must be positive")
  }

  if (meanRatio <= 0) {
    stop("meanRatio must be positive")
  }

  if (CV <= 0) {
    stop("CV must be positive")
  }

  if (allocationRatioPlanned <= 0) {
    stop("allocationRatioPlanned must be positive")
  }

  if (alpha < 0.00001 || alpha >= 1) {
    stop("alpha must lie in [0.00001, 1)")
  }


  des = getDesignMeanDiffXO(
    beta, n, meanDiffH0 = log(meanRatioH0),
    meanDiff = log(meanRatio),
    stDev = sqrt(log(1 + CV^2)),
    allocationRatioPlanned,
    normalApproximation, rounding,
    kMax, informationRates,
    efficacyStopping, futilityStopping,
    criticalValues, alpha, typeAlphaSpending,
    parameterAlphaSpending, userAlphaSpending,
    futilityBounds, typeBetaSpending,
    parameterBetaSpending, userBetaSpending,
    spendingTime)

  des$overallResults$meanRatioH0 = meanRatioH0
  des$overallResults$meanRatio = meanRatio
  des$overallResults$CV = CV

  des$overallResults$meanDiffH0 = NULL
  des$overallResults$meanDiff = NULL
  des$overallResults$stDev = NULL

  des$byStageResults$efficacyMeanRatio =
    exp(des$byStageResults$efficacyMeanDiff)
  des$byStageResults$futilityMeanRatio =
    exp(des$byStageResults$futilityMeanDiff)
  des$byStageResults$efficacyMeanDiff = NULL
  des$byStageResults$futilityMeanDiff = NULL

  attr(des, "class") = "designMeanRatioXO"

  des
}


#' @title Group Sequential Design for Equivalence in Paired Mean
#' Difference
#' @description Obtains the power given sample size or obtains the sample
#' size given power for a group sequential design for equivalence in
#' paired mean difference.
#'
#' @param beta The type II error.
#' @param n The total sample size.
#' @param pairedDiffLower The lower equivalence limit of paired difference.
#' @param pairedDiffUpper The upper equivalence limit of paired difference.
#' @param pairedDiff The paired difference under the alternative
#'   hypothesis.
#' @param stDev The standard deviation for paired difference.
#' @param normalApproximation The type of computation of the p-values.
#'   If \code{TRUE}, the variance is assumed to be known, otherwise
#'   the calculations are performed with the t distribution. The exact
#'   calculation using the t distribution is only implemented for the
#'   fixed design.
#' @param rounding Whether to round up sample size. Defaults to 1 for
#'   sample size rounding.
#' @inheritParams param_kMax
#' @param informationRates The information rates. Fixed prior to the trial.
#'   Defaults to \code{(1:kMax) / kMax} if left unspecified.
#' @inheritParams param_criticalValues
#' @param alpha The significance level for each of the two one-sided
#'   tests. Defaults to 0.05.
#' @inheritParams param_typeAlphaSpending
#' @inheritParams param_parameterAlphaSpending
#' @inheritParams param_userAlphaSpending
#' @param spendingTime A vector of length \code{kMax} for the error spending
#'   time at each analysis. Defaults to missing, in which case, it is the
#'   same as \code{informationRates}.
#'
#' @return An S3 class \code{designPairedMeanDiffEquiv} object with three
#' components:
#'
#' * \code{overallResults}: A data frame containing the following variables:
#'
#'     - \code{overallReject}: The overall rejection probability.
#'
#'     - \code{alpha}: The significance level for each of the two one-sided
#'       tests. Defaults to 0.05.
#'
#'     - \code{attainedAlpha}: The attained significance level under H0.
#'
#'     - \code{kMax}: The number of stages.
#'
#'     - \code{information}: The maximum information.
#'
#'     - \code{expectedInformationH1}: The expected information under H1.
#'
#'     - \code{expectedInformationH0}: The expected information under H0.
#'
#'     - \code{numberOfSubjects}: The maximum number of subjects.
#'
#'     - \code{expectedNumberOfSubjectsH1}: The expected number of subjects
#'       under H1.
#'
#'     - \code{expectedNumberOfSubjectsH0}: The expected number of subjects
#'       under H0.
#'
#'     - \code{pairedDiffLower}: The lower equivalence limit of paired
#'       difference.
#'
#'     - \code{pairedDiffUpper}: The upper equivalence limit of paired
#'       difference.
#'
#'     - \code{pairedDiff}: The paired difference under the alternative
#'       hypothesis.
#'
#'     - \code{stDev}: The standard deviation for paired difference.
#'
#' * \code{byStageResults}: A data frame containing the following variables:
#'
#'     - \code{informationRates}: The information rates.
#'
#'     - \code{efficacyBounds}: The efficacy boundaries on the Z-scale for
#'       each of the two one-sided tests.
#'
#'     - \code{rejectPerStage}: The probability for efficacy stopping.
#'
#'     - \code{cumulativeRejection}: The cumulative probability for efficacy
#'       stopping.
#'
#'     - \code{cumulativeAlphaSpent}: The cumulative alpha for each of
#'       the two one-sided tests.
#'
#'     - \code{cumulativeAttainedAlpha}: The cumulative probability for
#'       efficacy stopping under H0.
#'
#'     - \code{efficacyPairedDiffLower}: The efficacy boundaries on the
#'       paired difference scale for the one-sided null hypothesis on the
#'       lower equivalence limit.
#'
#'     - \code{efficacyPairedDiffUpper}: The efficacy boundaries on the
#'       paired difference scale for the one-sided null hypothesis on the
#'       upper equivalence limit.
#'
#'     - \code{efficacyP}: The efficacy bounds on the p-value scale for
#'       each of the two one-sided tests.
#'
#'     - \code{information}: The cumulative information.
#'
#'     - \code{numberOfSubjects}: The number of subjects.
#'
#' * \code{settings}: A list containing the following input parameters:
#'
#'     - \code{typeAlphaSpending}: The type of alpha spending.
#'
#'     - \code{parameterAlphaSpending}: The parameter value for alpha
#'       spending.
#'
#'     - \code{userAlphaSpending}: The user defined alpha spending.
#'
#'     - \code{spendingTime}: The error spending time at each analysis.
#'
#'     - \code{normalApproximation}: The type of computation of the p-values.
#'       If \code{TRUE}, the variance is assumed to be known, otherwise
#'       the calculations are performed with the t distribution. The exact
#'       calculation using the t distribution is only implemented for the
#'       fixed design.
#'
#'     - \code{rounding}: Whether to round up sample size.
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @examples
#'
#' # Example 1: group sequential trial power calculation
#' (design1 <- getDesignPairedMeanDiffEquiv(
#'   beta = 0.1, n = NA, pairedDiffLower = -1.3, pairedDiffUpper = 1.3,
#'   pairedDiff = 0, stDev = 2.2,
#'   kMax = 4, alpha = 0.05, typeAlphaSpending = "sfOF"))
#'
#' # Example 2: sample size calculation for t-test
#' (design2 <- getDesignPairedMeanDiffEquiv(
#'   beta = 0.1, n = NA, pairedDiffLower = -1.3, pairedDiffUpper = 1.3,
#'   pairedDiff = 0, stDev = 2.2,
#'   normalApproximation = FALSE, alpha = 0.05))
#'
#' @export
getDesignPairedMeanDiffEquiv <- function(
    beta = NA_real_,
    n = NA_real_,
    pairedDiffLower = NA_real_,
    pairedDiffUpper = NA_real_,
    pairedDiff = 0,
    stDev = 1,
    normalApproximation = TRUE,
    rounding = TRUE,
    kMax = 1L,
    informationRates = NA_real_,
    alpha = 0.05,
    typeAlphaSpending = "sfOF",
    parameterAlphaSpending = NA_real_,
    userAlphaSpending = NA_real_,
    spendingTime = NA_real_) {

  if (is.na(beta) && is.na(n)) {
    stop("beta and n cannot be both missing")
  }

  if (!is.na(beta) && !is.na(n)) {
    stop("Only one of beta and n should be provided")
  }

  if (!is.na(beta) && (beta >= 1-alpha || beta < 0.0001)) {
    stop("beta must lie in [0.0001, 1-alpha)")
  }

  if (!is.na(n) && n <= 0) {
    stop("n must be positive")
  }

  if (is.na(pairedDiffLower)) {
    stop("pairedDiffLower must be provided")
  }

  if (is.na(pairedDiffUpper)) {
    stop("pairedDiffUpper must be provided")
  }

  if (pairedDiffLower >= pairedDiff) {
    stop("pairedDiffLower must be less than pairedDiff")
  }

  if (pairedDiffUpper <= pairedDiff) {
    stop("pairedDiffUpper must be greater than pairedDiff")
  }

  if (stDev <= 0) {
    stop("stDev must be positive")
  }

  if (alpha < 0.00001 || alpha >= 1) {
    stop("alpha must lie in [0.00001, 1)")
  }

  if (any(is.na(informationRates))) {
    informationRates = (1:kMax)/kMax
  }


  # variance for one sampling unit
  v1 = stDev^2

  f <- function(n) { # power for two one-sided t-tests
    b = qt(1-alpha, n-1)
    ncpLower = (pairedDiff - pairedDiffLower)*sqrt(n/v1)
    powerLower = pt(b, n-1, ncpLower, lower.tail = FALSE)
    ncpUpper = (pairedDiffUpper - pairedDiff)*sqrt(n/v1)
    powerUpper = pt(b, n-1, ncpUpper, lower.tail = FALSE)
    power = powerLower + powerUpper - 1
    power
  }

  if (is.na(n)) { # calculate sample size
    des = getDesignEquiv(
      beta = beta, IMax = NA, thetaLower = pairedDiffLower,
      thetaUpper = pairedDiffUpper, theta = pairedDiff,
      kMax = kMax, informationRates = informationRates,
      alpha = alpha, typeAlphaSpending = typeAlphaSpending,
      parameterAlphaSpending = parameterAlphaSpending,
      userAlphaSpending = userAlphaSpending,
      spendingTime = spendingTime)

    n = des$overallResults$information*v1

    if (kMax == 1 && !normalApproximation) { # t-test for fixed design
      n = uniroot(function(n) f(n) - (1-beta), c(0.5*n, 1.5*n))$root
    }
  }

  if (rounding) {
    n = ceiling(n - 1.0e-12)
    informationRates = round(n*informationRates)/n
  }

  if (is.na(beta) || rounding) { # calculate power
    des = getDesignEquiv(
      beta = NA, IMax = n/v1, thetaLower = pairedDiffLower,
      thetaUpper = pairedDiffUpper, theta = pairedDiff,
      kMax = kMax, informationRates = informationRates,
      alpha = alpha, typeAlphaSpending = typeAlphaSpending,
      parameterAlphaSpending = parameterAlphaSpending,
      userAlphaSpending = userAlphaSpending,
      spendingTime = spendingTime)
  }

  if (kMax == 1 && !normalApproximation) { # t-test for fixed design
    power = f(n)

    b = qt(1-alpha, n-1)
    ncp = (pairedDiffUpper - pairedDiffLower)*sqrt(n/v1)
    attainedAlpha = alpha - pt(b, n-1, ncp)

    des$overallResults$overallReject = power
    des$overallResults$attainedAlpha = attainedAlpha
    des$overallResults$information = n/v1
    des$overallResults$expectedInformationH1 = n/v1
    des$overallResults$expectedInformationH0 = n/v1

    des$byStageResults$efficacyBounds = b
    des$byStageResults$rejectPerStage = power
    des$byStageResults$cumulativeRejection = power
    des$byStageResults$cumulativeAttainedAlpha = attainedAlpha
    des$byStageResults$efficacyPairedDiffLower = b*sqrt(v1/n) +
      pairedDiffLower
    des$byStageResults$efficacyPairedDiffUpper = -b*sqrt(v1/n) +
      pairedDiffUpper
    des$byStageResults$information = n/v1
  } else {
    des$overallResults$attainedAlpha =
      des$overallResults$attainedAlphaH10
    des$overallResults$expectedInformationH0 =
      des$overallResults$expectedInformationH10

    des$byStageResults$cumulativeAttainedAlpha =
      des$byStageResults$cumulativeAttainedAlphaH10
    des$byStageResults$efficacyPairedDiffLower =
      des$byStageResults$efficacyThetaLower
    des$byStageResults$efficacyPairedDiffUpper =
      des$byStageResults$efficacyThetaUpper
  }

  des$overallResults$numberOfSubjects = n
  des$overallResults$expectedNumberOfSubjectsH1 =
    des$overallResults$expectedInformationH1*v1
  des$overallResults$expectedNumberOfSubjectsH0 =
    des$overallResults$expectedInformationH0*v1
  des$overallResults$pairedDiffLower = pairedDiffLower
  des$overallResults$pairedDiffUpper = pairedDiffUpper
  des$overallResults$pairedDiff = pairedDiff
  des$overallResults$stDev = stDev
  des$overallResults <-
    des$overallResults[, c("overallReject", "alpha", "attainedAlpha",
                           "kMax", "information", "expectedInformationH1",
                           "expectedInformationH0", "numberOfSubjects",
                           "expectedNumberOfSubjectsH1",
                           "expectedNumberOfSubjectsH0", "pairedDiffLower",
                           "pairedDiffUpper", "pairedDiff", "stDev")]

  des$byStageResults$numberOfSubjects = n*informationRates
  des$byStageResults <-
    des$byStageResults[, c("informationRates", "efficacyBounds",
                           "rejectPerStage", "cumulativeRejection",
                           "cumulativeAlphaSpent", "cumulativeAttainedAlpha",
                           "efficacyPairedDiffLower",
                           "efficacyPairedDiffUpper", "efficacyP",
                           "information", "numberOfSubjects")]

  des$settings$normalApproximation = normalApproximation
  des$settings$rounding = rounding
  des$settings <-
    des$settings[c("typeAlphaSpending", "parameterAlphaSpending",
                   "userAlphaSpending", "spendingTime",
                   "normalApproximation", "rounding")]

  attr(des, "class") = "designPairedMeanDiffEquiv"

  des
}


#' @title Group Sequential Design for Equivalence in Paired Mean Ratio
#' @description Obtains the power given sample size or obtains the sample
#' size given power for a group sequential design for equivalence in
#' paired mean ratio.
#'
#' @param beta The type II error.
#' @param n The total sample size.
#' @param pairedRatioLower The lower equivalence limit of paired ratio.
#' @param pairedRatioUpper The upper equivalence limit of paired ratio.
#' @param pairedRatio The paired ratio under the alternative
#'   hypothesis.
#' @param CV The coefficient of variation for paired ratio.
#' @param normalApproximation The type of computation of the p-values.
#'   If \code{TRUE}, the variance is assumed to be known, otherwise
#'   the calculations are performed with the t distribution. The exact
#'   calculation using the t distribution is only implemented for the
#'   fixed design.
#' @param rounding Whether to round up sample size. Defaults to 1 for
#'   sample size rounding.
#' @inheritParams param_kMax
#' @param informationRates The information rates. Fixed prior to the trial.
#'   Defaults to \code{(1:kMax) / kMax} if left unspecified.
#' @inheritParams param_criticalValues
#' @param alpha The significance level for each of the two one-sided
#'   tests. Defaults to 0.05.
#' @inheritParams param_typeAlphaSpending
#' @inheritParams param_parameterAlphaSpending
#' @inheritParams param_userAlphaSpending
#' @param spendingTime A vector of length \code{kMax} for the error spending
#'   time at each analysis. Defaults to missing, in which case, it is the
#'   same as \code{informationRates}.
#'
#' @return An S3 class \code{designPairedMeanRatioEquiv} object with three
#' components:
#'
#' * \code{overallResults}: A data frame containing the following variables:
#'
#'     - \code{overallReject}: The overall rejection probability.
#'
#'     - \code{alpha}: The significance level for each of the two one-sided
#'       tests. Defaults to 0.05.
#'
#'     - \code{attainedAlpha}: The attained significance level under H0.
#'
#'     - \code{kMax}: The number of stages.
#'
#'     - \code{information}: The maximum information.
#'
#'     - \code{expectedInformationH1}: The expected information under H1.
#'
#'     - \code{expectedInformationH0}: The expected information under H0.
#'
#'     - \code{numberOfSubjects}: The maximum number of subjects.
#'
#'     - \code{expectedNumberOfSubjectsH1}: The expected number of subjects
#'       under H1.
#'
#'     - \code{expectedNumberOfSubjectsH0}: The expected number of subjects
#'       under H0.
#'
#'     - \code{pairedRatioLower}: The lower equivalence limit of paired
#'       ratio.
#'
#'     - \code{pairedRatioUpper}: The upper equivalence limit of paired
#'       ratio.
#'
#'     - \code{pairedRatio}: The paired ratio under the alternative
#'       hypothesis.
#'
#'     - \code{CV}: The coefficient of variation for paired ratios.
#'
#' * \code{byStageResults}: A data frame containing the following variables:
#'
#'     - \code{informationRates}: The information rates.
#'
#'     - \code{efficacyBounds}: The efficacy boundaries on the Z-scale for
#'       each of the two one-sided tests.
#'
#'     - \code{rejectPerStage}: The probability for efficacy stopping.
#'
#'     - \code{cumulativeRejection}: The cumulative probability for efficacy
#'       stopping.
#'
#'     - \code{cumulativeAlphaSpent}: The cumulative alpha for each of
#'       the two one-sided tests.
#'
#'     - \code{cumulativeAttainedAlpha}: The cumulative alpha attained under
#'       H0.
#'
#'     - \code{efficacyP}: The efficacy bounds on the p-value scale for
#'       each of the two one-sided tests.
#'
#'     - \code{information}:  The cumulative information.
#'
#'     - \code{numberOfSubjects}: The number of subjects.
#'
#'     - \code{efficacyPairedRatioLower}: The efficacy boundaries on the
#'       paired ratio scale for the one-sided null hypothesis on the
#'       lower equivalence limit.
#'
#'     - \code{efficacyPairedRatioUpper}: The efficacy boundaries on the
#'       paired ratio scale for the one-sided null hypothesis on the
#'       upper equivalence limit.
#'
#' * \code{settings}: A list containing the following input parameters:
#'
#'     - \code{typeAlphaSpending}: The type of alpha spending.
#'
#'     - \code{parameterAlphaSpending}: The parameter value for alpha
#'       spending.
#'
#'     - \code{userAlphaSpending}: The user defined alpha spending.
#'
#'     - \code{spendingTime}: The error spending time at each analysis.
#'
#'     - \code{normalApproximation}: The type of computation of the p-values.
#'       If \code{TRUE}, the variance is assumed to be known, otherwise
#'       the calculations are performed with the t distribution. The exact
#'       calculation using the t distribution is only implemented for the
#'       fixed design.
#'
#'     - \code{rounding}: Whether to round up sample size.
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @examples
#'
#' # Example 1: group sequential trial power calculation
#' (design1 <- getDesignPairedMeanRatioEquiv(
#'   beta = 0.1, n = NA, pairedRatioLower = 0.8, pairedRatioUpper = 1.25,
#'   pairedRatio = 1, CV = 0.35,
#'   kMax = 4, alpha = 0.05, typeAlphaSpending = "sfOF"))
#'
#' # Example 2: sample size calculation for t-test
#' (design2 <- getDesignPairedMeanRatioEquiv(
#'   beta = 0.1, n = NA, pairedRatioLower = 0.8, pairedRatioUpper = 1.25,
#'   pairedRatio = 1, CV = 0.35,
#'   normalApproximation = FALSE, alpha = 0.05))
#'
#' @export
getDesignPairedMeanRatioEquiv <- function(
    beta = NA_real_,
    n = NA_real_,
    pairedRatioLower = NA_real_,
    pairedRatioUpper = NA_real_,
    pairedRatio = 1,
    CV = 1,
    normalApproximation = TRUE,
    rounding = TRUE,
    kMax = 1L,
    informationRates = NA_real_,
    alpha = 0.05,
    typeAlphaSpending = "sfOF",
    parameterAlphaSpending = NA_real_,
    userAlphaSpending = NA_real_,
    spendingTime = NA_real_) {

  if (is.na(pairedRatioLower)) {
    stop("pairedRatioLower must be provided")
  }

  if (is.na(pairedRatioUpper)) {
    stop("pairedRatioUpper must be provided")
  }

  if (pairedRatioLower <= 0) {
    stop("pairedRatioLower must be positive")
  }

  if (pairedRatioLower >= pairedRatio) {
    stop("pairedRatioLower must be less than pairedRatio")
  }

  if (pairedRatioUpper <= pairedRatio) {
    stop("pairedRatioUpper must be greater than pairedRatio")
  }

  if (CV <= 0) {
    stop("CV must be positive")
  }

  if (alpha < 0.00001 || alpha >= 1) {
    stop("alpha must lie in [0.00001, 1)")
  }


  des = getDesignPairedMeanDiffEquiv(
    beta, n, pairedDiffLower = log(pairedRatioLower),
    pairedDiffUpper = log(pairedRatioUpper),
    pairedDiff = log(pairedRatio),
    stDev = sqrt(log(1 + CV^2)),
    normalApproximation, rounding,
    kMax, informationRates,
    alpha, typeAlphaSpending,
    parameterAlphaSpending, userAlphaSpending,
    spendingTime)

  des$overallResults$pairedRatioLower = pairedRatioLower
  des$overallResults$pairedRatioUpper = pairedRatioUpper
  des$overallResults$pairedRatio = pairedRatio
  des$overallResults$CV = CV

  des$overallResults$pairedDiffLower = NULL
  des$overallResults$pairedDiffUpper = NULL
  des$overallResults$pairedDiff = NULL
  des$overallResults$stDev = NULL

  des$byStageResults$efficacyPairedRatioLower =
    exp(des$byStageResults$efficacyPairedDiffLower)
  des$byStageResults$efficacyPairedRatioUpper =
    exp(des$byStageResults$efficacyPairedDiffUpper)

  des$byStageResults$efficacyPairedDiffLower = NULL
  des$byStageResults$efficacyPairedDiffUpper = NULL

  attr(des, "class") = "designPairedMeanRatioEquiv"

  des
}


#' @title Group Sequential Design for Equivalence in Two-Sample
#' Mean Difference
#' @description Obtains the power given sample size or obtains the sample
#' size given power for a group sequential design for equivalence in
#' two-sample mean difference.
#'
#' @param beta The type II error.
#' @param n The total sample size.
#' @param meanDiffLower The lower equivalence limit of mean difference.
#' @param meanDiffUpper The upper equivalence limit of mean difference.
#' @param meanDiff The mean difference under the alternative
#'   hypothesis.
#' @param stDev The standard deviation.
#' @param allocationRatioPlanned Allocation ratio for the active treatment
#'   versus control. Defaults to 1 for equal randomization.
#' @param normalApproximation The type of computation of the p-values.
#'   If \code{TRUE}, the variance is assumed to be known, otherwise
#'   the calculations are performed with the t distribution. The exact
#'   calculation using the t distribution is only implemented for the
#'   fixed design.
#' @param rounding Whether to round up sample size. Defaults to 1 for
#'   sample size rounding.
#' @inheritParams param_kMax
#' @param informationRates The information rates. Fixed prior to the trial.
#'   Defaults to \code{(1:kMax) / kMax} if left unspecified.
#' @inheritParams param_criticalValues
#' @param alpha The significance level for each of the two one-sided
#'   tests. Defaults to 0.05.
#' @inheritParams param_typeAlphaSpending
#' @inheritParams param_parameterAlphaSpending
#' @inheritParams param_userAlphaSpending
#' @param spendingTime A vector of length \code{kMax} for the error spending
#'   time at each analysis. Defaults to missing, in which case, it is the
#'   same as \code{informationRates}.
#'
#' @return An S3 class \code{designMeanDiffEquiv} object with three
#' components:
#'
#' * \code{overallResults}: A data frame containing the following variables:
#'
#'     - \code{overallReject}: The overall rejection probability.
#'
#'     - \code{alpha}: The significance level for each of the two one-sided
#'       tests. Defaults to 0.05.
#'
#'     - \code{attainedAlpha}: The attained significance level.
#'
#'     - \code{kMax}: The number of stages.
#'
#'     - \code{information}: The maximum information.
#'
#'     - \code{expectedInformationH1}: The expected information under H1.
#'
#'     - \code{expectedInformationH0}: The expected information under H0.
#'
#'     - \code{numberOfSubjects}: The maximum number of subjects.
#'
#'     - \code{expectedNumberOfSubjectsH1}: The expected number of subjects
#'       under H1.
#'
#'     - \code{expectedNumberOfSubjectsH0}: The expected number of subjects
#'       under H0.
#'
#'     - \code{meanDiffLower}: The lower equivalence limit of mean
#'       difference.
#'
#'     - \code{meanDiffUpper}: The upper equivalence limit of mean
#'       difference.
#'
#'     - \code{meanDiff}: The mean difference under the alternative
#'       hypothesis.
#'
#'     - \code{stDev}: The standard deviation.
#'
#' * \code{byStageResults}: A data frame containing the following variables:
#'
#'     - \code{informationRates}: The information rates.
#'
#'     - \code{efficacyBounds}: The efficacy boundaries on the Z-scale for
#'       each of the two one-sided tests.
#'
#'     - \code{rejectPerStage}: The probability for efficacy stopping.
#'
#'     - \code{cumulativeRejection}: The cumulative probability for efficacy
#'       stopping.
#'
#'     - \code{cumulativeAlphaSpent}: The cumulative alpha for each of
#'       the two one-sided tests.
#'
#'     - \code{cumulativeAttainedAlpha}: The cumulative probability for
#'       efficacy stopping under H0.
#'
#'     - \code{efficacyMeanDiffLower}: The efficacy boundaries on the
#'       mean difference scale for the one-sided null hypothesis on the
#'       lower equivalence limit.
#'
#'     - \code{efficacyMeanDiffUpper}: The efficacy boundaries on the
#'       mean difference scale for the one-sided null hypothesis on the
#'       upper equivalence limit.
#'
#'     - \code{efficacyP}: The efficacy bounds on the p-value scale for
#'       each of the two one-sided tests.
#'
#'     - \code{information}: The cumulative information.
#'
#'     - \code{numberOfSubjects}: The number of subjects.
#'
#' * \code{settings}: A list containing the following input parameters:
#'
#'     - \code{typeAlphaSpending}: The type of alpha spending.
#'
#'     - \code{parameterAlphaSpending}: The parameter value for alpha
#'       spending.
#'
#'     - \code{userAlphaSpending}: The user defined alpha spending.
#'
#'     - \code{spendingTime}: The error spending time at each analysis.
#'
#'     - \code{allocationRatioPlanned}: Allocation ratio for the active
#'       treatment versus control.
#'
#'     - \code{normalApproximation}: The type of computation of the p-values.
#'       If \code{TRUE}, the variance is assumed to be known, otherwise
#'       the calculations are performed with the t distribution. The exact
#'       calculation using the t distribution is only implemented for the
#'       fixed design.
#'
#'     - \code{rounding}: Whether to round up sample size.
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @examples
#'
#' # Example 1: group sequential trial power calculation
#' (design1 <- getDesignMeanDiffEquiv(
#'   beta = 0.1, n = NA, meanDiffLower = -1.3, meanDiffUpper = 1.3,
#'   meanDiff = 0, stDev = 2.2,
#'   kMax = 4, alpha = 0.05, typeAlphaSpending = "sfOF"))
#'
#' # Example 2: sample size calculation for t-test
#' (design2 <- getDesignMeanDiffEquiv(
#'   beta = 0.1, n = NA, meanDiffLower = -1.3, meanDiffUpper = 1.3,
#'   meanDiff = 0, stDev = 2.2,
#'   normalApproximation = FALSE, alpha = 0.05))
#'
#' @export
getDesignMeanDiffEquiv <- function(
    beta = NA_real_,
    n = NA_real_,
    meanDiffLower = NA_real_,
    meanDiffUpper = NA_real_,
    meanDiff = 0,
    stDev = 1,
    allocationRatioPlanned = 1,
    normalApproximation = TRUE,
    rounding = TRUE,
    kMax = 1L,
    informationRates = NA_real_,
    alpha = 0.05,
    typeAlphaSpending = "sfOF",
    parameterAlphaSpending = NA_real_,
    userAlphaSpending = NA_real_,
    spendingTime = NA_real_) {

  if (is.na(beta) && is.na(n)) {
    stop("beta and n cannot be both missing")
  }

  if (!is.na(beta) && !is.na(n)) {
    stop("Only one of beta and n should be provided")
  }

  if (!is.na(beta) && (beta >= 1-alpha || beta < 0.0001)) {
    stop("beta must lie in [0.0001, 1-alpha)")
  }

  if (!is.na(n) && n <= 0) {
    stop("n must be positive")
  }

  if (is.na(meanDiffLower)) {
    stop("meanDiffLower must be provided")
  }

  if (is.na(meanDiffUpper)) {
    stop("meanDiffUpper must be provided")
  }

  if (meanDiffLower >= meanDiff) {
    stop("meanDiffLower must be less than meanDiff")
  }

  if (meanDiffUpper <= meanDiff) {
    stop("meanDiffUpper must be greater than meanDiff")
  }

  if (stDev <= 0) {
    stop("stDev must be positive")
  }

  if (alpha < 0.00001 || alpha >= 1) {
    stop("alpha must lie in [0.00001, 1)")
  }

  if (allocationRatioPlanned <= 0) {
    stop("allocationRatioPlanned must be positive")
  }

  if (any(is.na(informationRates))) {
    informationRates = (1:kMax)/kMax
  }


  # variance for one sampling unit
  r = allocationRatioPlanned/(1 + allocationRatioPlanned)
  v1 = stDev^2/(r*(1-r))

  f <- function(n) { # power for two one-sided t-tests
    b = qt(1-alpha, n-2)
    ncpLower = (meanDiff - meanDiffLower)*sqrt(n/v1)
    powerLower = pt(b, n-2, ncpLower, lower.tail = FALSE)
    ncpUpper = (meanDiffUpper - meanDiff)*sqrt(n/v1)
    powerUpper = pt(b, n-2, ncpUpper, lower.tail = FALSE)
    power = powerLower + powerUpper - 1
    power
  }

  if (is.na(n)) { # calculate sample size
    des = getDesignEquiv(
      beta = beta, IMax = NA, thetaLower = meanDiffLower,
      thetaUpper = meanDiffUpper, theta = meanDiff,
      kMax = kMax, informationRates = informationRates,
      alpha = alpha, typeAlphaSpending = typeAlphaSpending,
      parameterAlphaSpending = parameterAlphaSpending,
      userAlphaSpending = userAlphaSpending,
      spendingTime = spendingTime)

    n = des$overallResults$information*v1

    if (kMax == 1 && !normalApproximation) { # t-test for fixed design
      n = uniroot(function(n) f(n) - (1-beta), c(0.5*n, 1.5*n))$root
    }
  }

  if (rounding) {
    n = ceiling(n - 1.0e-12)
    informationRates = round(n*informationRates)/n
  }

  if (is.na(beta) || rounding) { # calculate power
    des = getDesignEquiv(
      beta = NA, IMax = n/v1, thetaLower = meanDiffLower,
      thetaUpper = meanDiffUpper, theta = meanDiff,
      kMax = kMax, informationRates = informationRates,
      alpha = alpha, typeAlphaSpending = typeAlphaSpending,
      parameterAlphaSpending = parameterAlphaSpending,
      userAlphaSpending = userAlphaSpending,
      spendingTime = spendingTime)
  }

  if (kMax == 1 && !normalApproximation) { # t-test for fixed design
    power = f(n)

    b = qt(1-alpha, n-2)
    ncp = (meanDiffUpper - meanDiffLower)*sqrt(n/v1)
    attainedAlpha = alpha - pt(b, n-2, ncp)

    des$overallResults$overallReject = power
    des$overallResults$attainedAlpha = attainedAlpha
    des$overallResults$information = n/v1
    des$overallResults$expectedInformationH1 = n/v1
    des$overallResults$expectedInformationH0 = n/v1

    des$byStageResults$efficacyBounds = b
    des$byStageResults$rejectPerStage = power
    des$byStageResults$cumulativeRejection = power
    des$byStageResults$cumulativeAttainedAlpha = attainedAlpha
    des$byStageResults$efficacyMeanDiffLower = b*sqrt(v1/n) +
      meanDiffLower
    des$byStageResults$efficacyMeanDiffUpper = -b*sqrt(v1/n) +
      meanDiffUpper
    des$byStageResults$information = n/v1
  } else {
    des$overallResults$attainedAlpha =
      des$overallResults$attainedAlphaH10
    des$overallResults$expectedInformationH0 =
      des$overallResults$expectedInformationH10

    des$byStageResults$cumulativeAttainedAlpha =
      des$byStageResults$cumulativeAttainedAlphaH10
    des$byStageResults$efficacyMeanDiffLower =
      des$byStageResults$efficacyThetaLower
    des$byStageResults$efficacyMeanDiffUpper =
      des$byStageResults$efficacyThetaUpper
  }


  des$overallResults$numberOfSubjects = n
  des$overallResults$expectedNumberOfSubjectsH1 =
    des$overallResults$expectedInformationH1*v1
  des$overallResults$expectedNumberOfSubjectsH0 =
    des$overallResults$expectedInformationH0*v1
  des$overallResults$meanDiffLower = meanDiffLower
  des$overallResults$meanDiffUpper = meanDiffUpper
  des$overallResults$meanDiff = meanDiff
  des$overallResults$stDev = stDev
  des$overallResults <-
    des$overallResults[, c("overallReject", "alpha", "attainedAlpha",
                           "kMax", "information", "expectedInformationH1",
                           "expectedInformationH0", "numberOfSubjects",
                           "expectedNumberOfSubjectsH1",
                           "expectedNumberOfSubjectsH0", "meanDiffLower",
                           "meanDiffUpper", "meanDiff", "stDev")]

  des$byStageResults$numberOfSubjects = n*informationRates
  des$byStageResults <-
    des$byStageResults[, c("informationRates", "efficacyBounds",
                           "rejectPerStage", "cumulativeRejection",
                           "cumulativeAlphaSpent", "cumulativeAttainedAlpha",
                           "efficacyMeanDiffLower",
                           "efficacyMeanDiffUpper", "efficacyP",
                           "information", "numberOfSubjects")]

  des$settings$allocationRatioPlanned = allocationRatioPlanned
  des$settings$normalApproximation = normalApproximation
  des$settings$rounding = rounding
  des$settings <-
    des$settings[c("typeAlphaSpending", "parameterAlphaSpending",
                   "userAlphaSpending", "spendingTime",
                   "allocationRatioPlanned",
                   "normalApproximation", "rounding")]

  attr(des, "class") = "designMeanDiffEquiv"

  des
}


#' @title Group Sequential Design for Equivalence in Two-Sample
#' Mean Ratio
#' @description Obtains the power given sample size or obtains the sample
#' size given power for a group sequential design for equivalence in
#' two-sample mean ratio.
#'
#' @param beta The type II error.
#' @param n The total sample size.
#' @param meanRatioLower The lower equivalence limit of mean ratio.
#' @param meanRatioUpper The upper equivalence limit of mean ratio.
#' @param meanRatio The mean ratio under the alternative hypothesis.
#' @param CV The coefficient of variation.
#' @param allocationRatioPlanned Allocation ratio for the active treatment
#'   versus control. Defaults to 1 for equal randomization.
#' @param normalApproximation The type of computation of the p-values.
#'   If \code{TRUE}, the variance is assumed to be known, otherwise
#'   the calculations are performed with the t distribution. The exact
#'   calculation using the t distribution is only implemented for the
#'   fixed design.
#' @param rounding Whether to round up sample size. Defaults to 1 for
#'   sample size rounding.
#' @inheritParams param_kMax
#' @param informationRates The information rates. Fixed prior to the trial.
#'   Defaults to \code{(1:kMax) / kMax} if left unspecified.
#' @inheritParams param_criticalValues
#' @param alpha The significance level for each of the two one-sided
#'   tests. Defaults to 0.05.
#' @inheritParams param_typeAlphaSpending
#' @inheritParams param_parameterAlphaSpending
#' @inheritParams param_userAlphaSpending
#' @param spendingTime A vector of length \code{kMax} for the error spending
#'   time at each analysis. Defaults to missing, in which case, it is the
#'   same as \code{informationRates}.
#'
#' @return An S3 class \code{designMeanRatioEquiv} object with three
#' components:
#'
#' * \code{overallResults}: A data frame containing the following variables:
#'
#'     - \code{overallReject}: The overall rejection probability.
#'
#'     - \code{alpha}: The significance level for each of the two one-sided
#'       tests. Defaults to 0.05.
#'
#'     - \code{attainedAlpha}: The attained significance level.
#'
#'     - \code{kMax}: The number of stages.
#'
#'     - \code{information}: The maximum information.
#'
#'     - \code{expectedInformationH1}: The expected information under H1.
#'
#'     - \code{expectedInformationH0}: The expected information under H0.
#'
#'     - \code{numberOfSubjects}: The maximum number of subjects.
#'
#'     - \code{expectedNumberOfSubjectsH1}: The expected number of subjects
#'       under H1.
#'
#'     - \code{expectedNumberOfSubjectsH0}: The expected number of subjects
#'       under H0.
#'
#'     - \code{meanRatioLower}: The lower equivalence limit of mean ratio.
#'
#'     - \code{meanRatioUpper}: The upper equivalence limit of mean ratio.
#'
#'     - \code{meanRatio}: The mean ratio under the alternative hypothesis.
#'
#'     - \code{CV}: The coefficient of variation.
#'
#' * \code{byStageResults}: A data frame containing the following variables:
#'
#'     - \code{informationRates}: The information rates.
#'
#'     - \code{efficacyBounds}: The efficacy boundaries on the Z-scale for
#'       each of the two one-sided tests.
#'
#'     - \code{rejectPerStage}: The probability for efficacy stopping.
#'
#'     - \code{cumulativeRejection}: The cumulative probability for efficacy
#'       stopping.
#'
#'     - \code{cumulativeAlphaSpent}: The cumulative alpha for each of
#'       the two one-sided tests.
#'
#'     - \code{cumulativeAttainedAlpha}: The cumulative probability for
#'       efficacy stopping under H0.
#'
#'     - \code{efficacyP}: The efficacy bounds on the p-value scale for
#'       each of the two one-sided tests.
#'
#'     - \code{information}: The cumulative information.
#'
#'     - \code{numberOfSubjects}: The number of subjects.
#'
#'     - \code{efficacyMeanRatioLower}: The efficacy boundaries on the
#'       mean ratio scale for the one-sided null hypothesis on the
#'       lower equivalence limit.
#'
#'     - \code{efficacyMeanRatioUpper}: The efficacy boundaries on the
#'       mean ratio scale for the one-sided null hypothesis on the
#'       upper equivalence limit.
#'
#' * \code{settings}: A list containing the following input parameters:
#'
#'     - \code{typeAlphaSpending}: The type of alpha spending.
#'
#'     - \code{parameterAlphaSpending}: The parameter value for alpha
#'       spending.
#'
#'     - \code{userAlphaSpending}: The user defined alpha spending.
#'
#'     - \code{spendingTime}: The error spending time at each analysis.
#'
#'     - \code{allocationRatioPlanned}: Allocation ratio for the active
#'       treatment versus control.
#'
#'     - \code{normalApproximation}: The type of computation of the p-values.
#'       If \code{TRUE}, the variance is assumed to be known, otherwise
#'       the calculations are performed with the t distribution. The exact
#'       calculation using the t distribution is only implemented for the
#'       fixed design.
#'
#'     - \code{rounding}: Whether to round up sample size.
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @examples
#'
#' # Example 1: group sequential trial power calculation
#' (design1 <- getDesignMeanRatioEquiv(
#'   beta = 0.1, n = NA, meanRatioLower = 0.8, meanRatioUpper = 1.25,
#'   meanRatio = 1, CV = 0.35,
#'   kMax = 4, alpha = 0.05, typeAlphaSpending = "sfOF"))
#'
#' # Example 2: sample size calculation for t-test
#' (design2 <- getDesignMeanRatioEquiv(
#'   beta = 0.1, n = NA, meanRatioLower = 0.8, meanRatioUpper = 1.25,
#'   meanRatio = 1, CV = 0.35,
#'   normalApproximation = FALSE, alpha = 0.05))
#'
#' @export
getDesignMeanRatioEquiv <- function(
    beta = NA_real_,
    n = NA_real_,
    meanRatioLower = NA_real_,
    meanRatioUpper = NA_real_,
    meanRatio = 1,
    CV = 1,
    allocationRatioPlanned = 1,
    normalApproximation = TRUE,
    rounding = TRUE,
    kMax = 1L,
    informationRates = NA_real_,
    alpha = 0.05,
    typeAlphaSpending = "sfOF",
    parameterAlphaSpending = NA_real_,
    userAlphaSpending = NA_real_,
    spendingTime = NA_real_) {

  if (is.na(meanRatioLower)) {
    stop("meanRatioLower must be provided")
  }

  if (is.na(meanRatioUpper)) {
    stop("meanRatioUpper must be provided")
  }

  if (meanRatioLower <= 0) {
    stop("meanRatioLower must be positive")
  }

  if (meanRatioLower >= meanRatio) {
    stop("meanRatioLower must be less than meanRatio")
  }

  if (meanRatioUpper <= meanRatio) {
    stop("meanRatioUpper must be greater than meanRatio")
  }

  if (CV <= 0) {
    stop("CV must be positive")
  }

  if (allocationRatioPlanned <= 0) {
    stop("allocationRatioPlanned must be positive")
  }

  if (alpha < 0.00001 || alpha >= 1) {
    stop("alpha must lie in [0.00001, 1)")
  }


  des = getDesignMeanDiffEquiv(
    beta, n, meanDiffLower = log(meanRatioLower),
    meanDiffUpper = log(meanRatioUpper),
    meanDiff = log(meanRatio),
    stDev = sqrt(log(1 + CV^2)),
    allocationRatioPlanned,
    normalApproximation, rounding,
    kMax, informationRates,
    alpha, typeAlphaSpending,
    parameterAlphaSpending, userAlphaSpending,
    spendingTime)

  des$overallResults$meanRatioLower = meanRatioLower
  des$overallResults$meanRatioUpper = meanRatioUpper
  des$overallResults$meanRatio = meanRatio
  des$overallResults$CV = CV

  des$overallResults$meanDiffLower = NULL
  des$overallResults$meanDiffUpper = NULL
  des$overallResults$meanDiff = NULL
  des$overallResults$stDev = NULL

  des$byStageResults$efficacyMeanRatioLower =
    exp(des$byStageResults$efficacyMeanDiffLower)
  des$byStageResults$efficacyMeanRatioUpper =
    exp(des$byStageResults$efficacyMeanDiffUpper)

  des$byStageResults$efficacyMeanDiffLower = NULL
  des$byStageResults$efficacyMeanDiffUpper = NULL

  attr(des, "class") = "designMeanRatioEquiv"

  des
}


#' @title Group Sequential Design for Equivalence in Mean Difference
#' in 2x2 Crossover
#' @description Obtains the power given sample size or obtains the sample
#' size given power for a group sequential design for equivalence in
#' mean difference in 2x2 crossover.
#'
#' @param beta The type II error.
#' @param n The total sample size.
#' @param meanDiffLower The lower equivalence limit of mean difference.
#' @param meanDiffUpper The upper equivalence limit of mean difference.
#' @param meanDiff The mean difference under the alternative
#'   hypothesis.
#' @param stDev The standard deviation for within-subject random error.
#' @param allocationRatioPlanned Allocation ratio for sequence A/B
#'   versus sequence B/A. Defaults to 1 for equal randomization.
#' @param normalApproximation The type of computation of the p-values.
#'   If \code{TRUE}, the variance is assumed to be known, otherwise
#'   the calculations are performed with the t distribution. The exact
#'   calculation using the t distribution is only implemented for the
#'   fixed design.
#' @param rounding Whether to round up sample size. Defaults to 1 for
#'   sample size rounding.
#' @inheritParams param_kMax
#' @param informationRates The information rates. Fixed prior to the trial.
#'   Defaults to \code{(1:kMax) / kMax} if left unspecified.
#' @inheritParams param_criticalValues
#' @param alpha The significance level for each of the two one-sided
#'   tests. Defaults to 0.05.
#' @inheritParams param_typeAlphaSpending
#' @inheritParams param_parameterAlphaSpending
#' @inheritParams param_userAlphaSpending
#' @param spendingTime A vector of length \code{kMax} for the error spending
#'   time at each analysis. Defaults to missing, in which case, it is the
#'   same as \code{informationRates}.
#'
#' @return An S3 class \code{designMeanDiffXOEquiv} object with three
#' components:
#'
#' * \code{overallResults}: A data frame containing the following variables:
#'
#'     - \code{overallReject}: The overall rejection probability.
#'
#'     - \code{alpha}: The overall significance level.
#'
#'     - \code{attainedAlpha}: The attained significance level.
#'
#'     - \code{kMax}: The number of stages.
#'
#'     - \code{information}: The maximum information.
#'
#'     - \code{expectedInformationH1}: The expected information under H1.
#'
#'     - \code{expectedInformationH0}: The expected information under H0.
#'
#'     - \code{numberOfSubjects}: The maximum number of subjects.
#'
#'     - \code{expectedNumberOfSubjectsH1}: The expected number of subjects
#'       under H1.
#'
#'     - \code{expectedNumberOfSubjectsH0}: The expected number of subjects
#'       under H0.
#'
#'     - \code{meanDiffLower}: The lower equivalence limit of mean
#'       difference.
#'
#'     - \code{meanDiffUpper}: The upper equivalence limit of mean
#'       difference.
#'
#'     - \code{meanDiff}: The mean difference under the alternative
#'       hypothesis.
#'
#'     - \code{stDev}: The standard deviation for within-subject random
#'       error.
#'
#' * \code{byStageResults}: A data frame containing the following variables:
#'
#'     - \code{informationRates}: The information rates.
#'
#'     - \code{efficacyBounds}: The efficacy boundaries on the Z-scale for
#'       each of the two one-sided tests.
#'
#'     - \code{rejectPerStage}: The probability for efficacy stopping.
#'
#'     - \code{cumulativeRejection}: The cumulative probability for efficacy
#'       stopping.
#'
#'     - \code{cumulativeAlphaSpent}: The cumulative alpha for each of
#'       the two one-sided tests.
#'
#'     - \code{cumulativeAttainedAlpha}: The cumulative probability for
#'       efficacy stopping under H0.
#'
#'     - \code{efficacyMeanDiffLower}: The efficacy boundaries on the
#'       mean difference scale for the one-sided null hypothesis on the
#'       lower equivalence limit.
#'
#'     - \code{efficacyMeanDiffUpper}: The efficacy boundaries on the
#'       mean difference scale for the one-sided null hypothesis on the
#'       upper equivalence limit.
#'
#'     - \code{efficacyP}: The efficacy bounds on the p-value scale for
#'       each of the two one-sided tests.
#'
#'     - \code{information}: The cumulative information.
#'
#'     - \code{numberOfSubjects}: The number of subjects.
#'
#' * \code{settings}: A list containing the following input parameters:
#'
#'     - \code{typeAlphaSpending}: The type of alpha spending.
#'
#'     - \code{parameterAlphaSpending}: The parameter value for alpha
#'       spending.
#'
#'     - \code{userAlphaSpending}: The user defined alpha spending.
#'
#'     - \code{spendingTime}: The error spending time at each analysis.
#'
#'     - \code{allocationRatioPlanned}: Allocation ratio for sequence A/B
#'       versus sequence B/A.
#'
#'     - \code{normalApproximation}: The type of computation of the p-values.
#'       If \code{TRUE}, the variance is assumed to be known, otherwise
#'       the calculations are performed with the t distribution. The exact
#'       calculation using the t distribution is only implemented for the
#'       fixed design.
#'
#'     - \code{rounding}: Whether to round up sample size.
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @examples
#'
#' # Example 1: group sequential trial power calculation
#' (design1 <- getDesignMeanDiffXOEquiv(
#'   beta = 0.1, n = NA, meanDiffLower = -1.3, meanDiffUpper = 1.3,
#'   meanDiff = 0, stDev = 2.2,
#'   kMax = 4, alpha = 0.05, typeAlphaSpending = "sfOF"))
#'
#' # Example 2: sample size calculation for t-test
#' (design2 <- getDesignMeanDiffXOEquiv(
#'   beta = 0.1, n = NA, meanDiffLower = -1.3, meanDiffUpper = 1.3,
#'   meanDiff = 0, stDev = 2.2,
#'   normalApproximation = FALSE, alpha = 0.05))
#'
#' @export
getDesignMeanDiffXOEquiv <- function(
    beta = NA_real_,
    n = NA_real_,
    meanDiffLower = NA_real_,
    meanDiffUpper = NA_real_,
    meanDiff = 0,
    stDev = 1,
    allocationRatioPlanned = 1,
    normalApproximation = TRUE,
    rounding = TRUE,
    kMax = 1L,
    informationRates = NA_real_,
    alpha = 0.05,
    typeAlphaSpending = "sfOF",
    parameterAlphaSpending = NA_real_,
    userAlphaSpending = NA_real_,
    spendingTime = NA_real_) {

  if (is.na(beta) && is.na(n)) {
    stop("beta and n cannot be both missing")
  }

  if (!is.na(beta) && !is.na(n)) {
    stop("Only one of beta and n should be provided")
  }

  if (!is.na(beta) && (beta >= 1-alpha || beta < 0.0001)) {
    stop("beta must lie in [0.0001, 1-alpha)")
  }

  if (!is.na(n) && n <= 0) {
    stop("n must be positive")
  }

  if (is.na(meanDiffLower)) {
    stop("meanDiffLower must be provided")
  }

  if (is.na(meanDiffUpper)) {
    stop("meanDiffUpper must be provided")
  }

  if (meanDiffLower >= meanDiff) {
    stop("meanDiffLower must be less than meanDiff")
  }

  if (meanDiffUpper <= meanDiff) {
    stop("meanDiffUpper must be greater than meanDiff")
  }

  if (stDev <= 0) {
    stop("stDev must be positive")
  }

  if (allocationRatioPlanned <= 0) {
    stop("allocationRatioPlanned must be positive")
  }

  if (alpha < 0.00001 || alpha >= 1) {
    stop("alpha must lie in [0.00001, 1)")
  }

  if (any(is.na(informationRates))) {
    informationRates = (1:kMax)/kMax
  }


  # variance for one sampling unit
  r = allocationRatioPlanned/(1 + allocationRatioPlanned)
  v1 = stDev^2/(2*r*(1-r))

  f <- function(n) { # power for two one-sided t-tests
    b = qt(1-alpha, n-2)
    ncpLower = (meanDiff - meanDiffLower)*sqrt(n/v1)
    powerLower = pt(b, n-2, ncpLower, lower.tail = FALSE)
    ncpUpper = (meanDiffUpper - meanDiff)*sqrt(n/v1)
    powerUpper = pt(b, n-2, ncpUpper, lower.tail = FALSE)
    power = powerLower + powerUpper - 1
    power
  }

  if (is.na(n)) { # calculate sample size
    des = getDesignEquiv(
      beta = beta, IMax = NA, thetaLower = meanDiffLower,
      thetaUpper = meanDiffUpper, theta = meanDiff,
      kMax = kMax, informationRates = informationRates,
      alpha = alpha, typeAlphaSpending = typeAlphaSpending,
      parameterAlphaSpending = parameterAlphaSpending,
      userAlphaSpending = userAlphaSpending,
      spendingTime = spendingTime)

    n = des$overallResults$information*v1

    if (kMax == 1 && !normalApproximation) { # t-test for fixed design
      n = uniroot(function(n) f(n) - (1-beta), c(0.5*n, 1.5*n))$root
    }
  }

  if (rounding) {
    n = ceiling(n - 1.0e-12)
    informationRates = round(n*informationRates)/n
  }

  if (is.na(beta) || rounding) { # calculate power
    des = getDesignEquiv(
      beta = NA, IMax = n/v1, thetaLower = meanDiffLower,
      thetaUpper = meanDiffUpper, theta = meanDiff,
      kMax = kMax, informationRates = informationRates,
      alpha = alpha, typeAlphaSpending = typeAlphaSpending,
      parameterAlphaSpending = parameterAlphaSpending,
      userAlphaSpending = userAlphaSpending,
      spendingTime = spendingTime)
  }

  if (kMax == 1 && !normalApproximation) { # t-test for fixed design
    power = f(n)

    b = qt(1-alpha, n-2)
    ncp = (meanDiffUpper - meanDiffLower)*sqrt(n/v1)
    attainedAlpha = alpha - pt(b, n-2, ncp)

    des$overallResults$overallReject = power
    des$overallResults$attainedAlpha = attainedAlpha
    des$overallResults$information = n/v1
    des$overallResults$expectedInformationH1 = n/v1
    des$overallResults$expectedInformationH0 = n/v1

    des$byStageResults$efficacyBounds = b
    des$byStageResults$rejectPerStage = power
    des$byStageResults$cumulativeRejection = power
    des$byStageResults$cumulativeAttainedAlpha = attainedAlpha
    des$byStageResults$efficacyMeanDiffLower = b*sqrt(v1/n) +
      meanDiffLower
    des$byStageResults$efficacyMeanDiffUpper = -b*sqrt(v1/n) +
      meanDiffUpper
    des$byStageResults$information = n/v1
  } else {
    des$overallResults$attainedAlpha =
      des$overallResults$attainedAlphaH10
    des$overallResults$expectedInformationH0 =
      des$overallResults$expectedInformationH10

    des$byStageResults$cumulativeAttainedAlpha =
      des$byStageResults$cumulativeAttainedAlphaH10
    des$byStageResults$efficacyMeanDiffLower =
      des$byStageResults$efficacyThetaLower
    des$byStageResults$efficacyMeanDiffUpper =
      des$byStageResults$efficacyThetaUpper
  }


  des$overallResults$numberOfSubjects = n
  des$overallResults$expectedNumberOfSubjectsH1 =
    des$overallResults$expectedInformationH1*v1
  des$overallResults$expectedNumberOfSubjectsH0 =
    des$overallResults$expectedInformationH0*v1
  des$overallResults$meanDiffLower = meanDiffLower
  des$overallResults$meanDiffUpper = meanDiffUpper
  des$overallResults$meanDiff = meanDiff
  des$overallResults$stDev = stDev
  des$overallResults <-
    des$overallResults[, c("overallReject", "alpha", "attainedAlpha",
                           "kMax", "information", "expectedInformationH1",
                           "expectedInformationH0", "numberOfSubjects",
                           "expectedNumberOfSubjectsH1",
                           "expectedNumberOfSubjectsH0", "meanDiffLower",
                           "meanDiffUpper", "meanDiff", "stDev")]

  des$byStageResults$numberOfSubjects = n*informationRates
  des$byStageResults <-
    des$byStageResults[, c("informationRates", "efficacyBounds",
                           "rejectPerStage", "cumulativeRejection",
                           "cumulativeAlphaSpent", "cumulativeAttainedAlpha",
                           "efficacyMeanDiffLower",
                           "efficacyMeanDiffUpper", "efficacyP",
                           "information", "numberOfSubjects")]

  des$settings$allocationRatioPlanned = allocationRatioPlanned
  des$settings$normalApproximation = normalApproximation
  des$settings$rounding = rounding
  des$settings <-
    des$settings[c("typeAlphaSpending", "parameterAlphaSpending",
                   "userAlphaSpending", "spendingTime",
                   "allocationRatioPlanned",
                   "normalApproximation", "rounding")]

  attr(des, "class") = "designMeanDiffXOEquiv"

  des
}


#' @title Group Sequential Design for Equivalence in Mean Ratio
#' in 2x2 Crossover
#' @description Obtains the power given sample size or obtains the sample
#' size given power for a group sequential design for equivalence
#' mean ratio in 2x2 crossover.
#'
#' @param beta The type II error.
#' @param n The total sample size.
#' @param meanRatioLower The lower equivalence limit of mean ratio.
#' @param meanRatioUpper The upper equivalence limit of mean ratio.
#' @param meanRatio The mean ratio under the alternative hypothesis.
#' @param CV The coefficient of variation.
#' @param allocationRatioPlanned Allocation ratio for sequence A/B
#'   versus sequence B/A. Defaults to 1 for equal randomization.
#' @param normalApproximation The type of computation of the p-values.
#'   If \code{TRUE}, the variance is assumed to be known, otherwise
#'   the calculations are performed with the t distribution. The exact
#'   calculation using the t distribution is only implemented for the
#'   fixed design.
#' @param rounding Whether to round up sample size. Defaults to 1 for
#'   sample size rounding.
#' @inheritParams param_kMax
#' @param informationRates The information rates. Fixed prior to the trial.
#'   Defaults to \code{(1:kMax) / kMax} if left unspecified.
#' @inheritParams param_criticalValues
#' @param alpha The significance level for each of the two one-sided
#'   tests. Defaults to 0.05.
#' @inheritParams param_typeAlphaSpending
#' @inheritParams param_parameterAlphaSpending
#' @inheritParams param_userAlphaSpending
#' @param spendingTime A vector of length \code{kMax} for the error spending
#'   time at each analysis. Defaults to missing, in which case, it is the
#'   same as \code{informationRates}.
#'
#' @return An S3 class \code{designMeanRatioEquiv} object with three
#' components:
#'
#' * \code{overallResults}: A data frame containing the following variables:
#'
#'     - \code{overallReject}: The overall rejection probability.
#'
#'     - \code{alpha}: The overall significance level.
#'
#'     - \code{attainedAlpha}: The attained significance level.
#'
#'     - \code{kMax}: The number of stages.
#'
#'     - \code{information}: The maximum information.
#'
#'     - \code{expectedInformationH1}: The expected information under H1.
#'
#'     - \code{expectedInformationH0}: The expected information under H0.
#'
#'     - \code{numberOfSubjects}: The maximum number of subjects.
#'
#'     - \code{expectedNumberOfSubjectsH1}: The expected number of subjects
#'       under H1.
#'
#'     - \code{expectedNumberOfSubjectsH0}: The expected number of subjects
#'       under H0.
#'
#'     - \code{meanRatioLower}: The lower equivalence limit of mean ratio.
#'
#'     - \code{meanRatioUpper}: The upper equivalence limit of mean ratio.
#'
#'     - \code{meanRatio}: The mean ratio under the alternative hypothesis.
#'
#'     - \code{CV}: The coefficient of variation.
#'
#' * \code{byStageResults}: A data frame containing the following variables:
#'
#'     - \code{informationRates}: The information rates.
#'
#'     - \code{efficacyBounds}: The efficacy boundaries on the Z-scale for
#'       each of the two one-sided tests.
#'
#'     - \code{rejectPerStage}: The probability for efficacy stopping.
#'
#'     - \code{cumulativeRejection}: The cumulative probability for efficacy
#'       stopping.
#'
#'     - \code{cumulativeAlphaSpent}: The cumulative alpha for each of
#'       the two one-sided tests.
#'
#'     - \code{cumulativeAttainedAlpha}: The cumulative probability for
#'       efficacy stopping under H0.
#'
#'     - \code{efficacyMeanRatioLower}: The efficacy boundaries on the
#'       mean ratio scale for the one-sided null hypothesis on the
#'       lower equivalence limit.
#'
#'     - \code{efficacyMeanRatioUpper}: The efficacy boundaries on the
#'       mean ratio scale for the one-sided null hypothesis on the
#'       upper equivalence limit.
#'
#'     - \code{efficacyP}: The efficacy bounds on the p-value scale for
#'       each of the two one-sided tests.
#'
#'     - \code{information}: The cumulative information.
#'
#'     - \code{numberOfSubjects}: The number of subjects.
#'
#' * \code{settings}: A list containing the following input parameters:
#'
#'     - \code{typeAlphaSpending}: The type of alpha spending.
#'
#'     - \code{parameterAlphaSpending}: The parameter value for alpha
#'       spending.
#'
#'     - \code{userAlphaSpending}: The user defined alpha spending.
#'
#'     - \code{spendingTime}: The error spending time at each analysis.
#'
#'     - \code{allocationRatioPlanned}: Allocation ratio for sequence A/B
#'       versus sequence B/A.
#'
#'     - \code{normalApproximation}: The type of computation of the p-values.
#'       If \code{TRUE}, the variance is assumed to be known, otherwise
#'       the calculations are performed with the t distribution. The exact
#'       calculation using the t distribution is only implemented for the
#'       fixed design.
#'
#'     - \code{rounding}: Whether to round up sample size.
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @examples
#'
#' # Example 1: group sequential trial power calculation
#' (design1 <- getDesignMeanRatioXOEquiv(
#'   beta = 0.1, n = NA, meanRatioLower = 0.8, meanRatioUpper = 1.25,
#'   meanRatio = 1, CV = 0.35,
#'   kMax = 4, alpha = 0.05, typeAlphaSpending = "sfOF"))
#'
#' # Example 2: sample size calculation for t-test
#' (design2 <- getDesignMeanRatioXOEquiv(
#'   beta = 0.1, n = NA, meanRatioLower = 0.8, meanRatioUpper = 1.25,
#'   meanRatio = 1, CV = 0.35,
#'   normalApproximation = FALSE, alpha = 0.05))
#'
#' @export
getDesignMeanRatioXOEquiv <- function(
    beta = NA_real_,
    n = NA_real_,
    meanRatioLower = NA_real_,
    meanRatioUpper = NA_real_,
    meanRatio = 1,
    CV = 1,
    allocationRatioPlanned = 1,
    normalApproximation = TRUE,
    rounding = TRUE,
    kMax = 1L,
    informationRates = NA_real_,
    alpha = 0.05,
    typeAlphaSpending = "sfOF",
    parameterAlphaSpending = NA_real_,
    userAlphaSpending = NA_real_,
    spendingTime = NA_real_) {

  if (is.na(meanRatioLower)) {
    stop("meanRatioLower must be provided")
  }

  if (is.na(meanRatioUpper)) {
    stop("meanRatioUpper must be provided")
  }

  if (meanRatioLower <= 0) {
    stop("meanRatioLower must be positive")
  }

  if (meanRatioLower >= meanRatio) {
    stop("meanRatioLower must be less than meanRatio")
  }

  if (meanRatioUpper <= meanRatio) {
    stop("meanRatioUpper must be greater than meanRatio")
  }

  if (CV <= 0) {
    stop("CV must be positive")
  }

  if (allocationRatioPlanned <= 0) {
    stop("allocationRatioPlanned must be positive")
  }

  if (alpha < 0.00001 || alpha >= 1) {
    stop("alpha must lie in [0.00001, 1)")
  }


  des = getDesignMeanDiffXOEquiv(
    beta, n, meanDiffLower = log(meanRatioLower),
    meanDiffUpper = log(meanRatioUpper),
    meanDiff = log(meanRatio),
    stDev = sqrt(log(1 + CV^2)),
    allocationRatioPlanned,
    normalApproximation, rounding,
    kMax, informationRates,
    alpha, typeAlphaSpending,
    parameterAlphaSpending, userAlphaSpending,
    spendingTime)

  des$overallResults$meanRatioLower = meanRatioLower
  des$overallResults$meanRatioUpper = meanRatioUpper
  des$overallResults$meanRatio = meanRatio
  des$overallResults$CV = CV

  des$overallResults$meanDiffLower = NULL
  des$overallResults$meanDiffUpper = NULL
  des$overallResults$meanDiff = NULL
  des$overallResults$stDev = NULL

  des$byStageResults$efficacyMeanRatioLower =
    exp(des$byStageResults$efficacyMeanDiffLower)
  des$byStageResults$efficacyMeanRatioUpper =
    exp(des$byStageResults$efficacyMeanDiffUpper)

  des$byStageResults$efficacyMeanDiffLower = NULL
  des$byStageResults$efficacyMeanDiffUpper = NULL

  attr(des, "class") = "designMeanRatioXOEquiv"

  des
}


#' @title Group Sequential Design for Two-Sample Wilcoxon Test
#' @description Obtains the power given sample size or obtains the sample
#' size given power for a group sequential design for two-sample
#' Wilcoxon test.
#'
#' @param beta The type II error.
#' @param n The total sample size.
#' @param pLarger The probability that a randomly chosen sample from
#'   the treatment group is larger than a randomly chosen sample from the
#'   control group under the alternative hypothesis.
#' @param allocationRatioPlanned Allocation ratio for the active treatment
#'   versus control. Defaults to 1 for equal randomization.
#' @param rounding Whether to round up sample size. Defaults to 1 for
#'   sample size rounding.
#' @inheritParams param_kMax
#' @param informationRates The information rates. Fixed prior to the trial.
#'   Defaults to \code{(1:kMax) / kMax} if left unspecified.
#' @inheritParams param_efficacyStopping
#' @inheritParams param_futilityStopping
#' @inheritParams param_criticalValues
#' @inheritParams param_alpha
#' @inheritParams param_typeAlphaSpending
#' @inheritParams param_parameterAlphaSpending
#' @inheritParams param_userAlphaSpending
#' @inheritParams param_futilityBounds
#' @inheritParams param_typeBetaSpending
#' @inheritParams param_parameterBetaSpending
#' @inheritParams param_userBetaSpending
#' @param spendingTime A vector of length \code{kMax} for the error spending
#'   time at each analysis. Defaults to missing, in which case, it is the
#'   same as \code{informationRates}.
#'
#' @details
#' The Mann-Whitney U test is a non-parametric test for the
#' difference in location between two independent groups. It is
#' also known as the Wilcoxon rank-sum test. The test is based on
#' the ranks of the data rather than the actual values, making it
#' robust to outliers and non-normal distributions.
#' The test statistic is the number of times a randomly
#' chosen sample from the treatment group is larger than a
#' randomly chosen sample from the control group, i.e.,
#' \deqn{W_{XY} = \sum_{i=1}^{n_1}\sum_{j=1}^{n_2} I(X_i > Y_j)}
#' where \eqn{X_i} and \eqn{Y_j} are the samples from the
#' treatment and control groups, respectively.
#' The test is often used when the data
#' do not meet the assumptions of the t-test, such as
#' non-normality or unequal variances.
#' The test is also applicable to ordinal data.
#' The test is one-sided, meaning that it only tests
#' whether the treatment group is larger than the control group.
#' Asymptotically,
#' \deqn{\frac{W_{XY} - n_1 n_2/2}{\sqrt{n_1 n_2 (n+1)/12}} \sim N(0,1)
#' \quad \text{under} H_0}
#' where \eqn{n_1} and \eqn{n_2} are the sample sizes of the treatment
#' and control groups, respectively, and \eqn{n=n_1+n_2}.
#' Let \eqn{\theta = P(X > Y)}, and
#' \eqn{\hat{\theta} = \frac{1}{nm}W_{XY}}. It follows that
#' \deqn{\sqrt{n}(\hat{\theta} - 1/2) \sim
#' N\left(0, \frac{1}{12 r(1-r)}\right) \quad \text{under} H_0}
#' where \eqn{r = n_1/(n_1+n_2)} is the randomization probability
#' for the active treatment group.
#' This formulation allows for group sequential testing with
#' futility stopping and efficacy stopping.
#'
#' @return An S3 class \code{designWilcoxon} object with three
#' components:
#'
#' * \code{overallResults}: A data frame containing the following variables:
#'
#'     - \code{overallReject}: The overall rejection probability.
#'
#'     - \code{alpha}: The overall significance level.
#'
#'     - \code{attainedAlpha}: The attained significance level, which is
#'       different from the overall significance level in the presence of
#'       futility stopping..
#'
#'     - \code{kMax}: The number of stages.
#'
#'     - \code{theta}: The parameter value.
#'
#'     - \code{information}: The maximum information.
#'
#'     - \code{expectedInformationH1}: The expected information under H1.
#'
#'     - \code{expectedInformationH0}: The expected information under H0.
#'
#'     - \code{drift}: The drift parameter, equal to
#'       \code{theta*sqrt(information)}.
#'
#'     - \code{inflationFactor}: The inflation factor (relative to the
#'       fixed design).
#'
#'     - \code{numberOfSubjects}: The maximum number of subjects.
#'
#'     - \code{expectedNumberOfSubjectsH1}: The expected number of subjects
#'       under H1.
#'
#'     - \code{expectedNumberOfSubjectsH0}: The expected number of subjects
#'       under H0.
#'
#'     - \code{pLarger}: The probability that a randomly chosen sample from
#'       the treatment group is larger than a randomly chosen sample from the
#'       control group under the alternative hypothesis.
#'
#' * \code{byStageResults}: A data frame containing the following variables:
#'
#'     - \code{informationRates}: The information rates.
#'
#'     - \code{efficacyBounds}: The efficacy boundaries on the Z-scale.
#'
#'     - \code{futilityBounds}: The futility boundaries on the Z-scale.
#'
#'     - \code{rejectPerStage}: The probability for efficacy stopping.
#'
#'     - \code{futilityPerStage}: The probability for futility stopping.
#'
#'     - \code{cumulativeRejection}: The cumulative probability for efficacy
#'       stopping.
#'
#'     - \code{cumulativeFutility}: The cumulative probability for futility
#'       stopping.
#'
#'     - \code{cumulativeAlphaSpent}: The cumulative alpha spent.
#'
#'     - \code{efficacyP}: The efficacy boundaries on the p-value scale.
#'
#'     - \code{futilityP}: The futility boundaries on the p-value scale.
#'
#'     - \code{information}: The cumulative information.
#'
#'     - \code{efficacyStopping}: Whether to allow efficacy stopping.
#'
#'     - \code{futilityStopping}: Whether to allow futility stopping.
#'
#'     - \code{rejectPerStageH0}: The probability for efficacy stopping
#'       under H0.
#'
#'     - \code{futilityPerStageH0}: The probability for futility stopping
#'       under H0.
#'
#'     - \code{cumulativeRejectionH0}: The cumulative probability for
#'       efficacy stopping under H0.
#'
#'     - \code{cumulativeFutilityH0}: The cumulative probability for futility
#'       stopping under H0.
#'
#'     - \code{efficacyPLarger}: The efficacy boundaries on the proportion
#'       of pairs of samples from the two treatment groups with the sample
#'       from the treatment group greater than that from the control group.
#'
#'     - \code{futilityPLarger}: The futility boundaries on the proportion
#'       of pairs of samples from the two treatment groups with the sample
#'       from the treatment group greater than that from the control group.
#'
#'     - \code{numberOfSubjects}: The number of subjects.
#'
#' * \code{settings}: A list containing the following input parameters:
#'
#'     - \code{typeAlphaSpending}: The type of alpha spending.
#'
#'     - \code{parameterAlphaSpending}: The parameter value for alpha
#'       spending.
#'
#'     - \code{userAlphaSpending}: The user defined alpha spending.
#'
#'     - \code{typeBetaSpending}: The type of beta spending.
#'
#'     - \code{parameterBetaSpending}: The parameter value for beta spending.
#'
#'     - \code{userBetaSpending}: The user defined beta spending.
#'
#'     - \code{spendingTime}: The error spending time at each analysis.
#'
#'     - \code{allocationRatioPlanned}: Allocation ratio for the active
#'       treatment versus control.
#'
#'     - \code{rounding}: Whether to round up sample size.
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @examples
#'
#' # Example 1: fixed design
#' (design1 <- getDesignWilcoxon(
#'   beta = 0.1, n = NA,
#'   pLarger = pnorm((8 - 2)/sqrt(2*25^2)), alpha = 0.025))
#'
#' # Example 2: group sequential design
#' (design2 <- getDesignWilcoxon(
#'   beta = 0.1, n = NA,
#'   pLarger = pnorm((8 - 2)/sqrt(2*25^2)), alpha = 0.025,
#'   kMax = 3, typeAlphaSpending = "sfOF"))
#'
#' @export
getDesignWilcoxon <- function(
    beta = NA_real_,
    n = NA_real_,
    pLarger = 0.6,
    allocationRatioPlanned = 1,
    rounding = TRUE,
    kMax = 1L,
    informationRates = NA_real_,
    efficacyStopping = NA_integer_,
    futilityStopping = NA_integer_,
    criticalValues = NA_real_,
    alpha = 0.025,
    typeAlphaSpending = "sfOF",
    parameterAlphaSpending = NA_real_,
    userAlphaSpending = NA_real_,
    futilityBounds = NA_real_,
    typeBetaSpending = "none",
    parameterBetaSpending = NA_real_,
    userBetaSpending = NA_real_,
    spendingTime = NA_real_) {

  if (is.na(beta) && is.na(n)) {
    stop("beta and n cannot be both missing")
  }

  if (!is.na(beta) && !is.na(n)) {
    stop("Only one of beta and n should be provided")
  }

  if (!is.na(beta) && (beta >= 1-alpha || beta < 0.0001)) {
    stop("beta must lie in [0.0001, 1-alpha)")
  }

  if (!is.na(n) && n <= 0) {
    stop("n must be positive")
  }

  if (pLarger <= 0 || pLarger >= 1) {
    stop("pLarger must lie between 0 and 1")
  }

  if (allocationRatioPlanned <= 0) {
    stop("allocationRatioPlanned must be positive")
  }

  if (alpha < 0.00001 || alpha >= 1) {
    stop("alpha must lie in [0.00001, 1)")
  }

  if (any(is.na(informationRates))) {
    informationRates = (1:kMax)/kMax
  }


  r = allocationRatioPlanned/(1 + allocationRatioPlanned)

  directionUpper = pLarger > 0.5

  theta = ifelse(directionUpper, pLarger - 0.5, 0.5 - pLarger)

  # variance for one sampling unit
  v1 = 1/(12*r*(1-r))

  if (is.na(beta)) { # power calculation
    if (rounding) {
      n = ceiling(n - 1.0e-12)
      informationRates = round(n*informationRates)/n
    }

    des = getDesign(
      beta = NA, IMax = n/v1, theta,
      kMax, informationRates,
      efficacyStopping, futilityStopping,
      criticalValues, alpha, typeAlphaSpending,
      parameterAlphaSpending, userAlphaSpending,
      futilityBounds, typeBetaSpending,
      parameterBetaSpending, userBetaSpending,
      spendingTime, 1)
  } else { # sample size calculation
    des = getDesign(
      beta, IMax = NA, theta,
      kMax, informationRates,
      efficacyStopping, futilityStopping,
      criticalValues, alpha, typeAlphaSpending,
      parameterAlphaSpending, userAlphaSpending,
      futilityBounds, typeBetaSpending,
      parameterBetaSpending, userBetaSpending,
      spendingTime, 1)

    n = des$overallResults$information*v1

    if (rounding) {
      n = ceiling(n - 1.0e-12)
      informationRates = des$byStageResults$informationRates
      informationRates = round(n*informationRates)/n

      des = getDesign(
        beta = NA, IMax = n/v1, theta,
        kMax, informationRates,
        efficacyStopping, futilityStopping,
        criticalValues, alpha, typeAlphaSpending,
        parameterAlphaSpending, userAlphaSpending,
        futilityBounds, typeBetaSpending,
        parameterBetaSpending, userBetaSpending,
        spendingTime, 1)
    }
  }

  des$overallResults$numberOfSubjects = n
  des$overallResults$expectedNumberOfSubjectsH1 =
    des$overallResults$expectedInformationH1*v1
  des$overallResults$expectedNumberOfSubjectsH0 =
    des$overallResults$expectedInformationH0*v1
  des$overallResults$pLarger = pLarger

  if (directionUpper) {
    des$byStageResults$efficacyPLarger =
      des$byStageResults$efficacyTheta + 0.5
    des$byStageResults$futilityPLarger =
      des$byStageResults$futilityTheta + 0.5
  } else {
    des$byStageResults$efficacyPLarger =
      -des$byStageResults$efficacyTheta + 0.5
    des$byStageResults$futilityPLarger =
      -des$byStageResults$futilityTheta + 0.5
  }
  des$byStageResults$efficacyTheta = NULL
  des$byStageResults$futilityTheta = NULL
  des$byStageResults$numberOfSubjects =
    des$byStageResults$informationRates*n

  des$settings$allocationRatioPlanned = allocationRatioPlanned
  des$settings$rounding = rounding
  des$settings$varianceRatio = NULL

  attr(des, "class") = "designWilcoxon"

  des
}


#' @title Group Sequential Design for Two-Sample Mean Difference From the
#' MMRM Model
#' @description Obtains the power and sample size for two-sample
#' mean difference at the last time point from the mixed-model
#' for repeated measures (MMRM) model.
#'
#' @param beta The type II error.
#' @param meanDiffH0 The mean difference at the last time point
#'   under the null hypothesis. Defaults to 0.
#' @param meanDiff The mean difference at the last time point
#'   under the alternative hypothesis.
#' @param k The number of postbaseline time points.
#' @param t The postbaseline time points.
#' @param covar1 The covariance matrix for the repeated measures
#'   given baseline for the active treatment group.
#' @param covar2 The covariance matrix for the repeated measures
#'   given baseline for the control group. If missing, it will be
#'   set equal to the covariance matrix for the active treatment group.
#' @inheritParams param_accrualTime
#' @inheritParams param_accrualIntensity
#' @inheritParams param_piecewiseSurvivalTime
#' @inheritParams param_gamma1
#' @inheritParams param_gamma2
#' @inheritParams param_accrualDuration
#' @inheritParams param_allocationRatioPlanned
#' @param normalApproximation The type of computation of the p-values.
#'   If \code{TRUE}, the variance is assumed to be known, otherwise
#'   the calculations are performed with the t distribution. The
#'   degrees of freedom for the t-distribution is the total effective
#'   sample size minus 2.  The exact calculation using the t distribution
#'   is only implemented for the fixed design.
#' @param rounding Whether to round up sample size. Defaults to 1 for
#'   sample size rounding.
#' @inheritParams param_kMax
#' @param informationRates The information rates. Defaults to
#'   \code{(1:kMax) / kMax} if left unspecified.
#' @inheritParams param_efficacyStopping
#' @inheritParams param_futilityStopping
#' @inheritParams param_criticalValues
#' @inheritParams param_alpha
#' @inheritParams param_typeAlphaSpending
#' @inheritParams param_parameterAlphaSpending
#' @inheritParams param_userAlphaSpending
#' @inheritParams param_futilityBounds
#' @inheritParams param_typeBetaSpending
#' @inheritParams param_parameterBetaSpending
#' @inheritParams param_userBetaSpending
#' @param spendingTime A vector of length \code{kMax} for the error spending
#'   time at each analysis. Defaults to missing, in which case, it is the
#'   same as \code{informationRates}.
#'
#' @details
#' Consider a longitudinal study with two treatment groups.
#' The outcome is measured at baseline and at \eqn{k} postbaseline
#' time points. For each treatment group, the outcomes are assumed to
#' follow a multivariate normal distribution.
#' Conditional on baseline, the covariance matrix of the post-baseline
#' outcomes is denoted by \eqn{\Sigma_1} for the active treatment group
#' and \eqn{\Sigma_2} for the control group.
#' Let \eqn{\mu_1} and \eqn{\mu_2} denote the mean vectors of
#' post-baseline outcomes for the active and control groups, respectively.
#' We are interested in testing the null
#' hypothesis \eqn{H_0: \mu_{1,k} - \mu_{2,k} = \delta_0}
#' against the alternative \eqn{H_1: \mu_{1,k} - \mu_{2,k} = \delta}.
#'
#' The study design is based on the information for treatment difference
#' at the last postbaseline time point. This information is given by
#' \deqn{I = 1/\text{Var}(\hat{\mu}_{1,k} - \hat{\mu}_{2,k})}
#' In the presence of monotone missing data, let
#' \eqn{p_{g,1},\ldots,p_{g,k}} denote the proportions of subjects in
#' observed data patterns 1 through \eqn{k} for treatment group
#' \eqn{g=1} (active) or 2 (control).
#' A subject in pattern \eqn{j} has complete data up to time \eqn{t_j}, i.e.,
#' the observed outcomes are \eqn{y_{i,1},\ldots,y_{i,j}},
#' with missing values for \eqn{y_{i,j+1},\ldots,y_{i,k}}.
#'
#' According to Lu et al. (2008), the information
#' matrix for the post-baseline mean vector in group \eqn{g} is
#' \deqn{I_g = n \pi_g J_g}
#' where \eqn{\pi_g} is the proportion of subjects in group \eqn{g},
#' and \deqn{J_g = \sum_{j=1}^k p_{g,j} \left(
#' \begin{array}{cc} \Sigma_{g,j}^{-1} & 0 \\ 0 & 0 \end{array}\right)}
#' Here, \eqn{\Sigma_{g,j}} is the leading \eqn{j\times j}
#' principal submatrix of \eqn{\Sigma_g}. It follows that
#' \deqn{\text{Var}(\hat{\mu}_{1,k} - \hat{\mu}_{2,k}) =
#' \frac{1}{n}\left(\frac{1}{\pi_1} J_1^{-1}[k,k] +
#' \frac{1}{\pi_2} J_2^{-1}[k,k]\right)}
#'
#' The observed data pattern probabilities depend on the accrual and
#' dropout distributions.
#' Let \eqn{H(u)} denote the cumulative distribution function of
#' enrollment time \eqn{u},
#' \eqn{G_g(t)} denote the survival function of dropout time \eqn{t} for
#' treatment group \eqn{g},
#' and \eqn{\tau} denote the calendar time at interim or final analysis.
#' Then, for \eqn{j=1,\ldots,k-1}, the probability that a subject in group
#' \eqn{g} falls into observed data pattern \eqn{j} is
#' \deqn{p_{g,j} = H(\tau - t_j)G_g(t_j) - H(\tau - t_{j+1})G_g(t_{j+1})}
#' For the last pattern (\eqn{j=k}, i.e., completers),
#' \deqn{p_{g,k} = H(\tau - t_k)G_g(t_k)}
#' For the final analysis, \eqn{\tau} is the study duration,
#' so \eqn{H(\tau - t_j) = 1} for all \eqn{j}.
#' Therefore, the pattern probabilities depend only on the dropout
#' distribution:
#' \deqn{p_{g,j} = G_g(t_j) - G_g(t_{j+1}), \quad j=1,\ldots,k-1}
#' and \deqn{p_{g,k} = G_g(t_k)}
#'
#' Cumulative dropout probabilities at post-baseline time points
#' can be used to define a piecewise exponential dropout distribution.
#' Let \eqn{F_g(t_j)} denote the cumulative probability of dropout
#' by time \eqn{t_j} for treatment group \eqn{g}.
#' The left endpoints of the piecewise survival time intervals are given by
#' \eqn{t_0=0,t_1,\ldots,t_{k-1}}.
#' The hazard rate in the interval \eqn{(t_{j-1},t_j]} is given by
#' \deqn{\gamma_{g,j} = -\log\left(\frac{1 - F_g(t_j)}{1 - F_g(t_{j-1})}
#' \right) / (t_j - t_{j-1}), \quad j=1,\ldots,k}
#'
#' @return An S3 class \code{designMeanDiffMMRM} object with three
#' components:
#'
#' * \code{overallResults}: A data frame containing the following variables:
#'
#'     - \code{overallReject}: The overall rejection probability.
#'
#'     - \code{alpha}: The overall significance level.
#'
#'     - \code{attainedAlpha}: The attained significance level, which is
#'       different from the overall significance level in the presence of
#'       futility stopping.
#'
#'     - \code{kMax}: The number of stages.
#'
#'     - \code{theta}: The parameter value.
#'
#'     - \code{information}: The maximum information.
#'
#'     - \code{expectedInformationH1}: The expected information under H1.
#'
#'     - \code{expectedInformationH0}: The expected information under H0.
#'
#'     - \code{drift}: The drift parameter, equal to
#'       \code{theta*sqrt(information)}.
#'
#'     - \code{inflationFactor}: The inflation factor (relative to the
#'       fixed design).
#'
#'     - \code{numberOfSubjects}: The maximum number of subjects.
#'
#'     - \code{studyDuration}: The maximum study duration.
#'
#'     - \code{expectedNumberOfSubjectsH1}: The expected number of subjects
#'       under H1.
#'
#'     - \code{expectedNumberOfSubjectsH0}: The expected number of subjects
#'       under H0.
#'
#'     - \code{expectedStudyDurationH1}: The expected study duration
#'       under H1.
#'
#'     - \code{expectedStudyDurationH0}: The expected study duration
#'       under H0.
#'
#'     - \code{accrualDuration}: The accrual duration.
#'
#'     - \code{followupTime}: The follow-up time.
#'
#'     - \code{fixedFollowup}: Whether a fixed follow-up design is used.
#'
#'     - \code{meanDiffH0}: The mean difference under H0.
#'
#'     - \code{meanDiff}: The mean difference under H1.
#'
#' * \code{byStageResults}: A data frame containing the following variables:
#'
#'     - \code{informationRates}: The information rates.
#'
#'     - \code{efficacyBounds}: The efficacy boundaries on the Z-scale.
#'
#'     - \code{futilityBounds}: The futility boundaries on the Z-scale.
#'
#'     - \code{rejectPerStage}: The probability for efficacy stopping.
#'
#'     - \code{futilityPerStage}: The probability for futility stopping.
#'
#'     - \code{cumulativeRejection}: The cumulative probability for efficacy
#'       stopping.
#'
#'     - \code{cumulativeFutility}: The cumulative probability for futility
#'       stopping.
#'
#'     - \code{cumulativeAlphaSpent}: The cumulative alpha spent.
#'
#'     - \code{efficacyP}: The efficacy boundaries on the p-value scale.
#'
#'     - \code{futilityP}: The futility boundaries on the p-value scale.
#'
#'     - \code{information}: The cumulative information.
#'
#'     - \code{efficacyStopping}: Whether to allow efficacy stopping.
#'
#'     - \code{futilityStopping}: Whether to allow futility stopping.
#'
#'     - \code{rejectPerStageH0}: The probability for efficacy stopping
#'       under H0.
#'
#'     - \code{futilityPerStageH0}: The probability for futility stopping
#'       under H0.
#'
#'     - \code{cumulativeRejectionH0}: The cumulative probability for
#'       efficacy stopping under H0.
#'
#'     - \code{cumulativeFutilityH0}: The cumulative probability for futility
#'       stopping under H0.
#'
#'     - \code{efficacyMeanDiff}: The efficacy boundaries on the mean
#'       difference scale.
#'
#'     - \code{futilityMeanDiff}: The futility boundaries on the mean
#'       difference scale.
#'
#'     - \code{numberOfSubjects}: The number of subjects.
#'
#'     - \code{numberOfCompleters}: The number of completers.
#'
#'     - \code{analysisTime}: The average time since trial start.
#'
#' * \code{settings}: A list containing the following input parameters:
#'
#'     - \code{typeAlphaSpending}: The type of alpha spending.
#'
#'     - \code{parameterAlphaSpending}: The parameter value for alpha
#'       spending.
#'
#'     - \code{userAlphaSpending}: The user defined alpha spending.
#'
#'     - \code{typeBetaSpending}: The type of beta spending.
#'
#'     - \code{parameterBetaSpending}: The parameter value for beta spending.
#'
#'     - \code{userBetaSpending}: The user defined beta spending.
#'
#'     - \code{spendingTime}: The error spending time at each analysis.
#'
#'     - \code{allocationRatioPlanned}: The allocation ratio for the active
#'       treatment versus control.
#'
#'     - \code{accrualTime}: A vector that specifies the starting time of
#'       piecewise Poisson enrollment time intervals.
#'
#'     - \code{accrualIntensity}: A vector of accrual intensities.
#'        One for each accrual time interval.
#'
#'     - \code{piecewiseSurvivalTime}: A vector that specifies the
#'       starting time of piecewise exponential survival time intervals.
#'
#'     - \code{gamma1}: The hazard rate for exponential dropout or
#'       a vector of hazard rates for piecewise exponential dropout
#'       for the active treatment group.
#'
#'     - \code{gamma2}: The hazard rate for exponential dropout or
#'       a vector of hazard rates for piecewise exponential dropout
#'       for the control group.
#'
#'     - \code{k}: The number of postbaseline time points.
#'
#'     - \code{t}: The postbaseline time points.
#'
#'     - \code{covar1}: The covariance matrix for the repeated measures
#'       given baseline for the active treatment group.
#'
#'     - \code{covar2}: The covariance matrix for the repeated measures
#'       given baseline for the control group.
#'
#'     - \code{normalApproximation}: The type of computation of the p-values.
#'       If \code{TRUE}, the variance is assumed to be known, otherwise
#'       the calculations are performed with the t distribution.
#'
#'     - \code{rounding}: Whether to round up sample size.
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @references
#'
#' Kaifeng Lu, Xiaohui Luo, and Pei-Yun Chen.
#' Sample size estimation for repeated measures analysis in randomized
#' clinical trials with missing data.
#' The International Journal of Biostatistics 2008; 14(1), Article 9.
#'
#' @examples
#'
#' # function to generate the AR(1) correlation matrix
#' ar1_cor <- function(n, corr) {
#'   exponent <- abs(matrix((1:n) - 1, n, n, byrow = TRUE) - ((1:n) - 1))
#'   corr^exponent
#' }
#'
#' (design1 = getDesignMeanDiffMMRM(
#'   beta = 0.2,
#'   meanDiffH0 = 0,
#'   meanDiff = 0.5,
#'   k = 4,
#'   t = c(1,2,3,4),
#'   covar1 = ar1_cor(4, 0.7),
#'   accrualIntensity = 10,
#'   gamma1 = 0.02634013,
#'   gamma2 = 0.02634013,
#'   accrualDuration = NA,
#'   allocationRatioPlanned = 1,
#'   kMax = 3,
#'   alpha = 0.025,
#'   typeAlphaSpending = "sfOF"))
#'
#' @export
#'
getDesignMeanDiffMMRM <- function(
    beta = NA_real_,
    meanDiffH0 = 0,
    meanDiff = 0.5,
    k = 1,
    t = NA_real_,
    covar1 = diag(k),
    covar2 = NA_real_,
    accrualTime = 0,
    accrualIntensity = NA_real_,
    piecewiseSurvivalTime = 0,
    gamma1 = 0,
    gamma2 = 0,
    accrualDuration = NA_real_,
    allocationRatioPlanned = 1,
    normalApproximation = TRUE,
    rounding = TRUE,
    kMax = 1L,
    informationRates = NA_real_,
    efficacyStopping = NA_integer_,
    futilityStopping = NA_integer_,
    criticalValues = NA_real_,
    alpha = 0.025,
    typeAlphaSpending = "sfOF",
    parameterAlphaSpending = NA_real_,
    userAlphaSpending = NA_real_,
    futilityBounds = NA_real_,
    typeBetaSpending = "none",
    parameterBetaSpending = NA_real_,
    userBetaSpending = NA_real_,
    spendingTime = NA_real_) {

  nintervals = length(piecewiseSurvivalTime)

  if (is.na(beta) && is.na(accrualDuration)) {
    stop("beta and accrualDuration cannot be both missing")
  }

  if (!is.na(beta) && (beta >= 1-alpha || beta < 0.0001)) {
    stop("beta must lie in [0.0001, 1-alpha)")
  }

  if (k <= 0 || k != round(k)) {
    stop("k must be a positive integer")
  }

  if (any(is.na(t))) {
    stop("t must be provided")
  } else if (length(t) != k) {
    stop("Invalid length for t")
  } else if (any(t <= 0)) {
    stop("Elements of t must be positive")
  } else if (k > 1 && any(diff(t) <= 0)) {
    stop("Elements of t must be increasing")
  }

  if (!all(eigen(covar1)$values > 0)) {
    stop("covar1 must be positive definite")
  }

  if (any(is.na(covar2))) {
    covar2 = covar1
  } else if (!all(eigen(covar2)$values > 0)) {
    stop("covar2 must be positive definite")
  }

  if (accrualTime[1] != 0) {
    stop("accrualTime must start with 0")
  }

  if (length(accrualTime) > 1 && any(diff(accrualTime) <= 0)) {
    stop("accrualTime should be increasing")
  }

  if (any(is.na(accrualIntensity))) {
    stop("accrualIntensity must be provided")
  }

  if (length(accrualTime) != length(accrualIntensity)) {
    stop("accrualTime and accrualIntensity must have the same length")
  }

  if (any(accrualIntensity < 0)) {
    stop("accrualIntensity must be non-negative")
  }

  if (piecewiseSurvivalTime[1] != 0) {
    stop("piecewiseSurvivalTime must start with 0")
  }

  if (nintervals > 1 && any(diff(piecewiseSurvivalTime) <= 0)) {
    stop("piecewiseSurvivalTime should be increasing")
  }

  if (any(gamma1 < 0)) {
    stop("gamma1 must be non-negative")
  }

  if (any(gamma2 < 0)) {
    stop("gamma2 must be non-negative")
  }

  if (length(gamma1) != 1 && length(gamma1) != nintervals) {
    stop("Invalid length for gamma1")
  }

  if (length(gamma2) != 1 && length(gamma2) != nintervals) {
    stop("Invalid length for gamma2")
  }

  if (length(gamma1) == 1) {
    gamma1 = rep(gamma1, nintervals)
  }

  if (length(gamma2) == 1) {
    gamma2 = rep(gamma2, nintervals)
  }

  if (!is.na(accrualDuration) && accrualDuration <= 0) {
    stop("accrualDuration must be positive")
  }

  if (allocationRatioPlanned <= 0) {
    stop("allocationRatioPlanned must be positive")
  }

  if (alpha < 0.00001 || alpha >= 1) {
    stop("alpha must lie in [0.00001, 1)")
  }

  if (any(is.na(informationRates))) {
    informationRates = (1:kMax)/kMax
  }


  r = allocationRatioPlanned/(1 + allocationRatioPlanned)

  directionUpper = meanDiff > meanDiffH0

  theta = ifelse(directionUpper, meanDiff - meanDiffH0,
                 meanDiffH0 - meanDiff)


  # function to obtain the info for mean difference
  f_info <- function(tau, k, t, covar1, covar2, accrualTime,
                     accrualIntensity, piecewiseSurvivalTime,
                     gamma1, gamma2, accrualDuration, r) {

    # total number of enrolled subjects at interim analysis
    n = accrual(tau, accrualTime, accrualIntensity, accrualDuration)

    # total number of enrolled subjects at each time point
    ns = accrual(tau - t, accrualTime, accrualIntensity, accrualDuration)

    # probability of not dropping out at each time point by treatment
    q1 = ptpwexp(t, piecewiseSurvivalTime, gamma1, lower.tail = FALSE)
    q2 = ptpwexp(t, piecewiseSurvivalTime, gamma2, lower.tail = FALSE)

    # number of subjects remaining at each time point by treatment
    m1 = r*ns*q1
    m2 = (1-r)*ns*q2

    # number of subjects dropping out at each time point by treatment
    n1 = c(-diff(m1), m1[k])
    n2 = c(-diff(m2), m2[k])

    # information matrix in each treatment group
    I1 = 0
    I2 = 0
    I = matrix(0, k, k)
    for (j in 1:k) {
      I[1:j, 1:j] = solve(covar1[1:j, 1:j])
      I1 = I1 + n1[j]*I
      I[1:j, 1:j] = solve(covar2[1:j, 1:j])
      I2 = I2 + n2[j]*I
    }

    # variance for the treatment mean at the last time point by treatment
    V1 = solve(I1)
    V2 = solve(I2)

    # information for treatment difference
    IMax = 1/(V1[k,k] + V2[k,k])

    # variance for treatment difference for one sampling unit
    v1 = (V1[k,k] + V2[k,k])*n

    # inflation factors for each treatment
    phi1 = V1[k,k]/(covar1[k,k]/(r*n))
    phi2 = V2[k,k]/(covar2[k,k]/((1-r)*n))

    # information for treatment mean difference at the last time point
    list(IMax = IMax, v1 = v1, phi1 = phi1, phi2 = phi2)
  }


  # power calculation
  if (is.na(beta)) {
    n = accrual(accrualDuration, accrualTime, accrualIntensity,
                accrualDuration)

    if (rounding) {
      n = ceiling(n - 1.0e-12)
      accrualDuration = getAccrualDurationFromN(n, accrualTime,
                                                accrualIntensity)
    }

    studyDuration = accrualDuration + t[k]
    out = f_info(studyDuration, k, t, covar1, covar2, accrualTime,
                 accrualIntensity, piecewiseSurvivalTime,
                 gamma1, gamma2, accrualDuration, r)
    IMax = out$IMax
    phi1 = out$phi1
    phi2 = out$phi2

    des = getDesign(
      beta = NA, IMax, theta,
      kMax, informationRates,
      efficacyStopping, futilityStopping,
      criticalValues, alpha, typeAlphaSpending,
      parameterAlphaSpending, userAlphaSpending,
      futilityBounds, typeBetaSpending,
      parameterBetaSpending, userBetaSpending,
      spendingTime, 1)

    if (kMax == 1 && !normalApproximation) { # t-test for fixed design
      if (!any(is.na(criticalValues))) {
        alpha = 1 - pnorm(criticalValues)
      }

      nu = n*r/phi1 + n*(1-r)/phi2 - 2
      b = qt(1-alpha, nu)
      ncp = theta*sqrt(IMax)
      power = pt(b, nu, ncp, lower.tail = FALSE)

      des$overallResults$overallReject = power
      des$byStageResults$rejectPerStage = power
      des$byStageResults$futilityPerStage = 1 - power
      des$byStageResults$cumulativeRejection = power
      des$byStageResults$cumulativeFutility = 1 - power
      des$byStageResults$efficacyBounds = b
      des$byStageResults$futilityBounds = b

      delta = b/sqrt(IMax)
      if (directionUpper) {
        des$byStageResults$efficacyMeanDiff = delta + meanDiffH0
        des$byStageResults$futilityMeanDiff = delta + meanDiffH0
      } else {
        des$byStageResults$efficacyMeanDiff = -delta + meanDiffH0
        des$byStageResults$futilityMeanDiff = -delta + meanDiffH0
      }
    } else {
      if (directionUpper) {
        des$byStageResults$efficacyMeanDiff =
          des$byStageResults$efficacyTheta + meanDiffH0
        des$byStageResults$futilityMeanDiff =
          des$byStageResults$futilityTheta + meanDiffH0
      } else {
        des$byStageResults$efficacyMeanDiff =
          -des$byStageResults$efficacyTheta + meanDiffH0
        des$byStageResults$futilityMeanDiff =
          -des$byStageResults$futilityTheta + meanDiffH0
      }
    }
  } else { # sample size calculation
    des = getDesign(
      beta, IMax = NA, theta,
      kMax, informationRates,
      efficacyStopping, futilityStopping,
      criticalValues, alpha, typeAlphaSpending,
      parameterAlphaSpending, userAlphaSpending,
      futilityBounds, typeBetaSpending,
      parameterBetaSpending, userBetaSpending,
      spendingTime, 1)

    IMax = des$overallResults$information

    # v1, phi1, phi2 for the last time point are not affected by the
    # accrual duration if everyone has a chance to complete the study
    out = f_info(1 + t[k], k, t, covar1, covar2, accrualTime,
                 accrualIntensity, piecewiseSurvivalTime,
                 gamma1, gamma2, 1, r)
    v1 = out$v1
    phi1 = out$phi1
    phi2 = out$phi2

    n = IMax*v1

    if (kMax == 1 && !normalApproximation) { # t-test for fixed design
      if (!any(is.na(criticalValues))) {
        alpha = 1 - pnorm(criticalValues)
      }

      # power for t-test
      f = function(n, r, phi1, phi2, alpha, theta, v1) {
        nu = n*r/phi1 + n*(1-r)/phi2 - 2
        b = qt(1-alpha, nu)
        ncp = theta*sqrt(n/v1)
        pt(b, nu, ncp, lower.tail = FALSE)
      }

      n = uniroot(function(n) {
        f(n, r, phi1, phi2, alpha, theta, v1) - (1-beta)
      }, c(0.5*n, 1.5*n))$root

      if (rounding) {
        n = ceiling(n - 1.0e-12)
      }

      if (is.na(accrualDuration)) {
        accrualDuration = getAccrualDurationFromN(n, accrualTime,
                                                  accrualIntensity)
      } else { # calculate accrualIntensity
        accrualIntensity = uniroot(function(x) {
          accrual(accrualDuration, accrualTime, x*accrualIntensity,
                  accrualDuration) - n
        }, c(0.001, 240))$root*accrualIntensity
      }

      # final maximum information
      IMax = n/v1
      power = f(n, r, phi1, phi2, alpha, theta, v1)

      des$overallResults$overallReject = power
      des$byStageResults$rejectPerStage = power
      des$byStageResults$futilityPerStage = 1 - power
      des$byStageResults$cumulativeRejection = power
      des$byStageResults$cumulativeFutility = 1 - power
      des$byStageResults$efficacyBounds = b
      des$byStageResults$futilityBounds = b

      delta = b/sqrt(IMax)
      if (directionUpper) {
        des$byStageResults$efficacyMeanDiff = delta + meanDiffH0
        des$byStageResults$futilityMeanDiff = delta + meanDiffH0
      } else {
        des$byStageResults$efficacyMeanDiff = -delta + meanDiffH0
        des$byStageResults$futilityMeanDiff = -delta + meanDiffH0
      }
    } else {
      if (rounding) {
        n = ceiling(n - 1.0e-12)
      }

      if (is.na(accrualDuration)) {
        accrualDuration = getAccrualDurationFromN(n, accrualTime,
                                                  accrualIntensity)
      } else { # calculate accrualIntensity
        accrualIntensity = uniroot(function(x) {
          accrual(accrualDuration, accrualTime, x*accrualIntensity,
                  accrualDuration) - n
        }, c(0.001, 240))$root*accrualIntensity
      }

      # final maximum information
      IMax = n/v1

      des = getDesign(
        beta = NA, IMax, theta,
        kMax, informationRates,
        efficacyStopping, futilityStopping,
        criticalValues, alpha, typeAlphaSpending,
        parameterAlphaSpending, userAlphaSpending,
        futilityBounds, typeBetaSpending,
        parameterBetaSpending, userBetaSpending,
        spendingTime, 1)

      if (directionUpper) {
        des$byStageResults$efficacyMeanDiff =
          des$byStageResults$efficacyTheta + meanDiffH0
        des$byStageResults$futilityMeanDiff =
          des$byStageResults$futilityTheta + meanDiffH0
      } else {
        des$byStageResults$efficacyMeanDiff =
          -des$byStageResults$efficacyTheta + meanDiffH0
        des$byStageResults$futilityMeanDiff =
          -des$byStageResults$futilityTheta + meanDiffH0
      }
    }
  }


  # timing of interim analysis
  studyDuration = accrualDuration + t[k]
  information = IMax*informationRates
  analysisTime = rep(0, kMax)

  if (kMax > 1) {
    for (i in 1:(kMax-1)) {
      analysisTime[i] = uniroot(function(tau) {
        out = f_info(tau, k, t, covar1, covar2, accrualTime,
                     accrualIntensity, piecewiseSurvivalTime,
                     gamma1, gamma2, accrualDuration, r)
        out$IMax - information[i]
      }, c(t[k]+0.001, studyDuration))$root # at least 1 completer
    }
  }

  analysisTime[kMax] = studyDuration

  numberOfSubjects = accrual(analysisTime, accrualTime, accrualIntensity,
                             accrualDuration)
  numberOfCompleters = accrual(analysisTime - t[k], accrualTime,
                               accrualIntensity, accrualDuration)

  p = des$byStageResults$rejectPerStage +
    des$byStageResults$futilityPerStage

  pH0 = des$byStageResults$rejectPerStageH0 +
    des$byStageResults$futilityPerStageH0

  des$overallResults$numberOfSubjects = n
  des$overallResults$studyDuration = studyDuration
  des$overallResults$expectedNumberOfSubjectsH1 = sum(p*numberOfSubjects)
  des$overallResults$expectedNumberOfSubjectsH0 = sum(pH0*numberOfSubjects)
  des$overallResults$expectedStudyDurationH1 = sum(p*analysisTime)
  des$overallResults$expectedStudyDurationH0 = sum(pH0*analysisTime)
  des$overallResults$accrualDuration = accrualDuration
  des$overallResults$followupTime = t[k]
  des$overallResults$fixedFollowup = TRUE
  des$overallResults$meanDiffH0 = meanDiffH0
  des$overallResults$meanDiff = meanDiff

  des$byStageResults$numberOfSubjects = numberOfSubjects
  des$byStageResults$numberOfCompleters = numberOfCompleters
  des$byStageResults$analysisTime = analysisTime
  des$byStageResults$efficacyTheta = NULL
  des$byStageResults$futilityTheta = NULL

  des$settings$allocationRatioPlanned = allocationRatioPlanned
  des$settings$accrualTime = accrualTime
  des$settings$accrualIntensity = accrualIntensity
  des$settings$piecewiseSurvivalTime = piecewiseSurvivalTime
  des$settings$gamma1 = gamma1
  des$settings$gamma2 = gamma2
  des$settings$k = k
  des$settings$t = t
  des$settings$covar1 = covar1
  des$settings$covar2 = covar2
  des$settings$normalApproximation = normalApproximation
  des$settings$rounding = rounding
  des$settings$varianceRatio = NULL

  attr(des, "class") = "designMeanDiffMMRM"

  des
}


#' @title Power and Sample Size for Direct Treatment Effects in Crossover
#' Trials
#' @description Obtains the power and sample size for direct treatment
#' effects in crossover trials accounting or without accounting for
#' carryover effects.
#'
#' @param beta The type II error.
#' @param n The total sample size.
#' @param trtpair The treatment pair of interest to power the study.
#'   If not given, it defaults to comparing the first treatment to the
#'   last treatment.
#' @param carryover Whether to account for carryover effects in the
#'   power calculation. Defaults to TRUE.
#' @param meanDiffH0 The mean difference for the treatment pair of interest
#'   under the null hypothesis. Defaults to 0.
#' @param meanDiff The mean difference for the treatment pair of interest
#'   under the alternative hypothesis.
#' @param stDev The standard deviation for within-subject random error.
#' @param corr The intra-subject correlation due to subject random effect.
#' @param design The crossover design represented by a matrix with
#'   rows indexing the sequences, columns indexing the periods, and
#'   matrix entries indicating the treatments.
#' @param cumdrop The cumulative dropout rate over periods.
#' @param allocationRatioPlanned Allocation ratio for the sequences.
#'   Defaults to equal randomization if not provided.
#' @param normalApproximation The type of computation of the p-values.
#'   If \code{TRUE}, the variance is assumed to be known, otherwise
#'   the calculations are performed with the t distribution.
#' @param rounding Whether to round up the sample size. Defaults to TRUE for
#'   sample size rounding.
#' @param alpha The one-sided significance level. Defaults to 0.025.
#'
#' @details
#' The linear mixed-effects model to assess the direct treatment effects
#' in the presence of carryover treatment effects is given by
#' \deqn{y_{ijk} = \mu + \alpha_i + b_{ij} + \gamma_k + \tau_{d(i,k)}
#' + \lambda_{c(i,k-1)} + e_{ijk}}
#' \deqn{i=1,\ldots,n; j=1,\ldots,r_i; k = 1,\ldots,p; d,c = 1,\ldots,t}
#' where \eqn{\mu} is the general mean, \eqn{\alpha_i} is the effect of
#' the \eqn{i}th treatment sequence, \eqn{b_{ij}} is the random effect
#' with variance \eqn{\sigma_b^2} for the \eqn{j}th subject of the
#' \eqn{i}th treatment sequence, \eqn{\gamma_k} is the period effect,
#' and \eqn{e_{ijk}} is the random error with variance \eqn{\sigma^2}
#' for the subject in period \eqn{k}. The direct effect of the treatment
#' administered in period \eqn{k} of sequence \eqn{i} is
#' \eqn{\tau_{d(i,k)}}, and \eqn{\lambda_{c(i,k-1)}} is the carryover
#' effect of the treatment administered in period \eqn{k-1} of sequence
#' \eqn{i}. The value of the carryover effect for the observed
#' response in the first period is \eqn{\lambda_{c(i,0)} = 0} since
#' there is no carryover effect in the first period. The intra-subject
#' correlation due to the subject random effect is
#' \deqn{\rho = \frac{\sigma_b^2}{\sigma_b^2 + \sigma^2}.}
#' Therefore, \code{stDev} = \eqn{\sigma^2} and \code{corr} = \eqn{\rho}.
#' By constructing the design matrix \eqn{X} for the linear model with
#' a compound symmetry covariance matrix for the response vector of
#' a subject, we can obtain \deqn{Var(\hat{\beta}) = (X'V^{-1}X)^{-1}.}
#'
#' The covariance matrix for the direct treatment effects and
#' carryover treatment effects can be extracted from the appropriate
#' sub-matrices. The covariance matrix for the direct treatment effects
#' without accounting for the carryover treatment effects can be obtained
#' by omitting the carryover effect terms from the model.
#'
#' The power is for the direct treatment effect for the treatment pair of
#' interest with or without accounting for carryover effects as determined
#' by the input parameter \code{carryover}. The relative efficiency is
#' for the direct treatment effect for the treatment pair of interest
#' accounting for carryover effects relative to that without accounting
#' for carryover effects.
#'
#' The degrees of freedom for the t-test accounting for carryover effects
#' can be calculated as the total number of observations minus
#' the number of subjects minus \eqn{p-1} minus \eqn{2(t-1)}
#' to account for the subject effect, period effect, and direct and
#' carryover treatment effects. The degrees of freedom for the t-test
#' without accounting for carryover effects is equal to the total number
#' of observations minus the number of subjects minus \eqn{p-1} minus
#' \eqn{t-1} to account for the subject effect, period effect, and direct
#' treatment effects.
#'
#' @return An S3 class \code{designMeanDiffCarryover} object with the
#' following components:
#'
#' * \code{power}: The power to reject the null hypothesis.
#'
#' * \code{alpha}: The one-sided significance level.
#'
#' * \code{numberOfSubjects}: The maximum number of subjects.
#'
#' * \code{trtpair}: The treatment pair of interest to power the study.
#'
#' * \code{carryover}: Whether to account for carryover effects in
#'   the power calculation.
#'
#' * \code{meanDiffH0}: The mean difference for the treatment pair of
#'   interest under the null hypothesis.
#'
#' * \code{meanDiff}: The mean difference for the treatment pair of
#'   interest under the alternative hypothesis.
#'
#' * \code{stDev}: The standard deviation for within-subject random error.
#'
#' * \code{corr}: The intra-subject correlation due to subject random effect.
#'
#' * \code{design}: The crossover design represented by a matrix with
#'   rows indexing the sequences, columns indexing the periods, and
#'   matrix entries indicating the treatments.
#'
#' * \code{designMatrix}: The design matrix accounting for intercept,
#'   sequence, period, direct treatment effects and carryover treatment
#'   effects when \code{carryover = TRUE}, or the design matrix
#'   accounting for intercept, sequence, period, and direct treatment
#'   effects when \code{carryover = FALSE}.
#'
#' * \code{nseq}: The number of sequences.
#'
#' * \code{nprd}: The number of periods.
#'
#' * \code{ntrt}: The number of treatments.
#'
#' * \code{cumdrop}: The cumulative dropout rate over periods.
#'
#' * \code{V_direct_only}: The covariance matrix for direct treatment
#'   effects without accounting for carryover effects. The treatment
#'   comparisons for the covariance matrix are for the first \eqn{t-1}
#'   treatments relative to the last treatment.
#'
#' * \code{V_direct_carry}: The covariance matrix for direct and
#'   carryover treatment effects.
#'
#' * \code{v_direct_only}: The variance of the direct treatment effect for
#'   the treatment pair of interest without accounting for carryover
#'   effects.
#'
#' * \code{v_direct}: The variance of the direct treatment effect for
#'   the treatment pair of interest accounting for carryover effects.
#'
#' * \code{v_carry}: The variance of the carryover treatment effect for
#'   the treatment pair of interest.
#'
#' * \code{releff_direct}: The relative efficiency of the design
#'   for estimating the direct treatment effect for the treatment pair
#'   of interest after accounting for carryover effects with respect to
#'   that without accounting for carryover effects. This is equal to
#'   \code{v_direct_only/v_direct}.
#'
#' * \code{releff_carry}: The relative efficiency of the design
#'   for estimating the carryover effect for the treatment pair
#'   of interest. This is equal to \code{v_direct_only/v_carry}.
#'
#' * \code{half_width}: The half-width of the confidence interval for the
#'   direct treatment effect for the treatment pair of interest.
#'
#' * \code{nu}: Degrees of freedom for the t-test.
#'
#' * \code{allocationRatioPlanned}: Allocation ratio for the sequences.
#'
#' * \code{normalApproximation}: The type of computation of the p-values.
#'   If \code{TRUE}, the variance is assumed to be known, otherwise
#'   the calculations are performed with the t distribution.
#'
#' * \code{rounding}: Whether to round up the sample size.
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @references
#'
#' Robert O. Kuehl. Design of Experiments: Statistical Principles of
#' Research Design and Analysis. Brooks/Cole: Pacific Grove, CA. 2000.
#'
#' @examples
#'
#' # Williams design for 4 treatments
#'
#' (design1 = getDesignMeanDiffCarryover(
#'   beta = 0.2, n = NA,
#'   meanDiff = 0.5, stDev = 1,
#'   design = matrix(c(1, 4, 2, 3,
#'                     2, 1, 3, 4,
#'                     3, 2, 4, 1,
#'                     4, 3, 1, 2),
#'                   4, 4, byrow = TRUE),
#'   alpha = 0.025))
#'
#' @export
#'
getDesignMeanDiffCarryover <- function(
    beta = NA_real_,
    n = NA_real_,
    trtpair = NA_real_,
    carryover = TRUE,
    meanDiffH0 = 0,
    meanDiff = 0.5,
    stDev = 1,
    corr = 0.5,
    design = NA_real_,
    cumdrop = NA_real_,
    allocationRatioPlanned = NA_real_,
    normalApproximation = FALSE,
    rounding = TRUE,
    alpha = 0.025) {

  if (is.na(beta) && is.na(n)) {
    stop("beta and n cannot be both missing")
  }

  if (!is.na(beta) && !is.na(n)) {
    stop("Only one of beta and n should be provided")
  }

  if (!is.na(beta) && (beta >= 1-alpha || beta < 0.0001)) {
    stop("beta must lie in [0.0001, 1-alpha)")
  }

  if (!is.na(n) && n <= 0) {
    stop("n must be positive")
  }

  if (stDev <= 0) {
    stop("stDev must be positive")
  }

  if (corr <= -1 || corr >= 1) {
    stop("corr must lie between -1 and 1")
  }

  if (any(is.na(design))) {
    stop("design must be provided")
  }

  nseq = nrow(design)
  nprd = ncol(design)
  ntrt = length(unique(c(design)))

  if (any(design <= 0 | design != round(design)) ||
      !all.equal(sort(unique(c(design))), 1:ntrt)) {
    stop(paste("Elements of design must be positive integers",
               "ranging from 1 to", ntrt))
  }

  if (any(is.na(trtpair))) {
    trtpair = c(1, ntrt)
  } else if (length(trtpair) != 2 ||
             trtpair[1] == trtpair[2] ||
             !all(trtpair %in% 1:ntrt)) {
    stop(paste("trtpair must be a vector of two distinct integers",
               "selected from 1 to", ntrt))
  }

  # number of model parameters consisting of the intercept, sequence effects,
  # period effects, direct treatment effects, and carryover treatment effects
  q = 1 + (nseq-1) + (nprd-1) + (ntrt-1) + (ntrt-1)

  # number of model parameters consisting of the intercept, sequence effects,
  # period effects, and direct treatment effects
  q0 = 1 + (nseq-1) + (nprd-1) + (ntrt-1)

  # offset of direct treatment effect
  l = 1 + (nseq-1) + (nprd-1)

  # start of direct treatment effect
  l1 = l + 1

  # end of direct treatment effect
  m = 1 + (nseq-1) + (nprd-1) + (ntrt-1)

  if (any(is.na(cumdrop))) {
    cumdrop = rep(0, nprd)
  } else if (length(cumdrop) != nprd) {
    stop(paste("cumdrop should have", nprd, "elements"))
  }

  if (any(cumdrop < 0)) {
    stop("Elements of cumdrop must be nonnegative")
  } else if (nprd > 1 && any(diff(cumdrop) < 0)) {
    stop("Elements of cumdrop must be nondecreasing")
  } else if (any(cumdrop >= 1)) {
    stop("Elements of cumdrop must be less than 1")
  }

  # observed data pattern probabilities
  p = diff(c(cumdrop, 1))


  if (any(is.na(allocationRatioPlanned))) {
    allocationRatioPlanned = rep(1, nrow(design))
  }

  if (length(allocationRatioPlanned) != nseq-1 &&
      length(allocationRatioPlanned) != nseq) {
    stop(paste("allocationRatioPlanned should have", nseq-1,
               "or", nseq, "elements"))
  }

  if (length(allocationRatioPlanned) == nseq-1) {
    allocationRatioPlanned = c(allocationRatioPlanned, 1)
  }

  if (any(allocationRatioPlanned <= 0)) {
    stop("Elements of allocationRatioPlanned must be positive")
  }

  if (alpha < 0.00001 || alpha >= 1) {
    stop("alpha must lie in [0.00001, 1)")
  }


  # treatment sequence randomization probabilities
  r = allocationRatioPlanned/sum(allocationRatioPlanned)

  directionUpper = meanDiff > meanDiffH0

  theta = ifelse(directionUpper, meanDiff - meanDiffH0,
                 meanDiffH0 - meanDiff)

  # compound symmetry covariance matrix for repeated measures
  Sigma = stDev^2/(1-corr) * ((1-corr)*diag(nprd) + corr)

  # model design matrix with carryover effect
  X = matrix(0, nseq*nprd, q)
  X[,1] = 1
  for (i in 1:nseq) {
    for (j in 1:nprd) {
      k = (i-1)*nprd + j

      # sequence effects
      if (i < nseq) {
        X[k, 1 + i] = 1
      }

      # period effects
      if (j < nprd) {
        X[k, 1 + (nseq-1) + j] = 1
      }

      # direct treatment effects
      if (design[i,j] < ntrt) {
        X[k, 1 + (nseq-1) + (nprd-1) + design[i,j]] = 1
      }

      # carryover treatment effects
      if (j > 1 && design[i,j-1] < ntrt) {
        X[k, 1 + (nseq-1) + (nprd-1) + (ntrt-1) + design[i,j-1]] = 1
      }
    }
  }

  # singular values of X
  d = svdcpp(X, outtransform = FALSE, decreasing = FALSE)$d
  tol = max(dim(X))*.Machine$double.eps
  if (carryover && sum(d >= tol * max(d)) < q) {
    stop("The crossover design is overparameterized for carryover effects")
  }

  if (carryover) {
    # information matrix for model parameters with carryover effects
    I = 0
    for (i in 1:nseq) {
      offset = (i-1)*nprd
      J = 0
      for (j in 1:nprd) {
        idx = offset + (1:j)

        if (j==1) {
          x = t(as.matrix(X[idx,]))
        } else {
          x = X[idx,]
        }

        J = J + p[j]*t(x) %*% solve(Sigma[1:j,1:j]) %*% x
      }
      I = I + r[i]*J
    }

    # covariance matrix for model parameters with carryover effects
    V = solve(I)

    # variance for direct treatment effect for one sampling unit with carryover
    if (trtpair[2]==ntrt) {
      v1 = V[l+trtpair[1], l+trtpair[1]]
    } else if (trtpair[1]==ntrt) {
      v1 = V[l+trtpair[2], l+trtpair[2]]
    } else {
      v1 = V[l+trtpair[1], l+trtpair[1]] +
        V[l+trtpair[2], l+trtpair[2]] -
        2*V[l+trtpair[1], l+trtpair[2]]
    }


    # variance for carryover treatment effect for one sampling unit
    if (trtpair[2]==ntrt) {
      v2 = V[m+trtpair[1], m+trtpair[1]]
    } else if (trtpair[1]==ntrt) {
      v2 = V[m+trtpair[2], m+trtpair[2]]
    } else {
      v2 = V[m+trtpair[1], m+trtpair[1]] +
        V[m+trtpair[2], m+trtpair[2]] -
        2*V[m+trtpair[1], m+trtpair[2]]
    }
  }


  # design matrix for model parameters without carryover effects
  X0 = X[, 1:m]
  if (!carryover) {
    # singular values of X0
    d0 = svdcpp(X0, outtransform = FALSE, decreasing = FALSE)$d
    tol0 = max(dim(X0))*.Machine$double.eps
    if (sum(d0 >= tol0 * max(d0)) < q0) {
      stop("The crossover design is overparameterized")
    }
  }

  # information matrix for model parameters without carryover effects
  I0 = 0
  for (i in 1:nseq) {
    offset = (i-1)*nprd
    J = 0
    for (j in 1:nprd) {
      idx = offset+(1:j)

      if (j==1) {
        x = t(as.matrix(X0[idx,]))
      } else {
        x = X0[idx,]
      }

      J = J + p[j]*t(x) %*% solve(Sigma[1:j,1:j]) %*% x
    }
    I0 = I0 + r[i]*J
  }

  # covariance matrix for model parameters without carryover effects
  V0 = solve(I0)

  # variance for direct treatment effect for one sampling unit w/o carryover
  if (trtpair[2]==ntrt) {
    v0 = V0[l+trtpair[1], l+trtpair[1]]
  } else if (trtpair[1]==ntrt) {
    v0 = V0[l+trtpair[2], l+trtpair[2]]
  } else {
    v0 = V0[l+trtpair[1], l+trtpair[1]] +
      V0[l+trtpair[2], l+trtpair[2]] -
      2*V0[l+trtpair[1], l+trtpair[2]]
  }

  if (carryover) {
    v = v1
  } else {
    v = v0
  }

  # power for t test
  f = function(n) {
    # residual degrees of freedom after accounting for the subject effects
    mean_nprd = sum(p*(1:nprd))
    nu = n*mean_nprd - 1 - (n-1) - (nprd-1) - (ntrt-1) - (ntrt-1)*carryover
    b = qt(1-alpha, nu)
    ncp = theta*sqrt(n/v)
    pt(b, nu, ncp, lower.tail = FALSE)
  }


  if (is.na(n)) { # calculate the sample size
    if (normalApproximation) {
      n = (qnorm(1-alpha) + qnorm(1-beta))^2*v/theta^2
    } else {
      n0 = (qnorm(1-alpha) + qnorm(1-beta))^2*v/theta^2
      n = uniroot(function(n) f(n) - (1-beta), c(0.5*n0, 1.5*n0))$root
    }
  }

  if (rounding) {
    n = ceiling(n - 1.0e-12)
  }

  if (normalApproximation) {
    power = pnorm(theta*sqrt(n/v) - qnorm(1-alpha))
  } else {
    power = f(n)
  }

  rownames(design) = paste0("Seq", 1:nseq)
  colnames(design) = paste0("Prd", 1:nprd)

  if (!normalApproximation) {
    mean_nprd = sum(p*(1:nprd))
    nu = n*mean_nprd - 1 - (n-1) - (nprd-1) - (ntrt-1) - (ntrt-1)*carryover
    hw = qt(1-alpha, nu)*sqrt(v/n)
  } else {
    nu = Inf
    hw = qnorm(1-alpha)*sqrt(v/n)
  }

  if (carryover) {
    des = list(
      power = power, alpha = alpha, numberOfSubjects = n,
      trtpair = trtpair, carryover = carryover,
      meanDiffH0 = meanDiffH0, meanDiff = meanDiff,
      stDev = stDev, corr = corr,
      design = design, designMatrix = X,
      nseq = nseq, nprd = nprd, ntrt = ntrt,
      cumdrop = cumdrop,
      V_direct_only = V0[l1:m,l1:m]/n,
      V_direct_carry = V[l1:q, l1:q]/n,
      v_direct_only = v0/n,
      v_direct = v1/n,
      v_carry = v2/n,
      releff_direct = v0/v1,
      releff_carry = v0/v2,
      half_width = hw, nu = nu,
      allocationRatioPlanned = allocationRatioPlanned,
      normalApproximation = normalApproximation,
      rounding = rounding)
  } else {
    des = list(
      power = power, alpha = alpha, numberOfSubjects = n,
      trtpair = trtpair, carryover = carryover,
      meanDiffH0 = meanDiffH0, meanDiff = meanDiff,
      stDev = stDev, corr = corr,
      design = design, designMatrix = X0,
      nseq = nseq, nprd = nprd, ntrt = ntrt,
      cumdrop = cumdrop,
      V_direct_only = V0[l1:m,l1:m]/n,
      v_direct_only = v0/n,
      half_width = hw, nu = nu,
      allocationRatioPlanned = allocationRatioPlanned,
      normalApproximation = normalApproximation,
      rounding = rounding)
  }

  attr(des, "class") = "designMeanDiffCarryover"

  des
}


#' @title Power and Sample Size for Equivalence in Direct Treatment Effects
#' in Crossover Trials
#' @description Obtains the power and sample size for equivalence in
#' direct treatment effects in crossover trials accounting or without
#' accounting for carryover effects.
#'
#' @param beta The type II error.
#' @param n The total sample size.
#' @param trtpair The treatment pair of interest to power the study.
#'   If not given, it defaults to comparing the first treatment to the
#'   last treatment.
#' @param carryover Whether to account for carryover effects in the
#'   power calculation. Defaults to TRUE.
#' @param meanDiffLower The lower equivalence limit of mean difference
#'   for the treatment pair of interest.
#' @param meanDiffUpper The upper equivalence limit of mean difference
#'   for the treatment pair of interest.
#' @param meanDiff The mean difference under the alternative hypothesis,
#' @param stDev The standard deviation for within-subject random error.
#' @param corr The intra-subject correlation due to subject random effect.
#' @param design The crossover design represented by a matrix with
#'   rows indexing the sequences, columns indexing the periods, and
#'   matrix entries indicating the treatments.
#' @param cumdrop The cumulative dropout rate over periods.
#' @param allocationRatioPlanned Allocation ratio for the sequences.
#'   Defaults to equal randomization if not provided.
#' @param normalApproximation The type of computation of the p-values.
#'   If \code{TRUE}, the variance is assumed to be known, otherwise
#'   the calculations are performed with the t distribution.
#' @param rounding Whether to round up the sample size. Defaults to TRUE for
#'   sample size rounding.
#' @param alpha The one-sided significance level. Defaults to 0.025.
#'
#' @details
#' The linear mixed-effects model to assess the direct treatment effects
#' in the presence of carryover treatment effects is given by
#' \deqn{y_{ijk} = \mu + \alpha_i + b_{ij} + \gamma_k + \tau_{d(i,k)}
#' + \lambda_{c(i,k-1)} + e_{ijk}}
#' \deqn{i=1,\ldots,n; j=1,\ldots,r_i; k = 1,\ldots,p; d,c = 1,\ldots,t}
#' where \eqn{\mu} is the general mean, \eqn{\alpha_i} is the effect of
#' the \eqn{i}th treatment sequence, \eqn{b_{ij}} is the random effect
#' with variance \eqn{\sigma_b^2} for the \eqn{j}th subject of the
#' \eqn{i}th treatment sequence, \eqn{\gamma_k} is the period effect,
#' and \eqn{e_{ijk}} is the random error with variance \eqn{\sigma^2}
#' for the subject in period \eqn{k}. The direct effect of the treatment
#' administered in period \eqn{k} of sequence \eqn{i} is
#' \eqn{\tau_{d(i,k)}}, and \eqn{\lambda_{c(i,k-1)}} is the carryover
#' effect of the treatment administered in period \eqn{k-1} of sequence
#' \eqn{i}. The value of the carryover effect for the observed
#' response in the first period is \eqn{\lambda_{c(i,0)} = 0} since
#' there is no carryover effect in the first period. The intra-subject
#' correlation due to the subject random effect is
#' \deqn{\rho = \frac{\sigma_b^2}{\sigma_b^2 + \sigma^2}.}
#' Therefore, \code{stDev} = \eqn{\sigma^2} and \code{corr} = \eqn{\rho}.
#' By constructing the design matrix \eqn{X} for the linear model with
#' a compound symmetry covariance matrix for the response vector of
#' a subject, we can obtain \deqn{Var(\hat{\beta}) = (X'V^{-1}X)^{-1}.}
#'
#' The covariance matrix for the direct treatment effects and
#' carryover treatment effects can be extracted from the appropriate
#' sub-matrices. The covariance matrix for the direct treatment effects
#' without accounting for the carryover treatment effects can be obtained
#' by omitting the carryover effect terms from the model.
#'
#' The power is for the direct treatment effect for the treatment pair of
#' interest with or without accounting for carryover effects as determined
#' by the input parameter \code{carryover}. The relative efficiency is
#' for the direct treatment effect for the treatment pair of interest
#' accounting for carryover effects relative to that without accounting
#' for carryover effects.
#'
#' The degrees of freedom for the t-test accounting for carryover effects
#' can be calculated as the total number of observations minus
#' the number of subjects minus \eqn{p-1} minus \eqn{2(t-1)}
#' to account for the subject effect, period effect, and direct and
#' carryover treatment effects. The degrees of freedom for the t-test
#' without accounting for carryover effects is equal to the total number
#' of observations minus the number of subjects minus \eqn{p-1} minus
#' \eqn{t-1} to account for the subject effect, period effect, and direct
#' treatment effects.
#'
#' @return An S3 class \code{designMeanDiffCarryover} object with the
#' following components:
#'
#' * \code{power}: The power to reject the null hypothesis.
#'
#' * \code{alpha}: The one-sided significance level.
#'
#' * \code{numberOfSubjects}: The maximum number of subjects.
#'
#' * \code{trtpair}: The treatment pair of interest to power the study.
#'
#' * \code{carryover}: Whether to account for carryover effects in
#'   the power calculation.
#'
#' * \code{meanDiffLower}: The lower equivalence limit of mean difference
#'   for the treatment pair of interest.
#'
#' * \code{meanDiffUpper}: The upper equivalence limit of mean difference
#'   for the treatment pair of interest.
#'
#' * \code{meanDiff}: The mean difference for the treatment pair of interest
#'   under the alternative hypothesis.
#'
#' * \code{stDev}: The standard deviation for within-subject random error.
#'
#' * \code{corr}: The intra-subject correlation due to subject random effect.
#'
#' * \code{design}: The crossover design represented by a matrix with
#'   rows indexing the sequences, columns indexing the periods, and
#'   matrix entries indicating the treatments.
#'
#' * \code{designMatrix}: The design matrix accounting for intercept,
#'   sequence, period, direct treatment effects and carryover treatment
#'   effects when \code{carryover = TRUE}, or the design matrix
#'   accounting for intercept, sequence, period, and direct treatment
#'   effects when \code{carryover = FALSE}.
#'
#' * \code{nseq}: The number of sequences.
#'
#' * \code{nprd}: The number of periods.
#'
#' * \code{ntrt}: The number of treatments.
#'
#' * \code{cumdrop}: The cumulative dropout rate over periods.
#'
#' * \code{V_direct_only}: The covariance matrix for direct treatment
#'   effects without accounting for carryover effects. The treatment
#'   comparisons for the covariance matrix are for the first \eqn{t-1}
#'   treatments relative to the last treatment.
#'
#' * \code{V_direct_carry}: The covariance matrix for direct and
#'   carryover treatment effects.
#'
#' * \code{v_direct_only}: The variance of the direct treatment effect for
#'   the treatment pair of interest without accounting for carryover
#'   effects.
#'
#' * \code{v_direct}: The variance of the direct treatment effect for
#'   the treatment pair of interest accounting for carryover effects.
#'
#' * \code{v_carry}: The variance of the carryover treatment effect for
#'   the treatment pair of interest.
#'
#' * \code{releff_direct}: The relative efficiency of the design
#'   for estimating the direct treatment effect for the treatment pair
#'   of interest after accounting for carryover effects with respect to
#'   that without accounting for carryover effects. This is equal to
#'   \code{v_direct_only/v_direct}.
#'
#' * \code{releff_carry}: The relative efficiency of the design
#'   for estimating the carryover effect for the treatment pair
#'   of interest. This is equal to \code{v_direct_only/v_carry}.
#'
#' * \code{half_width}: The half-width of the confidence interval for the
#'   direct treatment effect for the treatment pair of interest.
#'
#' * \code{nu}: Degrees of freedom for the t-test.
#'
#' * \code{allocationRatioPlanned}: Allocation ratio for the sequences.
#'
#' * \code{normalApproximation}: The type of computation of the p-values.
#'   If \code{TRUE}, the variance is assumed to be known, otherwise
#'   the calculations are performed with the t distribution.
#'
#' * \code{rounding}: Whether to round up the sample size.
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @references
#'
#' Robert O. Kuehl. Design of Experiments: Statistical Principles of
#' Research Design and Analysis. Brooks/Cole: Pacific Grove, CA. 2000.
#'
#' @examples
#'
#' # Williams design for 4 treatments
#'
#' (design1 = getDesignMeanDiffCarryoverEquiv(
#'   beta = 0.2, n = NA,
#'   meanDiffLower = -1.3, meanDiffUpper = 1.3,
#'   meanDiff = 0, stDev = 2.2,
#'   design = matrix(c(1, 4, 2, 3,
#'                     2, 1, 3, 4,
#'                     3, 2, 4, 1,
#'                     4, 3, 1, 2),
#'                   4, 4, byrow = TRUE),
#'   alpha = 0.025))
#'
#' @export
#'
getDesignMeanDiffCarryoverEquiv <- function(
    beta = NA_real_,
    n = NA_real_,
    trtpair = NA_real_,
    carryover = TRUE,
    meanDiffLower = NA_real_,
    meanDiffUpper = NA_real_,
    meanDiff = 0,
    stDev = 1,
    corr = 0.5,
    design = NA_real_,
    cumdrop = NA_real_,
    allocationRatioPlanned = NA_real_,
    normalApproximation = FALSE,
    rounding = TRUE,
    alpha = 0.025) {

  if (is.na(beta) && is.na(n)) {
    stop("beta and n cannot be both missing")
  }

  if (!is.na(beta) && !is.na(n)) {
    stop("Only one of beta and n should be provided")
  }

  if (!is.na(beta) && (beta >= 1-alpha || beta < 0.0001)) {
    stop("beta must lie in [0.0001, 1-alpha)")
  }

  if (!is.na(n) && n <= 0) {
    stop("n must be positive")
  }

  if (is.na(meanDiffLower)) {
    stop("meanDiffLower must be provided")
  }

  if (is.na(meanDiffUpper)) {
    stop("meanDiffUpper must be provided")
  }

  if (meanDiffLower >= meanDiff) {
    stop("meanDiffLower must be less than meanDiff")
  }

  if (meanDiffUpper <= meanDiff) {
    stop("meanDiffUpper must be greater than meanDiff")
  }

  if (stDev <= 0) {
    stop("stDev must be positive")
  }

  if (corr <= -1 || corr >= 1) {
    stop("corr must lie between -1 and 1")
  }

  if (any(is.na(design))) {
    stop("design must be provided")
  }

  nseq = nrow(design)
  nprd = ncol(design)
  ntrt = length(unique(c(design)))

  if (any(design <= 0 | design != round(design)) ||
      !all.equal(sort(unique(c(design))), 1:ntrt)) {
    stop(paste("Elements of design must be positive integers",
               "ranging from 1 to", ntrt))
  }

  if (any(is.na(trtpair))) {
    trtpair = c(1, ntrt)
  } else if (length(trtpair) != 2 ||
             trtpair[1] == trtpair[2] ||
             !all(trtpair %in% 1:ntrt)) {
    stop(paste("trtpair must be a vector of two distinct integers",
               "selected from 1 to", ntrt))
  }

  # number of model parameters consisting of the intercept, sequence effects,
  # period effects, direct treatment effects, and carryover treatment effects
  q = 1 + (nseq-1) + (nprd-1) + (ntrt-1) + (ntrt-1)

  # number of model parameters consisting of the intercept, sequence effects,
  # period effects, and direct treatment effects
  q0 = 1 + (nseq-1) + (nprd-1) + (ntrt-1)

  # offset of direct treatment effect
  l = 1 + (nseq-1) + (nprd-1)

  # start of direct treatment effect
  l1 = l + 1

  # end of direct treatment effect
  m = 1 + (nseq-1) + (nprd-1) + (ntrt-1)

  if (any(is.na(cumdrop))) {
    cumdrop = rep(0, nprd)
  } else if (length(cumdrop) != nprd) {
    stop(paste("cumdrop should have", nprd, "elements"))
  }

  if (any(cumdrop < 0)) {
    stop("Elements of cumdrop must be nonnegative")
  } else if (nprd > 1 && any(diff(cumdrop) < 0)) {
    stop("Elements of cumdrop must be nondecreasing")
  } else if (any(cumdrop >= 1)) {
    stop("Elements of cumdrop must be less than 1")
  }

  # observed data pattern probabilities
  p = diff(c(cumdrop, 1))


  if (any(is.na(allocationRatioPlanned))) {
    allocationRatioPlanned = rep(1, nrow(design))
  }

  if (length(allocationRatioPlanned) != nseq-1 &&
      length(allocationRatioPlanned) != nseq) {
    stop(paste("allocationRatioPlanned should have", nseq-1,
               "or", nseq, "elements"))
  }

  if (length(allocationRatioPlanned) == nseq-1) {
    allocationRatioPlanned = c(allocationRatioPlanned, 1)
  }

  if (any(allocationRatioPlanned <= 0)) {
    stop("Elements of allocationRatioPlanned must be positive")
  }

  if (alpha < 0.00001 || alpha >= 1) {
    stop("alpha must lie in [0.00001, 1)")
  }


  # treatment sequence randomization probabilities
  r = allocationRatioPlanned/sum(allocationRatioPlanned)


  # compound symmetry covariance matrix for repeated measures
  Sigma = stDev^2/(1-corr) * ((1-corr)*diag(nprd) + corr)

  # model design matrix with carryover effect
  X = matrix(0, nseq*nprd, q)
  X[,1] = 1
  for (i in 1:nseq) {
    for (j in 1:nprd) {
      k = (i-1)*nprd + j

      # sequence effects
      if (i < nseq) {
        X[k, 1 + i] = 1
      }

      # period effects
      if (j < nprd) {
        X[k, 1 + (nseq-1) + j] = 1
      }

      # direct treatment effects
      if (design[i,j] < ntrt) {
        X[k, 1 + (nseq-1) + (nprd-1) + design[i,j]] = 1
      }

      # carryover treatment effects
      if (j > 1 && design[i,j-1] < ntrt) {
        X[k, 1 + (nseq-1) + (nprd-1) + (ntrt-1) + design[i,j-1]] = 1
      }
    }
  }

  # singular values of X
  d = svdcpp(X)$d
  tol = max(dim(X))*.Machine$double.eps
  if (carryover && sum(d >= tol * max(d)) < q) {
    stop("The crossover design is overparameterized for carryover effects")
  }

  if (carryover) {
    # information matrix for model parameters with carryover effects
    I = 0
    for (i in 1:nseq) {
      offset = (i-1)*nprd
      J = 0
      for (j in 1:nprd) {
        idx = offset + (1:j)

        if (j==1) {
          x = t(as.matrix(X[idx,]))
        } else {
          x = X[idx,]
        }

        J = J + p[j]*t(x) %*% solve(Sigma[1:j,1:j]) %*% x
      }
      I = I + r[i]*J
    }

    # covariance matrix for model parameters with carryover effects
    V = solve(I)

    # variance for direct treatment effect for one sampling unit with carryover
    if (trtpair[2]==ntrt) {
      v1 = V[l+trtpair[1], l+trtpair[1]]
    } else if (trtpair[1]==ntrt) {
      v1 = V[l+trtpair[2], l+trtpair[2]]
    } else {
      v1 = V[l+trtpair[1], l+trtpair[1]] +
        V[l+trtpair[2], l+trtpair[2]] -
        2*V[l+trtpair[1], l+trtpair[2]]
    }


    # variance for carryover treatment effect for one sampling unit
    if (trtpair[2]==ntrt) {
      v2 = V[m+trtpair[1], m+trtpair[1]]
    } else if (trtpair[1]==ntrt) {
      v2 = V[m+trtpair[2], m+trtpair[2]]
    } else {
      v2 = V[m+trtpair[1], m+trtpair[1]] +
        V[m+trtpair[2], m+trtpair[2]] -
        2*V[m+trtpair[1], m+trtpair[2]]
    }
  }


  # design matrix for model parameters without carryover effects
  X0 = X[, 1:m]
  if (!carryover) {
    # singular values of X0
    d0 = svdcpp(X0)$d
    tol0 = max(dim(X0))*.Machine$double.eps
    if (sum(d0 >= tol0 * max(d0)) < q0) {
      stop("The crossover design is overparameterized")
    }
  }

  # information matrix for model parameters without carryover effects
  I0 = 0
  for (i in 1:nseq) {
    offset = (i-1)*nprd
    J = 0
    for (j in 1:nprd) {
      idx = offset+(1:j)

      if (j==1) {
        x = t(as.matrix(X0[idx,]))
      } else {
        x = X0[idx,]
      }

      J = J + p[j]*t(x) %*% solve(Sigma[1:j,1:j]) %*% x
    }
    I0 = I0 + r[i]*J
  }

  # covariance matrix for model parameters without carryover effects
  V0 = solve(I0)

  # variance for direct treatment effect for one sampling unit w/o carryover
  if (trtpair[2]==ntrt) {
    v0 = V0[l+trtpair[1], l+trtpair[1]]
  } else if (trtpair[1]==ntrt) {
    v0 = V0[l+trtpair[2], l+trtpair[2]]
  } else {
    v0 = V0[l+trtpair[1], l+trtpair[1]] +
      V0[l+trtpair[2], l+trtpair[2]] -
      2*V0[l+trtpair[1], l+trtpair[2]]
  }

  if (carryover) {
    v = v1
  } else {
    v = v0
  }

  # power for two one-sided t-tests
  f = function(n) {
    # residual degrees of freedom after accounting for the subject effects
    mean_nprd = sum(p*(1:nprd))
    nu = n*mean_nprd - 1 - (n-1) - (nprd-1) - (ntrt-1) - (ntrt-1)*carryover
    b = qt(1-alpha, nu)
    ncpLower = (meanDiff - meanDiffLower)*sqrt(n/v)
    powerLower = pt(b, nu, ncpLower, lower.tail = FALSE)
    ncpUpper = (meanDiffUpper - meanDiff)*sqrt(n/v)
    powerUpper = pt(b, nu, ncpUpper, lower.tail = FALSE)
    power = powerLower + powerUpper - 1
    power
  }

  if (is.na(n)) { # calculate the sample size
    des = getDesignEquiv(
      beta = beta, IMax = NA, thetaLower = meanDiffLower,
      thetaUpper = meanDiffUpper, theta = meanDiff)

    n = des$overallResults$information*v

    if (!normalApproximation) {
      n = uniroot(function(n) f(n) - (1-beta), c(0.5*n, 1.5*n))$root
    }
  }

  if (rounding) {
    n = ceiling(n - 1.0e-12)
  }

  if (normalApproximation) {
    des = getDesignEquiv(
      beta = NA, IMax = n/v, thetaLower = meanDiffLower,
      thetaUpper = meanDiffUpper, theta = meanDiff)

    power = des$overallResults$overallReject

  } else {
    power = f(n)
  }

  rownames(design) = paste0("Seq", 1:nseq)
  colnames(design) = paste0("Prd", 1:nprd)

  if (!normalApproximation) {
    mean_nprd = sum(p*(1:nprd))
    nu = n*mean_nprd - 1 - (n-1) - (nprd-1) - (ntrt-1) - (ntrt-1)*carryover
    hw = qt(1-alpha, nu)*sqrt(v/n)
  } else {
    hw = qnorm(1-alpha)*sqrt(v/n)
  }

  if (carryover) {
    des = list(
      power = power, alpha = alpha, numberOfSubjects = n,
      trtpair = trtpair, carryover = carryover,
      meanDiffLower = meanDiffLower,
      meanDiffUpper = meanDiffUpper,
      meanDiff = meanDiff,
      stDev = stDev, corr = corr,
      design = design, designMatrix = X,
      nseq = nseq, nprd = nprd, ntrt = ntrt,
      cumdrop = cumdrop,
      V_direct_only = V0[l1:m,l1:m]/n,
      V_direct_carry = V[l1:q, l1:q]/n,
      v_direct_only = v0/n,
      v_direct = v1/n,
      v_carry = v2/n,
      releff_direct = v0/v1,
      releff_carry = v0/v2,
      half_width = hw, nu = nu,
      allocationRatioPlanned = allocationRatioPlanned,
      normalApproximation = normalApproximation,
      rounding = rounding)
  } else {
    des = list(
      power = power, alpha = alpha, numberOfSubjects = n,
      trtpair = trtpair, carryover = carryover,
      meanDiffLower = meanDiffLower,
      meanDiffUpper = meanDiffUpper,
      meanDiff = meanDiff,
      stDev = stDev, corr = corr,
      design = design, designMatrix = X0,
      nseq = nseq, nprd = nprd, ntrt = ntrt,
      cumdrop = cumdrop,
      V_direct_only = V0[l1:m,l1:m]/n,
      v_direct_only = v0/n,
      half_width = hw, nu = nu,
      allocationRatioPlanned = allocationRatioPlanned,
      normalApproximation = normalApproximation,
      rounding = rounding)
  }

  attr(des, "class") = "designMeanDiffCarryoverEquiv"

  des
}


#' @title Power and Sample Size for One-Way ANOVA
#' @description Obtains the power and sample size for one-way
#' analysis of variance.
#'
#' @param beta The type II error.
#' @param n The total sample size.
#' @param ngroups The number of treatment groups.
#' @param means The treatment group means.
#' @param stDev The common standard deviation.
#' @param allocationRatioPlanned Allocation ratio for the treatment
#'   groups. It has length \code{ngroups - 1} or \code{ngroups}. If it is
#'   of length \code{ngroups - 1}, then the last treatment group will
#'   assume value 1 for allocation ratio.
#' @param rounding Whether to round up sample size. Defaults to 1 for
#'   sample size rounding.
#' @param alpha The two-sided significance level. Defaults to 0.05.
#'
#' @details
#'
#' Let \eqn{\{\mu_i: i=1,\ldots,k\}} denote the group means, and
#' \eqn{\{r_i: i=1,\ldots,k\}} denote the randomization probabilities
#' to the \eqn{k} treatment groups. Let \eqn{\sigma} denote the
#' common standard deviation, and \eqn{n} denote the total sample
#' size. Then the \eqn{F}-statistic
#' \deqn{F = \frac{SSR/(k-1)}{SSE/(n-k)}
#' \sim F_{k-1, n-k, \lambda}} where
#' \deqn{\lambda = n \sum_{i=1}^k r_i (\mu_i - \bar{\mu})^2/\sigma^2}
#' is the noncentrality parameter, and
#' \eqn{\bar{\mu} = \sum_{i=1}^k r_i \mu_i}.
#'
#' @return An S3 class \code{designANOVA} object with the following
#' components:
#'
#' * \code{power}: The power to reject the null hypothesis that
#'   there is no difference among the treatment groups.
#'
#' * \code{alpha}: The two-sided significance level.
#'
#' * \code{n}: The number of subjects.
#'
#' * \code{ngroups}: The number of treatment groups.
#'
#' * \code{means}: The treatment group means.
#'
#' * \code{stDev}: The common standard deviation.
#'
#' * \code{effectsize}: The effect size.
#'
#' * \code{allocationRatioPlanned}: Allocation ratio for the treatment
#'   groups.
#'
#' * \code{rounding}: Whether to round up sample size.
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @examples
#'
#' (design1 <- getDesignANOVA(
#'   beta = 0.1, ngroups = 4, means = c(1.5, 2.5, 2, 0),
#'   stDev = 3.5, allocationRatioPlanned = c(2, 2, 2, 1),
#'   alpha = 0.05))
#'
#' @export
#'
getDesignANOVA <- function(
    beta = NA_real_,
    n = NA_real_,
    ngroups = 2,
    means = NA_real_,
    stDev = 1,
    allocationRatioPlanned = NA_real_,
    rounding = TRUE,
    alpha = 0.05) {

  if (is.na(beta) && is.na(n)) {
    stop("beta and n cannot be both missing")
  }

  if (!is.na(beta) && !is.na(n)) {
    stop("Only one of beta and n should be provided")
  }

  if (!is.na(beta) && (beta >= 1-alpha || beta < 0.0001)) {
    stop("beta must lie in [0.0001, 1-alpha)")
  }

  if (!is.na(n) && n <= 0) {
    stop("n must be positive")
  }

  if (!length(means) == ngroups) {
    stop(paste("means must have", ngroups, "elements"))
  }

  if (stDev <= 0) {
    stop("stDev must be positive")
  }

  if (any(is.na(allocationRatioPlanned))) {
    allocationRatioPlanned = rep(1, ngroups)
  }

  if (length(allocationRatioPlanned) != ngroups - 1 &&
      length(allocationRatioPlanned) != ngroups) {
    stop(paste("allocationRatioPlanned should have", ngroups - 1,
               "or", ngroups, "elements"))
  }

  if (length(allocationRatioPlanned) == ngroups - 1) {
    allocationRatioPlanned = c(allocationRatioPlanned, 1)
  }

  if (any(allocationRatioPlanned <= 0)) {
    stop("Elements of allocationRatioPlanned must be positive")
  }

  if (alpha < 0.00001 || alpha >= 1) {
    stop("alpha must lie in [0.00001, 1)")
  }


  r = allocationRatioPlanned/sum(allocationRatioPlanned)

  mubar = sum(r*means)
  vmu = sum(r*(means - mubar)^2)

  # power for F-test
  f <- function(n) {
    lambda = n*vmu/stDev^2
    b = qf(1 - alpha, ngroups - 1, n - ngroups)
    pf(b, ngroups - 1, n - ngroups, lambda, lower.tail = FALSE)
  }

  if (is.na(n)) {
    nu = ngroups - 1
    n0 = (qchisq(1-alpha, nu) - nu + qnorm(1-beta)*sqrt(2*nu))/(vmu/stDev^2)
    while (f(n0) < 1-beta) n0 <- 2*n0
    n = uniroot(function(n) f(n) - (1-beta), c(0.5*n0, n0))$root
  }

  if (rounding) {
    n = ceiling(n - 1.0e-12)
  }

  power = f(n)

  des = list(
    power = power, alpha = alpha, n = n,
    ngroups = ngroups, means = means, stDev = stDev,
    effectsize = vmu/stDev^2,
    allocationRatioPlanned = allocationRatioPlanned,
    rounding = rounding)

  attr(des, "class") = "designANOVA"

  des
}


#' @title Power and Sample Size for Two-Way ANOVA
#' @description Obtains the power and sample size for two-way
#' analysis of variance.
#'
#' @param beta The type II error.
#' @param n The total sample size.
#' @param nlevelsA The number of groups for Factor A.
#' @param nlevelsB The number of levels for Factor B.
#' @param means The matrix of treatment means for Factors A and B
#'   combination.
#' @param stDev The common standard deviation.
#' @param rounding Whether to round up sample size. Defaults to 1 for
#'   sample size rounding.
#' @param alpha The two-sided significance level. Defaults to 0.05.
#'
#' @return An S3 class \code{designTwoWayANOVA} object with the following
#' components:
#'
#' * \code{alpha}: The two-sided significance level.
#'
#' * \code{nlevelsA}: The number of levels for Factor A.
#'
#' * \code{nlevelsB}: The number of levels for Factor B.
#'
#' * \code{means}: The matrix of treatment group means.
#'
#' * \code{stDev}: The common standard deviation.
#'
#' * \code{effectsizeA}: The effect size for Factor A.
#'
#' * \code{effectsizeB}: The effect size for Factor B.
#'
#' * \code{effectsizeAB}: The effect size for Factor A and Factor B
#'   interaction.
#'
#' * \code{rounding}: Whether to round up sample size.
#'
#' * \code{powerdf}: The data frame containing the power and sample size
#'   results. It has the following variables:
#'
#'     - \code{n}: The sample size.
#'
#'     - \code{powerA}: The power to reject the null hypothesis that
#'       there is no difference among Factor A levels.
#'
#'     - \code{powerB}: The power to reject the null hypothesis that
#'       there is no difference among Factor B levels.
#'
#'     - \code{powerAB}: The power to reject the null hypothesis that
#'       there is no interaction between Factor A and Factor B.
#'
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @examples
#'
#' (design1 <- getDesignTwoWayANOVA(
#'   beta = 0.1, nlevelsA = 2, nlevelsB = 2,
#'   means = matrix(c(0.5, 4.7, 0.4, 6.9), 2, 2, byrow = TRUE),
#'   stDev = 2, alpha = 0.05))
#'
#' @export
#'
getDesignTwoWayANOVA <- function(
    beta = NA_real_,
    n = NA_real_,
    nlevelsA = 2,
    nlevelsB = 2,
    means = NA_real_,
    stDev = 1,
    rounding = TRUE,
    alpha = 0.05) {

  if (is.na(beta) && is.na(n)) {
    stop("beta and n cannot be both missing")
  }

  if (!is.na(beta) && !is.na(n)) {
    stop("Only one of beta and n should be provided")
  }

  if (!is.na(beta) && (beta >= 1-alpha || beta < 0.0001)) {
    stop("beta must lie in [0.0001, 1-alpha)")
  }

  if (!is.na(n) && n <= 0) {
    stop("n must be positive")
  }

  if (nlevelsA < 2 || nlevelsA != round(nlevelsA)) {
    stop("nlevelsA must be a positive integer >= 2")
  }

  if (nlevelsB < 2 || nlevelsB != round(nlevelsB)) {
    stop("nlevelsB must be a positive integer >= 2")
  }

  if (any(is.na(means))) {
    stop("means must be provided")
  }

  if (!is.matrix(means) || nrow(means) != nlevelsA ||
      ncol(means) !=  nlevelsB) {
    stop(paste("means must be a matrix with", nlevelsA, "rows and",
               nlevelsB, "columns"))
  }

  if (stDev <= 0) {
    stop("stDev must be positive")
  }

  if (alpha < 0.00001 || alpha >= 1) {
    stop("alpha must lie in [0.00001, 1)")
  }


  muA = as.numeric(rowMeans(means))
  muB = as.numeric(colMeans(means))
  mu = mean(means)

  vmuA = var(muA)*(nlevelsA - 1)/nlevelsA
  vmuB = var(muB)*(nlevelsB - 1)/nlevelsB

  mmuA = matrix(muA, nlevelsA, nlevelsB)
  mmuB = matrix(muB, nlevelsA, nlevelsB, byrow = TRUE)
  vmuAB = sum((means - mmuA - mmuB + mu)^2)/(nlevelsA*nlevelsB)

  effectsizeA = vmuA/stDev^2
  effectsizeB = vmuB/stDev^2
  effectsizeAB = vmuAB/stDev^2


  # power for F-test for Factor A
  fA <- function(n) {
    lambda = n*effectsizeA
    b = qf(1-alpha, nlevelsA - 1, n - nlevelsA*nlevelsB)
    pf(b, nlevelsA - 1, n - nlevelsA*nlevelsB, lambda, lower.tail = FALSE)
  }

  # power for F-test for Factor A
  fB <- function(n) {
    lambda = n*effectsizeB
    b = qf(1-alpha, nlevelsB - 1, n - nlevelsA*nlevelsB)
    pf(b, nlevelsB - 1, n - nlevelsA*nlevelsB, lambda, lower.tail = FALSE)
  }

  # power for F-test for Factor A and Factor B interaction
  fAB <- function(n) {
    lambda = n*effectsizeAB
    b = qf(1-alpha, (nlevelsA - 1)*(nlevelsB - 1), n - nlevelsA*nlevelsB)
    pf(b, (nlevelsA - 1)*(nlevelsB - 1), n - nlevelsA*nlevelsB,
       lambda, lower.tail = FALSE)
  }


  if (is.na(n)) {
    nu = nlevelsA - 1
    n0 = (qchisq(1-alpha, nu) - nu + qnorm(1-beta)*sqrt(2*nu))/effectsizeA
    n0 = max(n0, nlevelsA*nlevelsB + 1)
    while (fA(n0) < 1-beta) n0 <- 2*n0
    nA = uniroot(function(n) fA(n) - (1-beta), c(0.5*n0, n0))$root

    nu = nlevelsB - 1
    n0 = (qchisq(1-alpha, nu) - nu + qnorm(1-beta)*sqrt(2*nu))/effectsizeB
    n0 = max(n0, nlevelsA*nlevelsB + 1)
    while (fB(n0) < 1-beta) n0 <- 2*n0
    nB = uniroot(function(n) fB(n) - (1-beta), c(0.5*n0, n0))$root

    nu = (nlevelsA - 1)*(nlevelsB - 1)
    n0 = (qchisq(1-alpha, nu) - nu + qnorm(1-beta)*sqrt(2*nu))/effectsizeAB
    n0 = max(n0, nlevelsA*nlevelsB + 1)
    while (fAB(n0) < 1-beta) n0 <- 2*n0
    nAB = uniroot(function(n) fAB(n) - (1-beta), c(0.5*n0, n0))$root

    # ensure that nA, nB, and nAB are multiples of 4
    nA = ceiling(nA/4) * 4
    nB = ceiling(nB/4) * 4
    nAB = ceiling(nAB/4) * 4
  } else {
    nA = n
    nB = n
    nAB = n
  }

  if (rounding) {
    nA = ceiling(nA - 1.0e-12)
    nB = ceiling(nB - 1.0e-12)
    nAB = ceiling(nAB - 1.0e-12)
  }

  m = unique(c(nA, nB, nAB))
  powerdf = data.frame(n = m, powerA = sapply(m, fA),
                       powerB = sapply(m, fB), powerAB = sapply(m, fAB))
  des = list(
    alpha = alpha, nlevelsA = nlevelsA,
    nlevelsB = nlevelsB, means = means, stDev = stDev,
    effectsizeA = effectsizeA, effectsizeB = effectsizeB,
    effectsizeAB = effectsizeAB,
    rounding = rounding,
    powerdf = powerdf)

  attr(des, "class") = "designTwoWayANOVA"

  des
}


#' @title Power and Sample Size for One-Way ANOVA Contrast
#' @description Obtains the power and sample size for a single contrast
#' in one-way analysis of variance.
#'
#' @param beta The type II error.
#' @param n The total sample size.
#' @param ngroups The number of treatment groups.
#' @param means The treatment group means.
#' @param stDev The common standard deviation.
#' @param contrast The coefficients for the single contrast.
#' @param meanContrastH0 The mean of the contrast under the
#'   null hypothesis.
#' @param allocationRatioPlanned Allocation ratio for the treatment
#'   groups. It has length \code{ngroups - 1} or \code{ngroups}. If it is
#'   of length \code{ngroups - 1}, then the last treatment group will
#'   assume value 1 for allocation ratio.
#' @param rounding Whether to round up sample size. Defaults to 1 for
#'   sample size rounding.
#' @param alpha The one-sided significance level. Defaults to 0.025.
#'
#' @return An S3 class \code{designANOVAContrast} object with the following
#' components:
#'
#' * \code{power}: The power to reject the null hypothesis for the
#'   treatment contrast.
#'
#' * \code{alpha}: The one-sided significance level.
#'
#' * \code{n}: The number of subjects.
#'
#' * \code{ngroups}: The number of treatment groups.
#'
#' * \code{means}: The treatment group means.
#'
#' * \code{stDev}: The common standard deviation.
#'
#' * \code{contrast}: The coefficients for the single contrast.
#'
#' * \code{meanContrastH0}: The mean of the contrast under the null
#'   hypothesis.
#'
#' * \code{meanContrast}: The mean of the contrast under the alternative
#'   hypothesis.
#'
#' * \code{effectsize}: The effect size.
#'
#' * \code{allocationRatioPlanned}: Allocation ratio for the treatment
#'   groups.
#'
#' * \code{rounding}: Whether to round up sample size.
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @examples
#'
#' (design1 <- getDesignANOVAContrast(
#'   beta = 0.1, ngroups = 4, means = c(1.5, 2.5, 2, 0),
#'   stDev = 3.5, contrast = c(1, 1, 1, -3),
#'   allocationRatioPlanned = c(2, 2, 2, 1),
#'   alpha = 0.025))
#'
#' @export
#'
getDesignANOVAContrast <- function(
    beta = NA_real_,
    n = NA_real_,
    ngroups = 2,
    means = NA_real_,
    stDev = 1,
    contrast = NA_real_,
    meanContrastH0 = 0,
    allocationRatioPlanned = NA_real_,
    rounding = TRUE,
    alpha = 0.025) {

  if (is.na(beta) && is.na(n)) {
    stop("beta and n cannot be both missing")
  }

  if (!is.na(beta) && !is.na(n)) {
    stop("Only one of beta and n should be provided")
  }

  if (!is.na(beta) && (beta >= 1-alpha || beta < 0.0001)) {
    stop("beta must lie in [0.0001, 1-alpha)")
  }

  if (!is.na(n) && n <= 0) {
    stop("n must be positive")
  }

  if (length(means) != ngroups) {
    stop(paste("means must have", ngroups, "elements"))
  }

  if (stDev <= 0) {
    stop("stDev must be positive")
  }

  if (any(is.na(contrast))) {
    stop("contrast must be provided")
  }

  if (length(contrast) != ngroups) {
    stop(paste("contrast must have", ngroups, "elements"))
  }

  if (any(is.na(allocationRatioPlanned))) {
    allocationRatioPlanned = rep(1, ngroups)
  }

  if (length(allocationRatioPlanned) != ngroups - 1 &&
      length(allocationRatioPlanned) != ngroups) {
    stop(paste("allocationRatioPlanned should have", ngroups - 1,
               "or", ngroups, "elements"))
  }

  if (length(allocationRatioPlanned) == ngroups - 1) {
    allocationRatioPlanned = c(allocationRatioPlanned, 1)
  }

  if (any(allocationRatioPlanned <= 0)) {
    stop("Elements of allocationRatioPlanned must be positive")
  }

  if (alpha < 0.00001 || alpha >= 1) {
    stop("alpha must lie in [0.00001, 1)")
  }


  r = allocationRatioPlanned/sum(allocationRatioPlanned)

  meanContrast = sum(contrast*means)
  v1 = sum(contrast^2/r)*stDev^2

  directionUpper = meanContrast > meanContrastH0

  theta = ifelse(directionUpper, meanContrast - meanContrastH0,
                 meanContrastH0 - meanContrast)

  # power for t-test
  f <- function(n) {
    b = qt(1-alpha, n-ngroups)
    ncp = theta*sqrt(n/v1)
    power = pt(b, n-ngroups, ncp, lower.tail = FALSE)
  }

  if (is.na(n)) {
    n0 = (qnorm(1-alpha) + qnorm(1-beta))^2*v1/theta^2
    n = uniroot(function(n) f(n) - (1-beta), c(n0, 2*n0))$root
  }

  if (rounding) {
    n = ceiling(n - 1.0e-12)
  }

  power = f(n)

  des = list(
    power = power, alpha = alpha, n = n,
    ngroups = ngroups, means = means, stDev = stDev,
    contrast = contrast, meanContrastH0 = meanContrastH0,
    meanContrast = meanContrast, effectsize = theta^2/v1,
    allocationRatioPlanned = allocationRatioPlanned,
    rounding = rounding)

  attr(des, "class") = "designANOVAContrast"

  des
}


#' @title Power and Sample Size for Repeated-Measures ANOVA
#' @description Obtains the power and sample size for one-way repeated
#' measures analysis of variance. Each subject takes all treatments
#' in the longitudinal study.
#'
#' @param beta The type II error.
#' @param n The total sample size.
#' @param ngroups The number of treatment groups.
#' @param means The treatment group means.
#' @param stDev The total standard deviation.
#' @param corr The correlation among the repeated measures.
#' @param rounding Whether to round up sample size. Defaults to 1 for
#'   sample size rounding.
#' @param alpha The two-sided significance level. Defaults to 0.05.
#'
#' @details
#' Let \eqn{y_{ij}} denote the measurement under treatment condition
#' \eqn{j (j=1,\ldots,k)} for subject \eqn{i (i=1,\ldots,n)}. Then
#' \deqn{y_{ij} = \alpha + \beta_j + b_i + e_{ij}} where \eqn{b_i}
#' denotes the subject random effect, \eqn{b_i \sim N(0, \sigma_b^2)}
#' and \eqn{e_{ij} \sim N(0, \sigma_e^2)} denotes the within-subject
#' residual. If we set \eqn{\beta_k = 0}, then \eqn{\alpha} is the
#' mean of the last treatment (control), and \eqn{\beta_j} is the
#' difference in means between the \eqn{j}th treatment and the control
#' for \eqn{j=1,\ldots,k-1}.
#'
#' The repeated measures have a compound symmetry covariance structure.
#' Let \eqn{\sigma^2 = \sigma_b^2 + \sigma_e^2}, and
#' \eqn{\rho = \frac{\sigma_b^2}{\sigma_b^2 + \sigma_e^2}}. Then
#' \eqn{Var(y_i) = \sigma^2 \{(1-\rho) I_k + \rho 1_k 1_k^T\}}.
#' Let \eqn{X_i} denote the design matrix for subject \eqn{i}.
#' Let \eqn{\theta = (\alpha, \beta_1, \ldots, \beta_{k-1})^T}.
#' It follows that
#' \deqn{Var(\hat{\theta}) = \left(\sum_{i=1}^{n} X_i^T V_i^{-1}
#' X_i\right)^{-1}.} It can be shown that
#' \deqn{Var(\hat{\beta}) = \frac{\sigma^2 (1-\rho)}{n} (I_{k-1} +
#' 1_{k-1} 1_{k-1}^T).} It follows that
#' \eqn{\hat{\beta}^T \hat{V}_{\hat{\beta}}^{-1} \hat{\beta} \sim
#' F_{k-1,(n-1)(k-1), \lambda}} where the noncentrality parameter for
#' the \eqn{F} distribution is \deqn{\lambda =
#' \beta^T V_{\hat{\beta}}^{-1} \beta = \frac{n \sum_{j=1}^{k}
#' (\mu_j - \bar{\mu})^2}{\sigma^2(1-\rho)}.}
#'
#' @return An S3 class \code{designRepeatedANOVA} object with the
#' following components:
#'
#' * \code{power}: The power to reject the null hypothesis that
#'   there is no difference among the treatment groups.
#'
#' * \code{alpha}: The two-sided significance level.
#'
#' * \code{n}: The number of subjects.
#'
#' * \code{ngroups}: The number of treatment groups.
#'
#' * \code{means}: The treatment group means.
#'
#' * \code{stDev}: The total standard deviation.
#'
#' * \code{corr}: The correlation among the repeated measures.
#'
#' * \code{effectsize}: The effect size.
#'
#' * \code{rounding}: Whether to round up sample size.
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @examples
#'
#' (design1 <- getDesignRepeatedANOVA(
#'   beta = 0.1, ngroups = 4, means = c(1.5, 2.5, 2, 0),
#'   stDev = 5, corr = 0.2, alpha = 0.05))
#'
#' @export
#'
getDesignRepeatedANOVA <- function(
    beta = NA_real_,
    n = NA_real_,
    ngroups = 2,
    means = NA_real_,
    stDev = 1,
    corr = 0,
    rounding = TRUE,
    alpha = 0.05) {

  if (is.na(beta) && is.na(n)) {
    stop("beta and n cannot be both missing")
  }

  if (!is.na(beta) && !is.na(n)) {
    stop("Only one of beta and n should be provided")
  }

  if (!is.na(beta) && (beta >= 1-alpha || beta < 0.0001)) {
    stop("beta must lie in [0.0001, 1-alpha)")
  }

  if (!is.na(n) && n <= 0) {
    stop("n must be positive")
  }

  if (!length(means) == ngroups) {
    stop(paste("means must have", ngroups, "elements"))
  }

  if (stDev <= 0) {
    stop("stDev must be positive")
  }

  if (corr <= -1 || corr >= 1) {
    stop("corr must lie between 0 and 1")
  }

  if (alpha < 0.00001 || alpha >= 1) {
    stop("alpha must lie in [0.00001, 1)")
  }


  vmu = var(means)*(ngroups-1)/ngroups

  # power for F-test
  f <- function(n) {
    lambda = n*ngroups*vmu/(stDev^2*(1-corr))
    b = qf(1 - alpha, ngroups - 1, (n-1)*(ngroups-1))
    pf(b, ngroups - 1, (n-1)*(ngroups-1), lambda, lower.tail = FALSE)
  }

  if (is.na(n)) {
    nu = ngroups - 1
    n0 = (qchisq(1-alpha, nu) - nu + qnorm(1-beta)*sqrt(2*nu))/
      (ngroups*vmu/(stDev^2*(1-corr)))
    while (f(n0) < 1-beta) n0 <- 2*n0
    n = uniroot(function(n) f(n) - (1-beta), c(0.5*n0, n0))$root
  }

  if (rounding) {
    n = ceiling(n - 1.0e-12)
  }

  power = f(n)

  des = list(
    power = power, alpha = alpha, n = n,
    ngroups = ngroups, means = means, stDev = stDev,
    corr = corr, effectsize = ngroups*vmu/(stDev^2*(1-corr)),
    rounding = rounding)

  attr(des, "class") = "designRepeatedANOVA"

  des
}


#' @title Power and Sample Size for One-Way Repeated Measures ANOVA Contrast
#' @description Obtains the power and sample size for a single contrast
#' in one-way repeated measures analysis of variance.
#'
#' @param beta The type II error.
#' @param n The total sample size.
#' @param ngroups The number of treatment groups.
#' @param means The treatment group means.
#' @param stDev The total standard deviation.
#' @param corr The correlation among the repeated measures.
#' @param contrast The coefficients for the single contrast.
#' @param meanContrastH0 The mean of the contrast under the
#'   null hypothesis.
#' @param rounding Whether to round up sample size. Defaults to 1 for
#'   sample size rounding.
#' @param alpha The one-sided significance level. Defaults to 0.025.
#'
#' @return An S3 class \code{designRepeatedANOVAContrast} object with
#' the following components:
#'
#' * \code{power}: The power to reject the null hypothesis for the
#'   treatment contrast.
#'
#' * \code{alpha}: The one-sided significance level.
#'
#' * \code{n}: The number of subjects.
#'
#' * \code{ngroups}: The number of treatment groups.
#'
#' * \code{means}: The treatment group means.
#'
#' * \code{stDev}: The total standard deviation.
#'
#' * \code{corr}: The correlation among the repeated measures.
#'
#' * \code{contrast}: The coefficients for the single contrast.
#'
#' * \code{meanContrastH0}: The mean of the contrast under the null
#'   hypothesis.
#'
#' * \code{meanContrast}: The mean of the contrast under the alternative
#'   hypothesis.
#'
#' * \code{effectsize}: The effect size.
#'
#' * \code{rounding}: Whether to round up sample size.
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @examples
#'
#' (design1 <- getDesignRepeatedANOVAContrast(
#'   beta = 0.1, ngroups = 4, means = c(1.5, 2.5, 2, 0),
#'   stDev = 5, corr = 0.2, contrast = c(1, 1, 1, -3)/3,
#'   alpha = 0.025))
#'
#' @export
#'
getDesignRepeatedANOVAContrast <- function(
    beta = NA_real_,
    n = NA_real_,
    ngroups = 2,
    means = NA_real_,
    stDev = 1,
    corr = 0,
    contrast = NA_real_,
    meanContrastH0 = 0,
    rounding = TRUE,
    alpha = 0.025) {

  if (is.na(beta) && is.na(n)) {
    stop("beta and n cannot be both missing")
  }

  if (!is.na(beta) && !is.na(n)) {
    stop("Only one of beta and n should be provided")
  }

  if (!is.na(beta) && (beta >= 1-alpha || beta < 0.0001)) {
    stop("beta must lie in [0.0001, 1-alpha)")
  }

  if (!is.na(n) && n <= 0) {
    stop("n must be positive")
  }

  if (length(means) != ngroups) {
    stop(paste("means must have", ngroups, "elements"))
  }

  if (stDev <= 0) {
    stop("stDev must be positive")
  }

  if (corr <= -1 || corr >= 1) {
    stop("corr must lie between -1 and 1")
  }

  if (any(is.na(contrast))) {
    stop("contrast must be provided")
  }

  if (length(contrast) != ngroups) {
    stop(paste("contrast must have", ngroups, "elements"))
  }

  if (alpha < 0.00001 || alpha >= 1) {
    stop("alpha must lie in [0.00001, 1)")
  }


  meanContrast = sum(contrast*means)
  v1 = sum(contrast^2)*stDev^2*(1-corr)

  directionUpper = meanContrast > meanContrastH0

  theta = ifelse(directionUpper, meanContrast - meanContrastH0,
                 meanContrastH0 - meanContrast)

  # power for t-test
  f <- function(n) {
    b = qt(1-alpha, (n-1)*(ngroups-1))
    ncp = theta*sqrt(n/v1)
    power = pt(b, (n-1)*(ngroups-1), ncp, lower.tail = FALSE)
  }

  if (is.na(n)) {
    n0 = (qnorm(1-alpha) + qnorm(1-beta))^2*v1/theta^2
    n = uniroot(function(n) f(n) - (1-beta), c(n0, 2*n0))$root
  }

  if (rounding) {
    n = ceiling(n - 1.0e-12)
  }

  power = f(n)

  des = list(
    power = power, alpha = alpha, n = n,
    ngroups = ngroups, means = means, stDev = stDev, corr = corr,
    contrast = contrast, meanContrastH0 = meanContrastH0,
    meanContrast = meanContrast, effectsize = theta^2/v1,
    rounding = rounding)

  attr(des, "class") = "designRepeatedANOVAContrast"

  des
}


#' @title Group Sequential Design for One-Sample Slope
#' @description Obtains the power given sample size or obtains the sample
#' size given power for a group sequential design for one-sample slope.
#'
#' @param beta The type II error.
#' @param n The total sample size.
#' @param slopeH0 The slope under the null hypothesis.
#'   Defaults to 0.
#' @param slope The slope under the alternative hypothesis.
#' @param stDev The standard deviation of the residual.
#' @param stDevCovariate The standard deviation of the covariate.
#' @param normalApproximation The type of computation of the p-values.
#'   If \code{TRUE}, the variance is assumed to be known, otherwise
#'   the calculations are performed with the t distribution. The exact
#'   calculation using the t distribution is only implemented for the
#'   fixed design.
#' @param rounding Whether to round up sample size. Defaults to 1 for
#'   sample size rounding.
#' @inheritParams param_kMax
#' @param informationRates The information rates. Fixed prior to the trial.
#'   Defaults to \code{(1:kMax) / kMax} if left unspecified.
#' @inheritParams param_efficacyStopping
#' @inheritParams param_futilityStopping
#' @inheritParams param_criticalValues
#' @inheritParams param_alpha
#' @inheritParams param_typeAlphaSpending
#' @inheritParams param_parameterAlphaSpending
#' @inheritParams param_userAlphaSpending
#' @inheritParams param_futilityBounds
#' @inheritParams param_typeBetaSpending
#' @inheritParams param_parameterBetaSpending
#' @inheritParams param_userBetaSpending
#' @param spendingTime A vector of length \code{kMax} for the error spending
#'   time at each analysis. Defaults to missing, in which case, it is the
#'   same as \code{informationRates}.
#'
#' @details
#' We assume a simple linear regression of the form
#' \deqn{y_i = \alpha + \beta x_i + \epsilon_i}
#' where \eqn{\epsilon_i} is the residual error, which is assumed to
#' be normally distributed with mean 0 and standard deviation
#' \eqn{\sigma_\epsilon}. The covariate \eqn{x_i} is assumed to
#' be normally distributed with mean 0 and standard deviation
#' \eqn{\sigma_x}. The slope under the null hypothesis is
#' \eqn{\beta_0}, and the slope under the alternative
#' hypothesis is \eqn{\beta}. Since
#' \deqn{\hat{\beta} = \frac{\sum_{i=1}^{n} (x_i-\bar{x}) y_i}
#' {\sum_{i=1}^{n}(x_i-\bar{x})^2}}
#' it follows that
#' \deqn{\hat{\beta} \sim N(\beta,
#'  \frac{\sigma_\epsilon^2}{\sum_{i=1}^{n}(x_i-\bar{x})^2}).}
#' Since the variance of \eqn{\hat{\beta}} is
#' \deqn{\frac{\sigma_\epsilon^2}{n\sigma_x^2}} we can use it to
#' calculate the power and sample size for the group sequential design.
#'
#' @return An S3 class \code{designOneSlope} object with three components:
#'
#' * \code{overallResults}: A data frame containing the following variables:
#'
#'     - \code{overallReject}: The overall rejection probability.
#'
#'     - \code{alpha}: The overall significance level.
#'
#'     - \code{attainedAlpha}: The attained significance level, which is
#'       different from the overall significance level in the presence of
#'       futility stopping.
#'
#'     - \code{kMax}: The number of stages.
#'
#'     - \code{theta}: The parameter value.
#'
#'     - \code{information}: The maximum information.
#'
#'     - \code{expectedInformationH1}: The expected information under H1.
#'
#'     - \code{expectedInformationH0}: The expected information under H0.
#'
#'     - \code{drift}: The drift parameter, equal to
#'       \code{theta*sqrt(information)}.
#'
#'     - \code{inflationFactor}: The inflation factor (relative to the
#'       fixed design).
#'
#'     - \code{numberOfSubjects}: The maximum number of subjects.
#'
#'     - \code{expectedNumberOfSubjectsH1}: The expected number of subjects
#'       under H1.
#'
#'     - \code{expectedNumberOfSubjectsH0}: The expected number of subjects
#'       under H0.
#'
#'     - \code{slopeH0}: The slope under the null hypothesis.
#'
#'     - \code{slope}: The slope under the alternative hypothesis.
#'
#'     - \code{stDev}: The standard deviation of the residual.
#'
#'     - \code{stDevCovariate}: The standard deviation of the covariate.
#'
#' * \code{byStageResults}: A data frame containing the following variables:
#'
#'     - \code{informationRates}: The information rates.
#'
#'     - \code{efficacyBounds}: The efficacy boundaries on the Z-scale.
#'
#'     - \code{futilityBounds}: The futility boundaries on the Z-scale.
#'
#'     - \code{rejectPerStage}: The probability for efficacy stopping.
#'
#'     - \code{futilityPerStage}: The probability for futility stopping.
#'
#'     - \code{cumulativeRejection}: The cumulative probability for efficacy
#'       stopping.
#'
#'     - \code{cumulativeFutility}: The cumulative probability for futility
#'       stopping.
#'
#'     - \code{cumulativeAlphaSpent}: The cumulative alpha spent.
#'
#'     - \code{efficacyP}: The efficacy boundaries on the p-value scale.
#'
#'     - \code{futilityP}: The futility boundaries on the p-value scale.
#'
#'     - \code{information}: The cumulative information.
#'
#'     - \code{efficacyStopping}: Whether to allow efficacy stopping.
#'
#'     - \code{futilityStopping}: Whether to allow futility stopping.
#'
#'     - \code{rejectPerStageH0}: The probability for efficacy stopping
#'       under H0.
#'
#'     - \code{futilityPerStageH0}: The probability for futility stopping
#'       under H0.
#'
#'     - \code{cumulativeRejectionH0}: The cumulative probability for
#'       efficacy stopping under H0.
#'
#'     - \code{cumulativeFutilityH0}: The cumulative probability for futility
#'       stopping under H0.
#'
#'     - \code{efficacySlope}: The efficacy boundaries on the slope scale.
#'
#'     - \code{futilitySlope}: The futility boundaries on the slope scale.
#'
#'     - \code{numberOfSubjects}: The number of subjects.
#'
#' * \code{settings}: A list containing the following input parameters:
#'
#'     - \code{typeAlphaSpending}: The type of alpha spending.
#'
#'     - \code{parameterAlphaSpending}: The parameter value for alpha
#'       spending.
#'
#'     - \code{userAlphaSpending}: The user defined alpha spending.
#'
#'     - \code{typeBetaSpending}: The type of beta spending.
#'
#'     - \code{parameterBetaSpending}: The parameter value for beta spending.
#'
#'     - \code{userBetaSpending}: The user defined beta spending.
#'
#'     - \code{spendingTime}: The error spending time at each analysis.
#'
#'     - \code{normalApproximation}: The type of computation of the p-values.
#'       If \code{TRUE}, the variance is assumed to be known, otherwise
#'       the calculations are performed with the t distribution.
#'
#'     - \code{rounding}: Whether to round up sample size.
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @examples
#'
#' (design1 <- getDesignOneSlope(
#'   beta = 0.1, n = NA, slope = 0.5,
#'   stDev = 15, stDevCovariate = 9,
#'   normalApproximation = FALSE,
#'   alpha = 0.025))
#'
#' @export
getDesignOneSlope <- function(
    beta = NA_real_,
    n = NA_real_,
    slopeH0 = 0,
    slope = 0.5,
    stDev = 1,
    stDevCovariate = 1,
    normalApproximation = TRUE,
    rounding = TRUE,
    kMax = 1L,
    informationRates = NA_real_,
    efficacyStopping = NA_integer_,
    futilityStopping = NA_integer_,
    criticalValues = NA_real_,
    alpha = 0.025,
    typeAlphaSpending = "sfOF",
    parameterAlphaSpending = NA_real_,
    userAlphaSpending = NA_real_,
    futilityBounds = NA_real_,
    typeBetaSpending = "none",
    parameterBetaSpending = NA_real_,
    userBetaSpending = NA_real_,
    spendingTime = NA_real_) {

  if (is.na(beta) && is.na(n)) {
    stop("beta and n cannot be both missing")
  }

  if (!is.na(beta) && !is.na(n)) {
    stop("Only one of beta and n should be provided")
  }

  if (!is.na(beta) && (beta >= 1-alpha || beta < 0.0001)) {
    stop("beta must lie in [0.0001, 1-alpha)")
  }

  if (!is.na(n) && n <= 0) {
    stop("n must be positive")
  }

  if (stDev <= 0) {
    stop("stDev must be positive")
  }

  if (stDevCovariate <= 0) {
    stop("stDevCovariate must be positive")
  }

  if (alpha < 0.00001 || alpha >= 1) {
    stop("alpha must lie in [0.00001, 1)")
  }

  if (any(is.na(informationRates))) {
    informationRates = (1:kMax)/kMax
  }

  directionUpper = slope > slopeH0

  theta = ifelse(directionUpper, slope - slopeH0, slopeH0 - slope)

  # variance for one sampling unit
  v1 = stDev^2/stDevCovariate^2

  if (is.na(beta)) { # power calculation
    if (rounding) {
      n = ceiling(n - 1.0e-12)
      informationRates = round(n*informationRates)/n
    }

    des = getDesign(
      beta = NA, IMax = n/v1, theta,
      kMax, informationRates,
      efficacyStopping, futilityStopping,
      criticalValues, alpha, typeAlphaSpending,
      parameterAlphaSpending, userAlphaSpending,
      futilityBounds, typeBetaSpending,
      parameterBetaSpending, userBetaSpending,
      spendingTime, 1)

    if (kMax == 1 && !normalApproximation) { # t-test for fixed design
      if (!any(is.na(criticalValues))) {
        alpha = 1 - pnorm(criticalValues)
      }

      b = qt(1-alpha, n-2)
      ncp = theta*sqrt(n/v1)
      power = pt(b, n-2, ncp, lower.tail = FALSE)

      delta = b/sqrt(n/v1)

      des$overallResults$overallReject = power
      des$byStageResults$rejectPerStage = power
      des$byStageResults$futilityPerStage = 1 - power
      des$byStageResults$cumulativeRejection = power
      des$byStageResults$cumulativeFutility = 1 - power
      des$byStageResults$efficacyBounds = b
      des$byStageResults$futilityBounds = b

      if (directionUpper) {
        des$byStageResults$efficacySlope = delta + slopeH0
        des$byStageResults$futilitySlope = delta + slopeH0
      } else {
        des$byStageResults$efficacySlope = -delta + slopeH0
        des$byStageResults$futilitySlope = -delta + slopeH0
      }
    } else {
      if (directionUpper) {
        des$byStageResults$efficacySlope =
          des$byStageResults$efficacyTheta + slopeH0
        des$byStageResults$futilitySlope =
          des$byStageResults$futilityTheta + slopeH0
      } else {
        des$byStageResults$efficacySlope =
          -des$byStageResults$efficacyTheta + slopeH0
        des$byStageResults$futilitySlope =
          -des$byStageResults$futilityTheta + slopeH0
      }
    }
  } else { # sample size calculation
    if (kMax == 1 && !normalApproximation) { # t-test for fixed design
      if (!any(is.na(criticalValues))) {
        alpha = 1 - pnorm(criticalValues)
      }

      n0 = (qnorm(1-alpha) + qnorm(1-beta))^2*v1/theta^2

      n = uniroot(function(n) {
        b = qt(1-alpha, n-2)
        ncp = theta*sqrt(n/v1)
        pt(b, n-2, ncp, lower.tail = FALSE) - (1 - beta)
      }, c(n0, 2*n0))$root

      if (rounding) n = ceiling(n - 1.0e-12)

      des = getDesign(
        beta = NA, IMax = n/v1, theta,
        kMax, informationRates,
        efficacyStopping, futilityStopping,
        criticalValues, alpha, typeAlphaSpending,
        parameterAlphaSpending, userAlphaSpending,
        futilityBounds, typeBetaSpending,
        parameterBetaSpending, userBetaSpending,
        spendingTime, 1)

      b = qt(1-alpha, n-2)
      ncp = theta*sqrt(n/v1)
      power = pt(b, n-2, ncp, lower.tail = FALSE)

      delta = b/sqrt(n/v1)

      des$overallResults$overallReject = power
      des$byStageResults$rejectPerStage = power
      des$byStageResults$futilityPerStage = 1 - power
      des$byStageResults$cumulativeRejection = power
      des$byStageResults$cumulativeFutility = 1 - power
      des$byStageResults$efficacyBounds = b
      des$byStageResults$futilityBounds = b
      if (directionUpper) {
        des$byStageResults$efficacySlope = delta + slopeH0
        des$byStageResults$futilitySlope = delta + slopeH0
      } else {
        des$byStageResults$efficacySlope = -delta + slopeH0
        des$byStageResults$futilitySlope = -delta + slopeH0
      }
    } else {

      des = getDesign(
        beta, IMax = NA, theta,
        kMax, informationRates,
        efficacyStopping, futilityStopping,
        criticalValues, alpha, typeAlphaSpending,
        parameterAlphaSpending, userAlphaSpending,
        futilityBounds, typeBetaSpending,
        parameterBetaSpending, userBetaSpending,
        spendingTime, 1)

      n = des$overallResults$information*v1

      if (rounding) {
        n = ceiling(n - 1.0e-12)
        informationRates = des$byStageResults$informationRates
        informationRates = round(n*informationRates)/n

        des = getDesign(
          beta = NA, IMax = n/v1, theta,
          kMax, informationRates,
          efficacyStopping, futilityStopping,
          criticalValues, alpha, typeAlphaSpending,
          parameterAlphaSpending, userAlphaSpending,
          futilityBounds, typeBetaSpending,
          parameterBetaSpending, userBetaSpending,
          spendingTime, 1)
      }

      if (directionUpper) {
        des$byStageResults$efficacySlope =
          des$byStageResults$efficacyTheta + slopeH0
        des$byStageResults$futilitySlope =
          des$byStageResults$futilityTheta + slopeH0
      } else {
        des$byStageResults$efficacySlope =
          -des$byStageResults$efficacyTheta + slopeH0
        des$byStageResults$futilitySlope =
          -des$byStageResults$futilityTheta + slopeH0
      }
    }
  }

  des$overallResults$theta = theta
  des$overallResults$numberOfSubjects = n
  des$overallResults$expectedNumberOfSubjectsH1 =
    des$overallResults$expectedInformationH1*v1
  des$overallResults$expectedNumberOfSubjectsH0 =
    des$overallResults$expectedInformationH0*v1
  des$overallResults$slopeH0 = slopeH0
  des$overallResults$slope = slope
  des$overallResults$stDev = stDev
  des$overallResults$stDevCovariate = stDevCovariate

  des$byStageResults$efficacyTheta = NULL
  des$byStageResults$futilityTheta = NULL
  des$byStageResults$numberOfSubjects =
    des$byStageResults$informationRates*n

  des$settings$normalApproximation = normalApproximation
  des$settings$rounding = rounding
  des$settings$varianceRatio = NULL

  attr(des, "class") = "designOneSlope"

  des
}


#' @title Group Sequential Design for Two-Sample Slope Difference
#' @description Obtains the power given sample size or obtains the sample
#' size given power for a group sequential design for two-sample slope
#' difference.
#'
#' @param beta The type II error.
#' @param n The total sample size.
#' @param slopeDiffH0 The slope difference under the null hypothesis.
#'   Defaults to 0.
#' @param slopeDiff The slope difference under the alternative hypothesis.
#' @param stDev The standard deviation of the residual.
#' @param stDevCovariate The standard deviation of the covariate.
#' @param allocationRatioPlanned Allocation ratio for the active treatment
#'   versus control. Defaults to 1 for equal randomization.
#' @param normalApproximation The type of computation of the p-values.
#'   If \code{TRUE}, the variance is assumed to be known, otherwise
#'   the calculations are performed with the t distribution. The exact
#'   calculation using the t distribution is only implemented for the
#'   fixed design.
#' @param rounding Whether to round up sample size. Defaults to 1 for
#'   sample size rounding.
#' @inheritParams param_kMax
#' @param informationRates The information rates. Fixed prior to the trial.
#'   Defaults to \code{(1:kMax) / kMax} if left unspecified.
#' @inheritParams param_efficacyStopping
#' @inheritParams param_futilityStopping
#' @inheritParams param_criticalValues
#' @inheritParams param_alpha
#' @inheritParams param_typeAlphaSpending
#' @inheritParams param_parameterAlphaSpending
#' @inheritParams param_userAlphaSpending
#' @inheritParams param_futilityBounds
#' @inheritParams param_typeBetaSpending
#' @inheritParams param_parameterBetaSpending
#' @inheritParams param_userBetaSpending
#' @param spendingTime A vector of length \code{kMax} for the error spending
#'   time at each analysis. Defaults to missing, in which case, it is the
#'   same as \code{informationRates}.
#'
#' @details
#' We assume a simple linear regression of the form
#' \deqn{y_{g,i} = \alpha_g + \beta_g x_{g,i} + \epsilon_{g,i}}
#' for treatment group \eqn{g}, where \eqn{\epsilon_{g,i}}
#' is the residual error for subject \eqn{i} in group \eqn{g},
#' which is assumed to be normally distributed with mean 0 and
#' standard deviation \eqn{\sigma_\epsilon}.
#' The covariate \eqn{x_{g,i}} is assumed to
#' be normally distributed with mean 0 and standard deviation
#' \eqn{\sigma_x}. Since
#' \deqn{\hat{\beta}_g = \frac{\sum_{i=1}^{n_g} (x_{g,i}-\bar{x}_g) y_{g,i}}
#' {\sum_{i=1}^{n_g}(x_{g,i}-\bar{x}_g)^2}}
#' where \eqn{n_g} is the sample size for group \eqn{g},
#' it follows that
#' \deqn{\hat{\beta}_g \sim N(\beta_g,
#' \frac{\sigma_\epsilon^2}{\sum_{i=1}^{n_g}(x_{g,i}-\bar{x}_g)^2}).}
#' The slope difference is defined as
#' \deqn{\hat{\beta}_1 - \hat{\beta}_2}
#' where \eqn{\hat{\beta}_1} and \eqn{\hat{\beta}_2}
#' are the estimated slopes for treatment groups 1 and 2, respectively.
#' Since the variance of \eqn{\hat{\beta}_g} is
#' \deqn{\frac{\sigma_\epsilon^2}{n_g \sigma_x^2}}
#' we have
#' \deqn{\hat{\beta}_1 - \hat{\beta}_2 \sim N(\beta_1 - \beta_2,
#' (n_1^{-1} + n_2^{-1})\sigma_\epsilon^2/\sigma_x^2)}
#' which can be used to
#' calculate the power and sample size for the group sequential design.
#'
#' @return An S3 class \code{designSlopeDiff} object with three components:
#'
#' * \code{overallResults}: A data frame containing the following variables:
#'
#'     - \code{overallReject}: The overall rejection probability.
#'
#'     - \code{alpha}: The overall significance level.
#'
#'     - \code{attainedAlpha}: The attained significance level, which is
#'       different from the overall significance level in the presence of
#'       futility stopping.
#'
#'     - \code{kMax}: The number of stages.
#'
#'     - \code{theta}: The parameter value.
#'
#'     - \code{information}: The maximum information.
#'
#'     - \code{expectedInformationH1}: The expected information under H1.
#'
#'     - \code{expectedInformationH0}: The expected information under H0.
#'
#'     - \code{drift}: The drift parameter, equal to
#'       \code{theta*sqrt(information)}.
#'
#'     - \code{inflationFactor}: The inflation factor (relative to the
#'       fixed design).
#'
#'     - \code{numberOfSubjects}: The maximum number of subjects.
#'
#'     - \code{expectedNumberOfSubjectsH1}: The expected number of subjects
#'       under H1.
#'
#'     - \code{expectedNumberOfSubjectsH0}: The expected number of subjects
#'       under H0.
#'
#'     - \code{slopeDiffH0}: The slope difference under the null hypothesis.
#'
#'     - \code{slopeDiff}: The slope difference under the alternative
#'       hypothesis.
#'
#'     - \code{stDev}: The standard deviation of the residual.
#'
#'     - \code{stDevCovariate}: The standard deviation of the covariate.
#'
#' * \code{byStageResults}: A data frame containing the following variables:
#'
#'     - \code{informationRates}: The information rates.
#'
#'     - \code{efficacyBounds}: The efficacy boundaries on the Z-scale.
#'
#'     - \code{futilityBounds}: The futility boundaries on the Z-scale.
#'
#'     - \code{rejectPerStage}: The probability for efficacy stopping.
#'
#'     - \code{futilityPerStage}: The probability for futility stopping.
#'
#'     - \code{cumulativeRejection}: The cumulative probability for efficacy
#'       stopping.
#'
#'     - \code{cumulativeFutility}: The cumulative probability for futility
#'       stopping.
#'
#'     - \code{cumulativeAlphaSpent}: The cumulative alpha spent.
#'
#'     - \code{efficacyP}: The efficacy boundaries on the p-value scale.
#'
#'     - \code{futilityP}: The futility boundaries on the p-value scale.
#'
#'     - \code{information}: The cumulative information.
#'
#'     - \code{efficacyStopping}: Whether to allow efficacy stopping.
#'
#'     - \code{futilityStopping}: Whether to allow futility stopping.
#'
#'     - \code{rejectPerStageH0}: The probability for efficacy stopping
#'       under H0.
#'
#'     - \code{futilityPerStageH0}: The probability for futility stopping
#'       under H0.
#'
#'     - \code{cumulativeRejectionH0}: The cumulative probability for
#'       efficacy stopping under H0.
#'
#'     - \code{cumulativeFutilityH0}: The cumulative probability for futility
#'       stopping under H0.
#'
#'     - \code{efficacySlopeDiff}: The efficacy boundaries on the slope
#'       difference scale.
#'
#'     - \code{futilitySlopeDiff}: The futility boundaries on the slope
#'       difference scale.
#'
#'     - \code{numberOfSubjects}: The number of subjects.
#'
#' * \code{settings}: A list containing the following input parameters:
#'
#'     - \code{typeAlphaSpending}: The type of alpha spending.
#'
#'     - \code{parameterAlphaSpending}: The parameter value for alpha
#'       spending.
#'
#'     - \code{userAlphaSpending}: The user defined alpha spending.
#'
#'     - \code{typeBetaSpending}: The type of beta spending.
#'
#'     - \code{parameterBetaSpending}: The parameter value for beta spending.
#'
#'     - \code{userBetaSpending}: The user defined beta spending.
#'
#'     - \code{spendingTime}: The error spending time at each analysis.
#'
#'     - \code{allocationRatioPlanned}: Allocation ratio for the active
#'       treatment versus control.
#'
#'     - \code{normalApproximation}: The type of computation of the p-values.
#'       If \code{TRUE}, the variance is assumed to be known, otherwise
#'       the calculations are performed with the t distribution.
#'
#'     - \code{rounding}: Whether to round up sample size.
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @examples
#'
#' (design1 <- getDesignSlopeDiff(
#'   beta = 0.1, n = NA, slopeDiff = -0.5,
#'   stDev = 10, stDevCovariate = 6,
#'   normalApproximation = FALSE, alpha = 0.025))
#'
#' @export
getDesignSlopeDiff <- function(
    beta = NA_real_,
    n = NA_real_,
    slopeDiffH0 = 0,
    slopeDiff = 0.5,
    stDev = 1,
    stDevCovariate = 1,
    allocationRatioPlanned = 1,
    normalApproximation = TRUE,
    rounding = TRUE,
    kMax = 1L,
    informationRates = NA_real_,
    efficacyStopping = NA_integer_,
    futilityStopping = NA_integer_,
    criticalValues = NA_real_,
    alpha = 0.025,
    typeAlphaSpending = "sfOF",
    parameterAlphaSpending = NA_real_,
    userAlphaSpending = NA_real_,
    futilityBounds = NA_real_,
    typeBetaSpending = "none",
    parameterBetaSpending = NA_real_,
    userBetaSpending = NA_real_,
    spendingTime = NA_real_) {

  if (is.na(beta) && is.na(n)) {
    stop("beta and n cannot be both missing")
  }

  if (!is.na(beta) && !is.na(n)) {
    stop("Only one of beta and n should be provided")
  }

  if (!is.na(beta) && (beta >= 1-alpha || beta < 0.0001)) {
    stop("beta must lie in [0.0001, 1-alpha)")
  }

  if (!is.na(n) && n <= 0) {
    stop("n must be positive")
  }

  if (stDev <= 0) {
    stop("stDev must be positive")
  }

  if (allocationRatioPlanned <= 0) {
    stop("allocationRatioPlanned must be positive")
  }

  if (alpha < 0.00001 || alpha >= 1) {
    stop("alpha must lie in [0.00001, 1)")
  }

  if (any(is.na(informationRates))) {
    informationRates = (1:kMax)/kMax
  }


  r = allocationRatioPlanned/(1 + allocationRatioPlanned)

  directionUpper = slopeDiff > slopeDiffH0

  theta = ifelse(directionUpper, slopeDiff - slopeDiffH0,
                 slopeDiffH0 - slopeDiff)

  # variance for one sampling unit
  v1 = stDev^2/stDevCovariate^2/(r*(1-r))

  if (is.na(beta)) { # power calculation
    if (rounding) {
      n = ceiling(n - 1.0e-12)
      informationRates = round(n*informationRates)/n
    }

    des = getDesign(
      beta = NA, IMax = n/v1, theta,
      kMax, informationRates,
      efficacyStopping, futilityStopping,
      criticalValues, alpha, typeAlphaSpending,
      parameterAlphaSpending, userAlphaSpending,
      futilityBounds, typeBetaSpending,
      parameterBetaSpending, userBetaSpending,
      spendingTime, 1)

    if (kMax == 1 && !normalApproximation) { # t-test for fixed design
      if (!any(is.na(criticalValues))) {
        alpha = 1 - pnorm(criticalValues)
      }

      b = qt(1-alpha, n-4)
      ncp = theta*sqrt(n/v1)
      power = pt(b, n-4, ncp, lower.tail = FALSE)

      delta = b/sqrt(n/v1)

      des$overallResults$overallReject = power
      des$byStageResults$rejectPerStage = power
      des$byStageResults$futilityPerStage = 1 - power
      des$byStageResults$cumulativeRejection = power
      des$byStageResults$cumulativeFutility = 1 - power
      des$byStageResults$efficacyBounds = b
      des$byStageResults$futilityBounds = b

      if (directionUpper) {
        des$byStageResults$efficacySlopeDiff = delta + slopeDiffH0
        des$byStageResults$futilitySlopeDiff = delta + slopeDiffH0
      } else {
        des$byStageResults$efficacySlopeDiff = -delta + slopeDiffH0
        des$byStageResults$futilitySlopeDiff = -delta + slopeDiffH0
      }
    } else {
      if (directionUpper) {
        des$byStageResults$efficacySlopeDiff =
          des$byStageResults$efficacyTheta + slopeDiffH0
        des$byStageResults$futilitySlopeDiff =
          des$byStageResults$futilityTheta + slopeDiffH0
      } else {
        des$byStageResults$efficacySlopeDiff =
          -des$byStageResults$efficacyTheta + slopeDiffH0
        des$byStageResults$futilitySlopeDiff =
          -des$byStageResults$futilityTheta + slopeDiffH0
      }
    }
  } else { # sample size calculation
    if (kMax == 1 && !normalApproximation) { # t-test for fixed design
      if (!any(is.na(criticalValues))) {
        alpha = 1 - pnorm(criticalValues)
      }

      n0 = (qnorm(1-alpha) + qnorm(1-beta))^2*v1/theta^2

      n = uniroot(function(n) {
        b = qt(1-alpha, n-4)
        ncp = theta*sqrt(n/v1)
        pt(b, n-4, ncp, lower.tail = FALSE) - (1 - beta)
      }, c(n0, 2*n0))$root

      if (rounding) n = ceiling(n - 1.0e-12)

      des = getDesign(
        beta = NA, IMax = n/v1, theta,
        kMax, informationRates,
        efficacyStopping, futilityStopping,
        criticalValues, alpha, typeAlphaSpending,
        parameterAlphaSpending, userAlphaSpending,
        futilityBounds, typeBetaSpending,
        parameterBetaSpending, userBetaSpending,
        spendingTime, 1)

      b = qt(1-alpha, n-4)
      ncp = theta*sqrt(n/v1)
      power = pt(b, n-4, ncp, lower.tail = FALSE)

      delta = b/sqrt(n/v1)

      des$overallResults$overallReject = power
      des$byStageResults$rejectPerStage = power
      des$byStageResults$futilityPerStage = 1 - power
      des$byStageResults$cumulativeRejection = power
      des$byStageResults$cumulativeFutility = 1 - power
      des$byStageResults$efficacyBounds = b
      des$byStageResults$futilityBounds = b
      if (directionUpper) {
        des$byStageResults$efficacySlopeDiff = delta + slopeDiffH0
        des$byStageResults$futilitySlopeDiff = delta + slopeDiffH0
      } else {
        des$byStageResults$efficacySlopeDiff = -delta + slopeDiffH0
        des$byStageResults$futilitySlopeDiff = -delta + slopeDiffH0
      }
    } else {
      des = getDesign(
        beta, IMax = NA, theta,
        kMax, informationRates,
        efficacyStopping, futilityStopping,
        criticalValues, alpha, typeAlphaSpending,
        parameterAlphaSpending, userAlphaSpending,
        futilityBounds, typeBetaSpending,
        parameterBetaSpending, userBetaSpending,
        spendingTime, 1)

      n = des$overallResults$information*v1

      if (rounding) {
        n = ceiling(n - 1.0e-12)
        informationRates = des$byStageResults$informationRates
        informationRates = round(n*informationRates)/n

        des = getDesign(
          beta = NA, IMax = n/v1, theta,
          kMax, informationRates,
          efficacyStopping, futilityStopping,
          criticalValues, alpha, typeAlphaSpending,
          parameterAlphaSpending, userAlphaSpending,
          futilityBounds, typeBetaSpending,
          parameterBetaSpending, userBetaSpending,
          spendingTime, 1)
      }

      if (directionUpper) {
        des$byStageResults$efficacySlopeDiff =
          des$byStageResults$efficacyTheta + slopeDiffH0
        des$byStageResults$futilitySlopeDiff =
          des$byStageResults$futilityTheta + slopeDiffH0
      } else {
        des$byStageResults$efficacySlopeDiff =
          -des$byStageResults$efficacyTheta + slopeDiffH0
        des$byStageResults$futilitySlopeDiff =
          -des$byStageResults$futilityTheta + slopeDiffH0
      }
    }
  }

  des$overallResults$theta = theta
  des$overallResults$numberOfSubjects = n
  des$overallResults$expectedNumberOfSubjectsH1 =
    des$overallResults$expectedInformationH1*v1
  des$overallResults$expectedNumberOfSubjectsH0 =
    des$overallResults$expectedInformationH0*v1
  des$overallResults$slopeDiffH0 = slopeDiffH0
  des$overallResults$slopeDiff = slopeDiff
  des$overallResults$stDev = stDev
  des$overallResults$stDevCovariate = stDevCovariate

  des$byStageResults$efficacyTheta = NULL
  des$byStageResults$futilityTheta = NULL
  des$byStageResults$numberOfSubjects =
    des$byStageResults$informationRates*n

  des$settings$allocationRatioPlanned = allocationRatioPlanned
  des$settings$normalApproximation = normalApproximation
  des$settings$rounding = rounding
  des$settings$varianceRatio = NULL

  attr(des, "class") = "designSlopeDiff"

  des
}


#' @title Group Sequential Design for Two-Sample Slope Difference
#' From the MMRM Model
#' @description Obtains the power given sample size or obtains the sample
#' size given power for two-sample slope difference from the growth curve
#' MMRM model.
#'
#' @param beta The type II error.
#' @param slopeDiffH0 The slope difference under the null hypothesis.
#'   Defaults to 0.
#' @param slopeDiff The slope difference under the alternative
#'   hypothesis.
#' @param stDev The standard deviation of the residual.
#' @param stDevIntercept The standard deviation of the random intercept.
#' @param stDevSlope The standard deviation of the random slope.
#' @param corrInterceptSlope The correlation between the random
#'   intercept and random slope.
#' @param w The number of time units (e.g. weeks) per measurement visit
#'   in a period. In general, visits are more frequent in the beginning
#'   of the study and less frequent towards the end.
#' @param N The number of measurement visits in a period. For example,
#'   \code{w = c(8, 16)} and \code{N = c(2, Inf)} means that the
#'   response variable will be collected at baseline, week 8,
#'   week 16, and every 16 weeks thereafter.
#' @inheritParams param_accrualTime
#' @inheritParams param_accrualIntensity
#' @inheritParams param_piecewiseSurvivalTime
#' @inheritParams param_gamma1
#' @inheritParams param_gamma2
#' @inheritParams param_accrualDuration
#' @inheritParams param_followupTime
#' @inheritParams param_fixedFollowup
#' @inheritParams param_allocationRatioPlanned
#' @param normalApproximation The type of computation of the p-values.
#'   If \code{TRUE}, the variance is assumed to be known, otherwise
#'   the calculations are performed with the t distribution. The
#'   degrees of freedom for the t-distribution for testing the slope
#'   difference is calculated using the containment method, and
#'   is equal to the total number of observations minus two times
#'   the total number of subjects. The exact calculation using the
#'   t distribution is only implemented for the fixed design.
#' @param rounding Whether to round up sample size. Defaults to 1 for
#'   sample size rounding.
#' @inheritParams param_kMax
#' @param informationRates The information rates. Defaults to
#'   \code{(1:kMax) / kMax} if left unspecified.
#' @inheritParams param_efficacyStopping
#' @inheritParams param_futilityStopping
#' @inheritParams param_criticalValues
#' @inheritParams param_alpha
#' @inheritParams param_typeAlphaSpending
#' @inheritParams param_parameterAlphaSpending
#' @inheritParams param_userAlphaSpending
#' @inheritParams param_futilityBounds
#' @inheritParams param_typeBetaSpending
#' @inheritParams param_parameterBetaSpending
#' @inheritParams param_userBetaSpending
#' @param spendingTime A vector of length \code{kMax} for the error spending
#'   time at each analysis. Defaults to missing, in which case, it is the
#'   same as \code{informationRates}.
#'
#' @details
#'
#' We use the following random-effects model to compare two slopes:
#' \deqn{y_{ij} = \alpha + (\beta + \gamma x_i) t_j + a_i + b_i t_j
#' + e_{ij}} where
#'
#' * \eqn{\alpha}: overall intercept common across treatment groups
#'   due to randomization
#'
#' * \eqn{\beta}: slope for the control group
#'
#' * \eqn{\gamma}: difference in slopes between the active treatment and
#'   control groups
#'
#' * \eqn{x_i}: treatment indicator for subject \eqn{i},
#'   1 for the active treatment and 0 for the control
#'
#' * \eqn{t_j}: time point \eqn{j} for repeated measurements,
#'   \eqn{t_1 = 0 < t_2 < \ldots < t_k}
#'
#' * \eqn{(a_i, b_i)}: random intercept and random slope
#'   for subject \eqn{i}, \eqn{Var(a_i) = \sigma_a^2},
#'   \eqn{Var(b_i) = \sigma_b^2}, \eqn{Corr(a_i, b_i) = \rho}
#'
#' * \eqn{e_{ij}}: within-subject residual with variance \eqn{\sigma_e^2}
#'
#' By accounting for randomization, we improve the efficiency for
#' estimating the difference in slopes. The model also accommodates
#' unequally spaced time points and missing data.
#' Specifically, given a calendar time \eqn{\tau} for an interim or
#' final analysis, let \eqn{k} be the number of scheduled time points
#' up to and including \eqn{\tau}, subject to the follow-up duration for
#' fixed follow-up designs. Let the observed time points be
#' \eqn{t_1, t_2, \ldots, t_k}, where \eqn{t_1 = 0} denotes baseline.
#'
#' For a subject in treatment group \eqn{g} with observed data
#' pattern \eqn{j}, the design matrix for the fixed effects
#' \eqn{(\alpha, \beta, \gamma)'} is given by
#' \deqn{X_{g,j} = (\bm{1}_j, \vec{t}_j, I(g=1)\vec{t}_j)}
#' where \eqn{\bm{1}_j} is a \eqn{j}-vector of ones, and
#' \eqn{\vec{t}_j = (t_1,\ldots,t_j)'} is the column vector of
#' observed time points. The design matrix for the
#' random effects \eqn{(a_i, b_i)'} is \deqn{Z_j = (\bm{1}_j, \vec{t}_j)}
#' The variance-covariance matrix of the random effects is
#' \deqn{D = \left(\begin{array}{cc} \sigma_a^2 & \rho \sigma_a \sigma_b \\
#' \rho \sigma_a \sigma_b & \sigma_b^2 \end{array}\right)}
#' Therefore, the variance-covariance matrix for the observed data
#' for the subject is
#' \deqn{V_{j} = Z_j D Z_j' + \sigma_e^2 I_j}
#' where \eqn{I_j} is the \eqn{j\times j} identity matrix.
#' Let \eqn{\pi_g} denote the proportion of subjects in group \eqn{g}.
#' The information matrix for the fixed effects is \deqn{I = nJ}
#' where
#' \deqn{J = \sum_{g=1}^{2} \pi_g \sum_{j=1}^{k} p_{g,j}
#' X_{g,j}' V_j^{-1} X_{g,j}}
#' and \eqn{p_{g,j}} is the proportion of subjects in group \eqn{g}
#' with observed data pattern \eqn{j}.
#'
#' The variance of the estimator for the slope difference
#' \eqn{\hat{\gamma}} is given by
#' \deqn{\text{Var}(\hat{\gamma}) = \frac{1}{n} J^{-1}[3,3]}
#' which can be used to calculate the power and sample size for the
#' group sequential design to detect a slope difference.
#'
#' @return An S3 class \code{designSlopeDiffMMRM} object with three
#' components:
#'
#' * \code{overallResults}: A data frame containing the following variables:
#'
#'     - \code{overallReject}: The overall rejection probability.
#'
#'     - \code{alpha}: The overall significance level.
#'
#'     - \code{attainedAlpha}: The attained significance level, which is
#'       different from the overall significance level in the presence of
#'       futility stopping.
#'
#'     - \code{kMax}: The number of stages.
#'
#'     - \code{theta}: The parameter value.
#'
#'     - \code{information}: The maximum information.
#'
#'     - \code{expectedInformationH1}: The expected information under H1.
#'
#'     - \code{expectedInformationH0}: The expected information under H0.
#'
#'     - \code{drift}: The drift parameter, equal to
#'       \code{theta*sqrt(information)}.
#'
#'     - \code{inflationFactor}: The inflation factor (relative to the
#'       fixed design).
#'
#'     - \code{numberOfSubjects}: The maximum number of subjects.
#'
#'     - \code{studyDuration}: The maximum study duration.
#'
#'     - \code{expectedNumberOfSubjectsH1}: The expected number of subjects
#'       under H1.
#'
#'     - \code{expectedNumberOfSubjectsH0}: The expected number of subjects
#'       under H0.
#'
#'     - \code{expectedStudyDurationH1}: The expected study duration
#'       under H1.
#'
#'     - \code{expectedStudyDurationH0}: The expected study duration
#'       under H0.
#'
#'     - \code{accrualDuration}: The accrual duration.
#'
#'     - \code{followupTime}: The follow-up time.
#'
#'     - \code{fixedFollowup}: Whether a fixed follow-up design is used.
#'
#'     - \code{slopeDiffH0}: The slope difference under H0.
#'
#'     - \code{slopeDiff}: The slope difference under H1.
#'
#' * \code{byStageResults}: A data frame containing the following variables:
#'
#'     - \code{informationRates}: The information rates.
#'
#'     - \code{efficacyBounds}: The efficacy boundaries on the Z-scale.
#'
#'     - \code{futilityBounds}: The futility boundaries on the Z-scale.
#'
#'     - \code{rejectPerStage}: The probability for efficacy stopping.
#'
#'     - \code{futilityPerStage}: The probability for futility stopping.
#'
#'     - \code{cumulativeRejection}: The cumulative probability for efficacy
#'       stopping.
#'
#'     - \code{cumulativeFutility}: The cumulative probability for futility
#'       stopping.
#'
#'     - \code{cumulativeAlphaSpent}: The cumulative alpha spent.
#'
#'     - \code{efficacyP}: The efficacy boundaries on the p-value scale.
#'
#'     - \code{futilityP}: The futility boundaries on the p-value scale.
#'
#'     - \code{information}: The cumulative information.
#'
#'     - \code{efficacyStopping}: Whether to allow efficacy stopping.
#'
#'     - \code{futilityStopping}: Whether to allow futility stopping.
#'
#'     - \code{rejectPerStageH0}: The probability for efficacy stopping
#'       under H0.
#'
#'     - \code{futilityPerStageH0}: The probability for futility stopping
#'       under H0.
#'
#'     - \code{cumulativeRejectionH0}: The cumulative probability for
#'       efficacy stopping under H0.
#'
#'     - \code{cumulativeFutilityH0}: The cumulative probability for futility
#'       stopping under H0.
#'
#'     - \code{efficacySlopeDiff}: The efficacy boundaries on the slope
#'       difference scale.
#'
#'     - \code{futilitySlopeDiff}: The futility boundaries on the slope
#'       difference scale.
#'
#'     - \code{numberOfSubjects}: The number of subjects.
#'
#'     - \code{analysisTime}: The average time since trial start.
#'
#' * \code{settings}: A list containing the following input parameters:
#'
#'     - \code{typeAlphaSpending}: The type of alpha spending.
#'
#'     - \code{parameterAlphaSpending}: The parameter value for alpha
#'       spending.
#'
#'     - \code{userAlphaSpending}: The user defined alpha spending.
#'
#'     - \code{typeBetaSpending}: The type of beta spending.
#'
#'     - \code{parameterBetaSpending}: The parameter value for beta spending.
#'
#'     - \code{userBetaSpending}: The user defined beta spending.
#'
#'     - \code{spendingTime}: The error spending time at each analysis.
#'
#'     - \code{allocationRatioPlanned}: The allocation ratio for the active
#'       treatment versus control.
#'
#'     - \code{accrualTime}: A vector that specifies the starting time of
#'       piecewise Poisson enrollment time intervals.
#'
#'     - \code{accrualIntensity}: A vector of accrual intensities.
#'        One for each accrual time interval.
#'
#'     - \code{piecewiseSurvivalTime}: A vector that specifies the
#'       starting time of piecewise exponential survival time intervals.
#'
#'     - \code{gamma1}: The hazard rate for exponential dropout or
#'       a vector of hazard rates for piecewise exponential dropout
#'       for the active treatment group.
#'
#'     - \code{gamma2}: The hazard rate for exponential dropout or
#'       a vector of hazard rates for piecewise exponential dropout
#'       for the control group.
#'
#'     - \code{w}: The number of time units per measurement visit in a
#'       period.
#'
#'     - \code{N}: The number of measurement visits in a period.
#'
#'     - \code{stdDev}: The standard deviation of the residual.
#'
#'     - \code{G}: The covariance matrix for the random intercept and
#'       random slope.
#'
#'     - \code{normalApproximation}: The type of computation of the p-values.
#'       If \code{TRUE}, the variance is assumed to be known, otherwise
#'       the calculations are performed with the t distribution.
#'
#'     - \code{rounding}: Whether to round up sample size.
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @references
#'
#' Daniel O. Scharfstein, Anastasios A. Tsiatis, and James M. Robins.
#' Semiparametric efficiency and its implication on the
#' design and analysis of group-sequential studies.
#' Journal of the American Statistical Association 1997; 92:1342-1350.
#'
#' @examples
#'
#' (design1 <- getDesignSlopeDiffMMRM(
#'   beta = 0.2, slopeDiff = log(1.15)/52,
#'   stDev = sqrt(.182),
#'   stDevIntercept = sqrt(.238960),
#'   stDevSlope = sqrt(.000057),
#'   corrInterceptSlope = .003688/sqrt(.238960*.000057),
#'   w = 8,
#'   N = 10000,
#'   accrualIntensity = 15,
#'   gamma1 = 1/(4.48*52),
#'   gamma2 = 1/(4.48*52),
#'   accrualDuration = NA,
#'   followupTime = 8,
#'   alpha = 0.025))
#'
#' @export
getDesignSlopeDiffMMRM <- function(
    beta = NA_real_,
    slopeDiffH0 = 0,
    slopeDiff = 0.5,
    stDev = 1,
    stDevIntercept = 1,
    stDevSlope = 1,
    corrInterceptSlope = 0.5,
    w = NA_real_,
    N = NA_real_,
    accrualTime = 0,
    accrualIntensity = NA_real_,
    piecewiseSurvivalTime = 0,
    gamma1 = 0,
    gamma2 = 0,
    accrualDuration = NA_real_,
    followupTime = NA_real_,
    fixedFollowup = FALSE,
    allocationRatioPlanned = 1,
    normalApproximation = TRUE,
    rounding = TRUE,
    kMax = 1L,
    informationRates = NA_real_,
    efficacyStopping = NA_integer_,
    futilityStopping = NA_integer_,
    criticalValues = NA_real_,
    alpha = 0.025,
    typeAlphaSpending = "sfOF",
    parameterAlphaSpending = NA_real_,
    userAlphaSpending = NA_real_,
    futilityBounds = NA_real_,
    typeBetaSpending = "none",
    parameterBetaSpending = NA_real_,
    userBetaSpending = NA_real_,
    spendingTime = NA_real_) {

  m = length(w)             # total number of periods
  cumN = c(0, cumsum(N))    # cumulative number of visits
  cumwN = c(0, cumsum(w*N)) # cumulative number of weeks

  nintervals = length(piecewiseSurvivalTime)

  if (is.na(beta) + is.na(accrualDuration) + is.na(followupTime) > 1) {
    stop(paste("At most one can be missing among beta, accrualDuration,",
               "and followupTime"))
  }

  if (is.na(beta)) {
    unknown = "beta"
  } else if (is.na(accrualDuration)) {
    unknown = "accrualDuration"
  } else if (is.na(followupTime)) {
    unknown = "followupTime"
  } else {
    unknown = "accrualIntensity"
  }

  if (!is.na(beta) && (beta >= 1-alpha || beta < 0.0001)) {
    stop("beta must lie in [0.0001, 1-alpha)")
  }

  if (any(is.na(w))) {
    stop("w must be provided")
  }

  if (any(w <= 0)) {
    stop("Elements of w must be positive")
  }

  if (any(is.na(N))) {
    stop("N must be provided")
  }

  if (length(N) != m) {
    stop("w and N must have the same length")
  }

  if (any(N <= 0)) {
    stop("Elements of N must be positive")
  }

  if (sum(N) < 10000) {
    stop("The max # of measurements must be >=10000 to enable root finding")
  }

  if (stDev <= 0) {
    stop("stDev must be positive")
  }

  if (stDevIntercept <= 0) {
    stop("stDevIntercept must be positive")
  }

  if (stDevSlope <= 0) {
    stop("stDevSlope must be positive")
  }

  if (corrInterceptSlope <= -1 || corrInterceptSlope >= 1) {
    stop("corrInterceptSlope must lie between -1 and 1")
  }

  if (accrualTime[1] != 0) {
    stop("accrualTime must start with 0")
  }

  if (length(accrualTime) > 1 && any(diff(accrualTime) <= 0)) {
    stop("accrualTime should be increasing")
  }

  if (any(is.na(accrualIntensity))) {
    stop("accrualIntensity must be provided")
  }

  if (length(accrualTime) != length(accrualIntensity)) {
    stop("accrualTime and accrualIntensity must have the same length")
  }

  if (any(accrualIntensity < 0)) {
    stop("accrualIntensity must be non-negative")
  }

  if (piecewiseSurvivalTime[1] != 0) {
    stop("piecewiseSurvivalTime must start with 0")
  }

  if (nintervals > 1 && any(diff(piecewiseSurvivalTime) <= 0)) {
    stop("piecewiseSurvivalTime should be increasing")
  }

  if (any(gamma1 < 0)) {
    stop("gamma1 must be non-negative")
  }

  if (any(gamma2 < 0)) {
    stop("gamma2 must be non-negative")
  }

  if (length(gamma1) != 1 && length(gamma1) != nintervals) {
    stop("Invalid length for gamma1")
  }

  if (length(gamma2) != 1 && length(gamma2) != nintervals) {
    stop("Invalid length for gamma2")
  }

  if (length(gamma1) == 1) {
    gamma1 = rep(gamma1, nintervals)
  }

  if (length(gamma2) == 1) {
    gamma2 = rep(gamma2, nintervals)
  }

  if (!is.na(accrualDuration) && accrualDuration <= 0) {
    stop("accrualDuration must be positive")
  }

  if (!is.na(followupTime) && followupTime < 0) {
    stop("followupTime must be nonnegative")
  }

  if (fixedFollowup && is.na(followupTime)) {
    stop("followupTime must be provided for fixedFollowup")
  }

  if (fixedFollowup) {
    i = findInterval(followupTime, cumwN)
    if ((followupTime - cumwN[i]) %% w[i] != 0) {
      stop("followupTime must match a scheduled visit for fixedFollowup")
    }
  }

  if (allocationRatioPlanned <= 0) {
    stop("allocationRatioPlanned must be positive")
  }

  if (alpha < 0.00001 || alpha >= 1) {
    stop("alpha must lie in [0.00001, 1)")
  }

  if (any(is.na(informationRates))) {
    informationRates = (1:kMax)/kMax
  }


  r = allocationRatioPlanned/(1 + allocationRatioPlanned)

  directionUpper = slopeDiff > slopeDiffH0

  theta = ifelse(directionUpper, slopeDiff - slopeDiffH0,
                 slopeDiffH0 - slopeDiff)


  v11 = stDevIntercept^2
  v12 = corrInterceptSlope*stDevIntercept*stDevSlope
  v22 = stDevSlope^2

  G = matrix(c(v11, v12, v12, v22), 2, 2)

  # function to obtain the info for slope difference
  f_info <- function(tau, w, m, cumN, cumwN, stDev, G,
                     accrualTime, accrualIntensity,
                     piecewiseSurvivalTime, gamma1, gamma2,
                     accrualDuration, followupTime, fixedFollowup, r) {

    # total number of enrolled subjects at interim analysis
    n = accrual(tau, accrualTime, accrualIntensity, accrualDuration)

    i = findInterval(tau, cumwN)  # period for tau
    k = floor((tau - cumwN[i])/w[i]) + cumN[i] + 1  # interval for tau

    # interval indicators
    j = 1:k
    # period indicators for the visits at the start of the intervals
    i = pmax(findInterval(j-2, cumN), 1)
    # time points for the visits up to time tau
    # 0, w[1], ..., N[1]*w[1], N[1]*w[1]+w[2], ..., N[1]*w[1]+N[2]*w[2], ...
    t = cumwN[i] + (j - cumN[i] - 1)*w[i]
    if (fixedFollowup) t = pmin(t, followupTime)

    # total number of enrolled subjects at each time point
    ns = accrual(tau - t, accrualTime, accrualIntensity, accrualDuration)

    # probability of not dropping out at each time point by treatment
    q1 = ptpwexp(t, piecewiseSurvivalTime, gamma1, lower.tail = FALSE)
    q2 = ptpwexp(t, piecewiseSurvivalTime, gamma2, lower.tail = FALSE)

    # number of subjects remaining at each time point by treatment
    m1 = r*ns*q1
    m2 = (1-r)*ns*q2

    # number of subjects dropping out at each time point by treatment
    n1 = c(-diff(m1), m1[k])
    n2 = c(-diff(m2), m2[k])

    Z = matrix(c(rep(1,k), t), k, 2)
    covar = Z %*% G %*% t(Z) + stDev^2*diag(k)

    X1 = matrix(c(rep(1,k), t, t), k, 3)
    X2 = matrix(c(rep(1,k), t, rep(0,k)), k, 3)

    # information matrix per subject
    I1 = 0
    I2 = 0
    for (j in 1:k) {
      x1 = X1[1:j, , drop=FALSE]
      x2 = X2[1:j, , drop=FALSE]
      I = solve(covar[1:j, 1:j])
      I1 = I1 + n1[j]*t(x1) %*% I %*% x1
      I2 = I2 + n2[j]*t(x2) %*% I %*% x2
    }

    # variance for common intercept, control slope, and slope difference
    V = solve(I1 + I2)

    # information for slope difference
    IMax = 1/V[3,3]

    # degrees of freedom for slope difference
    nu = sum((n1 + n2)*(1:k)) - 2*n

    list(IMax = IMax, nu = nu)
  }


  # power calculation
  if (is.na(beta)) {
    n = accrual(accrualDuration, accrualTime, accrualIntensity,
                accrualDuration)

    if (rounding) {
      n = ceiling(n - 1.0e-12)
      accrualDuration = getAccrualDurationFromN(n, accrualTime,
                                                accrualIntensity)
    }

    studyDuration = accrualDuration + followupTime
    out = f_info(studyDuration, w, m, cumN, cumwN, stDev, G,
                 accrualTime, accrualIntensity,
                 piecewiseSurvivalTime, gamma1, gamma2,
                 accrualDuration, followupTime, fixedFollowup, r)

    IMax = out$IMax
    nu = out$nu

    des = getDesign(
      beta = NA, IMax, theta,
      kMax, informationRates,
      efficacyStopping, futilityStopping,
      criticalValues, alpha, typeAlphaSpending,
      parameterAlphaSpending, userAlphaSpending,
      futilityBounds, typeBetaSpending,
      parameterBetaSpending, userBetaSpending,
      spendingTime, 1)

    if (kMax == 1 && !normalApproximation) { # t-test for fixed design
      if (!any(is.na(criticalValues))) {
        alpha = 1 - pnorm(criticalValues)
      }

      b = qt(1-alpha, nu)
      ncp = theta*sqrt(IMax)
      power = pt(b, nu, ncp, lower.tail = FALSE)

      des$overallResults$overallReject = power
      des$byStageResults$rejectPerStage = power
      des$byStageResults$futilityPerStage = 1 - power
      des$byStageResults$cumulativeRejection = power
      des$byStageResults$cumulativeFutility = 1 - power
      des$byStageResults$efficacyBounds = b
      des$byStageResults$futilityBounds = b

      delta = b/sqrt(IMax)
      if (directionUpper) {
        des$byStageResults$efficacySlopeDiff = delta + slopeDiffH0
        des$byStageResults$futilitySlopeDiff = delta + slopeDiffH0
      } else {
        des$byStageResults$efficacySlopeDiff = -delta + slopeDiffH0
        des$byStageResults$futilitySlopeDiff = -delta + slopeDiffH0
      }
    } else {
      if (directionUpper) {
        des$byStageResults$efficacySlopeDiff =
          des$byStageResults$efficacyTheta + slopeDiffH0
        des$byStageResults$futilitySlopeDiff =
          des$byStageResults$futilityTheta + slopeDiffH0
      } else {
        des$byStageResults$efficacySlopeDiff =
          -des$byStageResults$efficacyTheta + slopeDiffH0
        des$byStageResults$futilitySlopeDiff =
          -des$byStageResults$futilityTheta + slopeDiffH0
      }
    }
  } else { # sample size calculation
    des = getDesign(
      beta, IMax = NA, theta,
      kMax, informationRates,
      efficacyStopping, futilityStopping,
      criticalValues, alpha, typeAlphaSpending,
      parameterAlphaSpending, userAlphaSpending,
      futilityBounds, typeBetaSpending,
      parameterBetaSpending, userBetaSpending,
      spendingTime, 1)

    IMax = des$overallResults$information

    if (kMax == 1 && !normalApproximation) { # t-test for fixed design
      if (unknown == "accrualDuration") {
        accrualDuration = uniroot(function(x) {
          studyDuration = x + followupTime
          out = f_info(studyDuration, w, m, cumN, cumwN, stDev, G,
                       accrualTime, accrualIntensity,
                       piecewiseSurvivalTime, gamma1, gamma2,
                       x, followupTime, fixedFollowup, r)
          out$IMax - IMax
        }, c(0.001, 240))$root
      } else if (unknown == "followupTime") {
        followupTime = uniroot(function(x) {
          studyDuration = accrualDuration + x
          out = f_info(studyDuration, w, m, cumN, cumwN, stDev, G,
                       accrualTime, accrualIntensity,
                       piecewiseSurvivalTime, gamma1, gamma2,
                       accrualDuration, followupTime, fixedFollowup, r)
          out$IMax - IMax
        }, c(0.001, 240))$root
      } else {
        accrualIntensity = uniroot(function(x) {
          studyDuration = accrualDuration + followupTime
          out = f_info(studyDuration, w, m, cumN, cumwN, stDev, G,
                       accrualTime, x*accrualIntensity,
                       piecewiseSurvivalTime, gamma1, gamma2,
                       accrualDuration, followupTime, fixedFollowup, r)
          out$IMax - IMax
        }, c(0.001, 240))$root*accrualIntensity
      }

      studyDuration = accrualDuration + followupTime
      out = f_info(studyDuration, w, m, cumN, cumwN, stDev, G,
                   accrualTime, accrualIntensity,
                   piecewiseSurvivalTime, gamma1, gamma2,
                   accrualDuration, followupTime, fixedFollowup, r)

      nu = out$nu # this serves as the lower bound for degrees of freedom

      if (!any(is.na(criticalValues))) {
        alpha = 1 - pnorm(criticalValues)
      }

      # power for t-test
      b = qt(1-alpha, nu)

      info = uniroot(function(info) {
        ncp = theta*sqrt(info)
        pt(b, nu, ncp, lower.tail = FALSE) - (1-beta)
      }, c(0.5*IMax, 1.5*IMax))$root


      if (unknown == "accrualDuration") {
        accrualDuration = uniroot(function(x) {
          studyDuration = x + followupTime
          out = f_info(studyDuration, w, m, cumN, cumwN, stDev, G,
                       accrualTime, accrualIntensity,
                       piecewiseSurvivalTime, gamma1, gamma2,
                       x, followupTime, fixedFollowup, r)
          out$IMax - info
        }, c(0.5*accrualDuration, 1.5*accrualDuration))$root
      } else if (unknown == "followupTime") {
        followupTime = uniroot(function(x) {
          studyDuration = accrualDuration + x
          out = f_info(studyDuration, w, m, cumN, cumwN, stDev, G,
                       accrualTime, accrualIntensity,
                       piecewiseSurvivalTime, gamma1, gamma2,
                       accrualDuration, followupTime, fixedFollowup, r)
          out$IMax - info
        }, c(0.5*followupTime, 1.5*followupTime))$root
      } else {
        accrualIntensity = uniroot(function(x) {
          studyDuration = accrualDuration + followupTime
          out = f_info(studyDuration, w, m, cumN, cumwN, stDev, G,
                       accrualTime, x*accrualIntensity,
                       piecewiseSurvivalTime, gamma1, gamma2,
                       accrualDuration, followupTime, fixedFollowup, r)
          out$IMax - info
        }, c(0.5, 1.5))$root*accrualIntensity
      }

      n = accrual(accrualDuration, accrualTime, accrualIntensity,
                  accrualDuration)

      if (rounding) {
        n = ceiling(n - 1.0e-12)

        if (unknown == "accrualDuration" || unknown == "followupTime") {
          accrualDuration = getAccrualDurationFromN(n, accrualTime,
                                                    accrualIntensity)
        } else {
          accrualIntensity = uniroot(function(x) {
            accrual(accrualDuration, accrualTime, x*accrualIntensity,
                    accrualDuration) - n
          }, c(0.5, 1.5))$root*accrualIntensity
        }
      }

      # update maximum information and degrees of freedom
      studyDuration = accrualDuration + followupTime
      out = f_info(studyDuration, w, m, cumN, cumwN, stDev, G,
                   accrualTime, accrualIntensity,
                   piecewiseSurvivalTime, gamma1, gamma2,
                   accrualDuration, followupTime, fixedFollowup, r)

      IMax = out$IMax
      nu = out$nu
      b = qt(1-alpha, nu)
      ncp = theta*sqrt(IMax)
      power = pt(b, nu, ncp, lower.tail = FALSE)

      des$overallResults$overallReject = power
      des$byStageResults$rejectPerStage = power
      des$byStageResults$futilityPerStage = 1 - power
      des$byStageResults$cumulativeRejection = power
      des$byStageResults$cumulativeFutility = 1 - power
      des$byStageResults$efficacyBounds = b
      des$byStageResults$futilityBounds = b

      delta = b/sqrt(IMax)
      if (directionUpper) {
        des$byStageResults$efficacySlopeDiff = delta + slopeDiffH0
        des$byStageResults$futilitySlopeDiff = delta + slopeDiffH0
      } else {
        des$byStageResults$efficacySlopeDiff = -delta + slopeDiffH0
        des$byStageResults$futilitySlopeDiff = -delta + slopeDiffH0
      }
    } else {
      if (unknown == "accrualDuration") {
        accrualDuration = uniroot(function(x) {
          studyDuration = x + followupTime
          out = f_info(studyDuration, w, m, cumN, cumwN, stDev, G,
                       accrualTime, accrualIntensity,
                       piecewiseSurvivalTime, gamma1, gamma2,
                       x, followupTime, fixedFollowup, r)
          out$IMax - IMax
        }, c(0.001, 240))$root
      } else if (unknown == "followupTime") {
        followupTime = uniroot(function(x) {
          studyDuration = accrualDuration + x
          out = f_info(studyDuration, w, m, cumN, cumwN, stDev, G,
                       accrualTime, accrualIntensity,
                       piecewiseSurvivalTime, gamma1, gamma2,
                       accrualDuration, followupTime, fixedFollowup, r)
          out$IMax - IMax
        }, c(0.001, 240))$root
      } else {
        accrualIntensity = uniroot(function(x) {
          studyDuration = accrualDuration + followupTime
          out = f_info(studyDuration, w, m, cumN, cumwN, stDev, G,
                       accrualTime, x*accrualIntensity,
                       piecewiseSurvivalTime, gamma1, gamma2,
                       accrualDuration, followupTime, fixedFollowup, r)
          out$IMax - IMax
        }, c(0.001, 240))$root*accrualIntensity
      }

      n = accrual(accrualDuration, accrualTime, accrualIntensity,
                  accrualDuration)

      if (rounding) {
        n = ceiling(n - 1.0e-12)

        if (unknown == "accrualDuration" || unknown == "followupTime") {
          accrualDuration = getAccrualDurationFromN(n, accrualTime,
                                                    accrualIntensity)
        } else {
          accrualIntensity = uniroot(function(x) {
            accrual(accrualDuration, accrualTime, x*accrualIntensity,
                    accrualDuration) - n
          }, c(0.5, 1.5))$root*accrualIntensity
        }


        # final maximum information
        studyDuration = accrualDuration + followupTime
        out = f_info(studyDuration, w, m, cumN, cumwN, stDev, G,
                     accrualTime, accrualIntensity,
                     piecewiseSurvivalTime, gamma1, gamma2,
                     accrualDuration, followupTime, fixedFollowup, r)

        IMax = out$IMax

        des = getDesign(
          beta = NA, IMax, theta,
          kMax, informationRates,
          efficacyStopping, futilityStopping,
          criticalValues, alpha, typeAlphaSpending,
          parameterAlphaSpending, userAlphaSpending,
          futilityBounds, typeBetaSpending,
          parameterBetaSpending, userBetaSpending,
          spendingTime, 1)
      }

      if (directionUpper) {
        des$byStageResults$efficacySlopeDiff =
          des$byStageResults$efficacyTheta + slopeDiffH0
        des$byStageResults$futilitySlopeDiff =
          des$byStageResults$futilityTheta + slopeDiffH0
      } else {
        des$byStageResults$efficacySlopeDiff =
          -des$byStageResults$efficacyTheta + slopeDiffH0
        des$byStageResults$futilitySlopeDiff =
          -des$byStageResults$futilityTheta + slopeDiffH0
      }
    }
  }


  # timing of interim analysis
  studyDuration = accrualDuration + followupTime
  information = IMax*informationRates
  analysisTime = rep(0, kMax)

  if (kMax > 1) {
    for (i in 1:(kMax-1)) {
      analysisTime[i] = uniroot(function(tau) {
        out = f_info(tau, w, m, cumN, cumwN, stDev, G,
                     accrualTime, accrualIntensity,
                     piecewiseSurvivalTime, gamma1, gamma2,
                     accrualDuration, followupTime, fixedFollowup, r)
        out$IMax - information[i]
      }, c(w[1]+0.001, studyDuration))$root # at least one postbaseline visit
    }
  }

  analysisTime[kMax] = studyDuration

  numberOfSubjects = accrual(analysisTime, accrualTime, accrualIntensity,
                             accrualDuration)

  p = des$byStageResults$rejectPerStage +
    des$byStageResults$futilityPerStage

  pH0 = des$byStageResults$rejectPerStageH0 +
    des$byStageResults$futilityPerStageH0

  des$overallResults$numberOfSubjects = n
  des$overallResults$studyDuration = studyDuration
  des$overallResults$expectedNumberOfSubjectsH1 = sum(p*numberOfSubjects)
  des$overallResults$expectedNumberOfSubjectsH0 = sum(pH0*numberOfSubjects)
  des$overallResults$expectedStudyDurationH1 = sum(p*analysisTime)
  des$overallResults$expectedStudyDurationH0 = sum(pH0*analysisTime)
  des$overallResults$accrualDuration = accrualDuration
  des$overallResults$followupTime = followupTime
  des$overallResults$fixedFollowup = fixedFollowup
  des$overallResults$slopeDiffH0 = slopeDiffH0
  des$overallResults$slopeDiff = slopeDiff

  des$byStageResults$numberOfSubjects = numberOfSubjects
  if (fixedFollowup) {
    numberOfCompleters = accrual(analysisTime - followupTime, accrualTime,
                                 accrualIntensity, accrualDuration)
    des$byStageResults$numberOfCompleters = numberOfCompleters
  }
  des$byStageResults$analysisTime = analysisTime
  des$byStageResults$efficacyTheta = NULL
  des$byStageResults$futilityTheta = NULL

  des$settings$allocationRatioPlanned = allocationRatioPlanned
  des$settings$accrualTime = accrualTime
  des$settings$accrualIntensity = accrualIntensity
  des$settings$piecewiseSurvivalTime = piecewiseSurvivalTime
  des$settings$gamma1 = gamma1
  des$settings$gamma2 = gamma2
  des$settings$w = w
  des$settings$N = N
  des$settings$stDev = stDev
  des$settings$G = G
  des$settings$normalApproximation = normalApproximation
  des$settings$rounding = rounding
  des$settings$varianceRatio = NULL

  attr(des, "class") = "designSlopeDiffMMRM"

  des
}


#' @title Hedges' g Effect Size
#' @description Obtains Hedges' g estimate and confidence interval of
#' effect size.
#'
#' @param tstat The value of the t-test statistic for comparing two
#'   treatment conditions.
#' @param m The degrees of freedom for the t-test.
#' @param ntilde The normalizing sample size to convert the
#'   standardized treatment difference to the t-test statistic, i.e.,
#'   \code{tstat = sqrt(ntilde)*meanDiff/stDev}.
#' @param cilevel The confidence interval level. Defaults to 0.95.
#'
#' @details
#'
#' Hedges' \eqn{g} is an effect size measure commonly used in meta-analysis
#' to quantify the difference between two groups. It's an improvement
#' over Cohen's \eqn{d}, particularly when dealing with small sample sizes.
#'
#' The formula for Hedges' \eqn{g} is \deqn{g = c(m) d} where \eqn{d}
#' is Cohen's \eqn{d} effect size estimate, and \eqn{c(m)} is the bias
#' correction factor, \deqn{d = (\hat{\mu}_1 - \hat{\mu}_2)/\hat{\sigma}}
#' \deqn{c(m) = 1 - \frac{3}{4m-1}.}
#' Since \eqn{c(m) < 1}, Cohen's \eqn{d} overestimates the true effect size,
#' \eqn{\delta = (\mu_1 - \mu_2)/\sigma.}
#' Since \deqn{t = \sqrt{\tilde{n}} d} we have
#' \deqn{g = \frac{c(m)}{\sqrt{\tilde{n}}} t} where \eqn{t}
#' has a noncentral \eqn{t} distribution with \eqn{m} degrees of freedom
#' and noncentrality parameter \eqn{\sqrt{\tilde{n}} \delta}.
#'
#' The asymptotic variance of \eqn{g} can be approximated by
#' \deqn{Var(g) = \frac{1}{\tilde{n}} + \frac{g^2}{2m}.}
#' The confidence interval for \eqn{\delta}
#' can be constructed using normal approximation.
#'
#' For two-sample mean difference with sample size \eqn{n_1} for the
#' treatment group and \eqn{n_2} for the control group, we have
#' \eqn{\tilde{n} = \frac{n_1n_2}{n_1+n_2}} and \eqn{m=n_1+n_2-2}
#' for pooled variance estimate.
#'
#' @return A data frame with the following variables:
#'
#' * \code{tstat}: The value of the \code{t} test statistic.
#'
#' * \code{m}: The degrees of freedom for the t-test.
#'
#' * \code{ntilde}: The normalizing sample size to convert the
#'   standardized treatment difference to the t-test statistic.
#'
#' * \code{g}: Hedges' \code{g} effect size estimate.
#'
#' * \code{varg}: Variance of \code{g}.
#'
#' * \code{lower}: The lower confidence limit for effect size.
#'
#' * \code{upper}: The upper confidence limit for effect size.
#'
#' * \code{cilevel}: The confidence interval level.
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @references
#'
#' Larry V. Hedges. Distribution theory for Glass's estimator of
#' effect size and related estimators.
#' Journal of Educational Statistics 1981; 6:107-128.
#'
#' @examples
#'
#' n1 = 7
#' n2 = 8
#' meanDiff = 0.444
#' stDev = 1.201
#' m = n1+n2-2
#' ntilde = n1*n2/(n1+n2)
#' tstat = sqrt(ntilde)*meanDiff/stDev
#'
#' hedgesg(tstat, m, ntilde)
#'
#' @export
hedgesg <- function(tstat, m, ntilde, cilevel = 0.95) {
  d = 1/sqrt(ntilde)*tstat
  g = (1 - 3/(4*m-1))*d
  varg = 1/ntilde + g^2/(2*m)
  lower = g - qnorm((1 + cilevel)/2)*sqrt(varg)
  upper = g + qnorm((1 + cilevel)/2)*sqrt(varg)
  data.frame(tstat, m, ntilde, g, varg, lower, upper, cilevel)
}


