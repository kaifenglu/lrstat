#' @title Group sequential design for one-sample proportion
#' @description Obtains the power given sample size or obtains the sample
#' size given power for a group sequential design for one-sample proportion.
#'
#' @param beta The type II error.
#' @param n The total sample size.
#' @param piH0 The response probability under the null hypothesis.
#' @param pi The response probability under the alternative hypothesis.
#' @param normalApproximation The type of computation of the p-values.
#'   If \code{TRUE}, the normal approximation will be used, otherwise
#'   the calculations are performed with the binomial distribution. The exact
#'   calculation using the binomial distribution is only implemented for the
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
#' @return An S3 class \code{designOneProportion} object with three
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
#'       futility stopping as well as for the binomial exact test in a fixed
#'       design.
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
#'     - \code{piH0}: The response probability under the null hypothesis.
#'
#'     - \code{pi}: The response probability under the alternative
#'       hypothesis.
#'
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
#'     - \code{efficacyResponses}: The efficacy boundaries on the number
#'       of responses scale.
#'
#'     - \code{futilityResponses}: The futility boundaries on the number
#'       of responses scale.
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
#'     - \code{calculationTarget}: The calculation target, \code{beta} or
#'       \code{n}.
#'
#'     - \code{normalApproximation}: The type of computation of the p-values.
#'       If \code{TRUE}, the variance is assumed to be known, otherwise
#'       the calculations are performed with the binomial distribution.
#'
#'     - \code{rounding}: Whether to round up sample size.
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @examples
#'
#' # Example 1: group sequential trial power calculation
#' (design1 <- getDesignOneProportion(
#'   beta = 0.2, n = NA, piH0 = 0.15, pi = 0.25,
#'   kMax = 3, alpha = 0.05, typeAlphaSpending = "sfOF"))
#'
#' # Example 2: sample size calculation for one-sample binomial exact test
#' (design2 <- getDesignOneProportion(
#'   beta = 0.2, n = NA, piH0 = 0.15, pi = 0.25,
#'   normalApproximation = FALSE, alpha = 0.05))
#'
#' @export
getDesignOneProportion <- function(
    beta = NA_real_,
    n = NA_real_,
    piH0 = 0.1,
    pi = 0.2,
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
    stop("beta must lie in [0.0001, 1-alpha)");
  }

  if (!is.na(n) && n <= 0) {
    stop("n must be positive")
  }

  if (piH0 <= 0 || piH0 >= 1) {
    stop("piH0 must lie between 0 and 1")
  }

  if (pi <= 0 || pi >= 1) {
    stop("pi must lie between 0 and 1")
  }

  if (alpha < 0.00001 || alpha >= 1) {
    stop("alpha must lie in [0.00001, 1)")
  }

  if (any(is.na(informationRates))) {
    informationRates = (1:kMax)/kMax
  }

  directionUpper = pi > piH0

  theta = ifelse(directionUpper, pi - piH0, piH0 - pi)
  v1 = piH0*(1-piH0)

  if (!is.na(n)) { # power calculation
    if (rounding) {
      n = ceiling(n)
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
      spendingTime)

    if (kMax == 1 && !normalApproximation) { # exact test for fixed design
      if (!any(is.na(criticalValues))) {
        alpha = 1 - pnorm(criticalValues)
      }

      a = powerOnePropExact(n, piH0, pi, alpha)
      r = a$r
      attainedAlpha = a$attainedAlpha
      power = a$power
      n = a$n

      if (directionUpper) {
        b = (r/n - piH0)*sqrt(n/v1)
      } else {
        b = (piH0 - r/n)*sqrt(n/v1)
      }

      des$overallResults$overallReject = power
      des$overallResults$attainedAlpha = attainedAlpha

      des$byStageResults$rejectPerStage = power
      des$byStageResults$futilityPerStage = 1 - power
      des$byStageResults$cumulativeRejection = power
      des$byStageResults$cumulativeFutility = 1 - power
      des$byStageResults$cumulativeAlphaSpent = attainedAlpha
      des$byStageResults$efficacyBounds = b
      des$byStageResults$futilityBounds = b
      des$byStageResults$efficacyP = attainedAlpha
      des$byStageResults$futilityP = attainedAlpha
      des$byStageResults$efficacyResponses = r
      des$byStageResults$futilityResponses = r
      des$byStageResults$rejectPerStageH0 = attainedAlpha
      des$byStageResults$futilityPerStageH0 = 1 - attainedAlpha
      des$byStageResults$cumulativeRejectionH0 = attainedAlpha
      des$byStageResults$cumulativeFutilityH0 = 1 - attainedAlpha
    } else {
      ns = des$byStageResults$informationRates*n

      if (directionUpper) {
        des$byStageResults$efficacyResponses =
          ceiling(ns*(des$byStageResults$efficacyTheta + piH0))
        des$byStageResults$futilityResponses =
          ceiling(ns*(des$byStageResults$futilityTheta + piH0))
      } else {
        des$byStageResults$efficacyResponses =
          floor(ns*(-des$byStageResults$efficacyTheta + piH0))
        des$byStageResults$futilityResponses =
          floor(ns*(-des$byStageResults$futilityTheta + piH0))
      }
    }
  } else { # sample size calculation
    if (kMax == 1 && !normalApproximation) { # exact test for fixed design
      if (!any(is.na(criticalValues))) {
        alpha = 1 - pnorm(criticalValues)
      }

      a = samplesizeOnePropExact(beta, piH0, pi, alpha)
      attainedAlpha = a$attainedAlpha
      r = a$r
      power = a$power
      n = a$n

      if (directionUpper) {
        b = (r/n - piH0)*sqrt(n/v1)
      } else {
        b = (piH0 - r/n)*sqrt(n/v1)
      }

      des = getDesign(
        beta = NA, IMax = n/v1, theta,
        kMax, informationRates,
        efficacyStopping, futilityStopping,
        criticalValues, alpha, typeAlphaSpending,
        parameterAlphaSpending, userAlphaSpending,
        futilityBounds, typeBetaSpending,
        parameterBetaSpending, userBetaSpending,
        spendingTime)

      des$overallResults$overallReject = power
      des$overallResults$attainedAlpha = attainedAlpha

      des$byStageResults$rejectPerStage = power
      des$byStageResults$futilityPerStage = 1 - power
      des$byStageResults$cumulativeRejection = power
      des$byStageResults$cumulativeFutility = 1 - power
      des$byStageResults$cumulativeAlphaSpent = attainedAlpha
      des$byStageResults$efficacyBounds = b
      des$byStageResults$futilityBounds = b
      des$byStageResults$efficacyP = attainedAlpha
      des$byStageResults$futilityP = attainedAlpha
      des$byStageResults$efficacyResponses = r
      des$byStageResults$futilityResponses = r
      des$byStageResults$rejectPerStageH0 = attainedAlpha
      des$byStageResults$futilityPerStageH0 = 1 - attainedAlpha
      des$byStageResults$cumulativeRejectionH0 = attainedAlpha
      des$byStageResults$cumulativeFutilityH0 = 1 - attainedAlpha
    } else {
      des = getDesign(
        beta, IMax = NA, theta,
        kMax, informationRates,
        efficacyStopping, futilityStopping,
        criticalValues, alpha, typeAlphaSpending,
        parameterAlphaSpending, userAlphaSpending,
        futilityBounds, typeBetaSpending,
        parameterBetaSpending, userBetaSpending,
        spendingTime)

      n = des$overallResults$information*v1

      if (rounding) {
        n = ceiling(n)
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
          spendingTime)
      }

      ns = des$byStageResults$informationRates*n

      if (directionUpper) {
        des$byStageResults$efficacyResponses =
          ceiling(ns*(des$byStageResults$efficacyTheta + piH0))
        des$byStageResults$futilityResponses =
          ceiling(ns*(des$byStageResults$futilityTheta + piH0))
      } else {
        des$byStageResults$efficacyResponses =
          floor(ns*(-des$byStageResults$efficacyTheta + piH0))
        des$byStageResults$futilityResponses =
          floor(ns*(-des$byStageResults$futilityTheta + piH0))
      }
    }
  }

  des$overallResults$theta = theta
  des$overallResults$numberOfSubjects = n
  des$overallResults$expectedNumberOfSubjectsH1 =
    des$overallResults$expectedInformationH1*v1
  des$overallResults$expectedNumberOfSubjectsH0 =
    des$overallResults$expectedInformationH0*v1
  des$overallResults$piH0 = piH0
  des$overallResults$pi = pi

  des$byStageResults$efficacyTheta = NULL
  des$byStageResults$futilityTheta = NULL
  des$byStageResults$numberOfSubjects =
    des$byStageResults$informationRates*n

  des$settings$calculationTarget = ifelse(is.na(beta), "beta", "n")
  des$settings$normalApproximation = normalApproximation
  des$settings$rounding = rounding

  attr(des, "class") = "designOneProportion"

  des
}


#' @title Group sequential design for McNemar's test for paired
#' proportions
#' @description Obtains the power given sample size or obtains the sample
#' size given power for a group sequential design for McNemar's test for
#' paired proportions.
#'
#' @param beta The type II error.
#' @param n The total sample size.
#' @param pDiscordant The proportion of discordant pairs
#'   (xi = pi01 + pi10).
#' @param riskDiff The risk difference between the active and control
#'   treatments (delta = pi_t - pi_c = pi01 - pi10)
#' @param nullVariance Whether to use the variance under the null
#'   or the variance under the alternative.
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
#' @return An S3 class \code{designPairedPropMcNemar} object with three
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
#'     - \code{pDiscordant}: The proportion of discordant pairs
#'       (xi = pi01 + pi10).
#'
#'     - \code{riskDiff}: The risk difference between the active and control
#'       treatments (delta = pi_t - pi_c = pi01 - pi10)
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
#'     - \code{efficacyRiskDiff}: The efficacy boundaries on the risk
#'       difference scale.
#'
#'     - \code{futilityRiskDiff}: The futility boundaries on the risk
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
#'     - \code{varianceRatio}: The ratio of the variance under H0 to the
#'       variance under H1.
#'
#'     - \code{calculationTarget}: The calculation target, \code{beta} or
#'       \code{n}.
#'
#'     - \code{rounding}: Whether to round up sample size.
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @examples
#'
#' # Example 1: fixed design
#' (design1 <- getDesignPairedPropMcNemar(
#'   beta = 0.1, n = NA, pDiscordant = 0.16, riskDiff = 0.1,
#'   alpha = 0.025))
#'
#' # Example 2: group sequential design
#' (design2 <- getDesignPairedPropMcNemar(
#'   beta = 0.1, n = NA, pDiscordant = 0.16, riskDiff = 0.1,
#'   alpha = 0.025, kMax = 3, typeAlphaSpending = "sfOF"))
#'
#' @export
getDesignPairedPropMcNemar <- function(
    beta = NA_real_,
    n = NA_real_,
    pDiscordant = NA_real_,
    riskDiff = NA_real_,
    nullVariance = TRUE,
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
    stop("beta must lie in [0.0001, 1-alpha)");
  }

  if (!is.na(n) && n <= 0) {
    stop("n must be positive")
  }

  if (is.na(pDiscordant)) {
    stop("pDiscordant must be provided")
  }

  if (pDiscordant <= 0 || pDiscordant >= 1) {
    stop("pDiscordant must lie between 0 and 1")
  }

  if (is.na(riskDiff)) {
    stop("riskDiff must be provided")
  }

  if (riskDiff <= -1 || riskDiff >= 1) {
    stop("riskDiff must lie between -1 and 1")
  }

  if (alpha < 0.00001 || alpha >= 1) {
    stop("alpha must lie in [0.00001, 1)")
  }

  if (any(is.na(informationRates))) {
    informationRates = (1:kMax)/kMax
  }

  directionUpper = riskDiff > 0

  theta = ifelse(directionUpper, riskDiff, -riskDiff)

  # variance for one sampling unit
  v0 = pDiscordant
  v1 = pDiscordant - riskDiff^2
  varianceRatio = ifelse(nullVariance, v0/v1, 1)

  if (!is.na(n)) { # power calculation
    if (rounding) {
      n = ceiling(n)
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
      spendingTime, varianceRatio)
  } else { # sample size calculation
    des = getDesign(
      beta, IMax = NA, theta,
      kMax, informationRates,
      efficacyStopping, futilityStopping,
      criticalValues, alpha, typeAlphaSpending,
      parameterAlphaSpending, userAlphaSpending,
      futilityBounds, typeBetaSpending,
      parameterBetaSpending, userBetaSpending,
      spendingTime, varianceRatio)

    n = des$overallResults$information*v1

    if (rounding) {
      n = ceiling(n)
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
        spendingTime, varianceRatio)
    }
  }

  des$overallResults$numberOfSubjects = n
  des$overallResults$expectedNumberOfSubjectsH1 =
    des$overallResults$expectedInformationH1*v1
  des$overallResults$expectedNumberOfSubjectsH0 =
    des$overallResults$expectedInformationH0*v1
  des$overallResults$pDiscordant = pDiscordant
  des$overallResults$riskDiff = riskDiff

  if (directionUpper) {
    des$byStageResults$efficacyRiskDiff =
      des$byStageResults$efficacyTheta
    des$byStageResults$futilityRiskDiff =
      des$byStageResults$futilityTheta
  } else {
    des$byStageResults$efficacyRiskDiff =
      -des$byStageResults$efficacyTheta
    des$byStageResults$futilityRiskDiff =
      -des$byStageResults$futilityTheta
  }
  des$byStageResults$efficacyTheta = NULL
  des$byStageResults$futilityTheta = NULL
  des$byStageResults$numberOfSubjects =
    des$byStageResults$informationRates*n

  des$settings$calculationTarget = ifelse(is.na(beta), "beta", "n")
  des$settings$nullVariance = nullVariance
  des$settings$rounding = rounding

  attr(des, "class") = "designPairedPropMcNemar"

  des
}


#' @title Group sequential design for two-sample risk difference
#' @description Obtains the power given sample size or obtains the sample
#' size given power for a group sequential design for two-sample risk
#' difference.
#'
#' @param beta The type II error.
#' @param n The total sample size.
#' @param riskDiffH0 The risk difference under the null hypothesis.
#'   Defaults to 0.
#' @param pi1 The assumed probability for the active treatment group.
#' @param pi2 The assumed probability for the control group.
#' @param nullVariance Whether to use the variance under the null or
#'   the empirical variance under the alternative.
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
#' @return An S3 class \code{designRiskDiff} object with three
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
#'     - \code{riskDiffH0}: The risk difference under the null hypothesis.
#'
#'     - \code{pi1}: The assumed probability for the active treatment group.
#'
#'     - \code{pi2}: The assumed probability for the control group.
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
#'     - \code{efficacyRiskDiff}: The efficacy boundaries on the risk
#'       difference scale.
#'
#'     - \code{futilityRiskDiff}: The futility boundaries on the risk
#'       difference scale.
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
#'     - \code{varianceRatio}: The ratio of the variance under H0 to
#'       the variance under H1.
#'
#'     - \code{calculationTarget}: The calculation target, \code{beta} or
#'       \code{n}.
#'
#'     - \code{nullVariance}: Whether to use the variance under the null or
#'       the empirical variance under the alternative.
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
#' (design1 <- getDesignRiskDiff(
#'   beta = 0.2, n = NA, pi1 = 0.1, pi2 = 0.15,
#'   kMax = 3, alpha = 0.025, typeAlphaSpending = "sfOF",
#'   nullVariance = 0))
#'
#' @export
getDesignRiskDiff <- function(
    beta = NA_real_,
    n = NA_real_,
    riskDiffH0 = 0,
    pi1 = NA_real_,
    pi2 = NA_real_,
    nullVariance = TRUE,
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
    stop("beta must lie in [0.0001, 1-alpha)");
  }

  if (!is.na(n) && n <= 0) {
    stop("n must be positive")
  }

  if (riskDiffH0 <= -1 || riskDiffH0 >= 1) {
    stop("riskDiffH0 must lie between -1 and 1")
  }

  if (pi1 <= 0 || pi1 >= 1) {
    stop("pi1 must lie between 0 and 1")
  }

  if (pi2 <= 0 || pi2 >= 1) {
    stop("pi2 must lie between 0 and 1")
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

  riskDiff = pi1 - pi2

  directionUpper = riskDiff > riskDiffH0

  theta = ifelse(directionUpper, riskDiff - riskDiffH0,
                 riskDiffH0 - riskDiff)

  # restricted maximum likelihood estimates
  mr = remlRiskDiff(riskDiffH0, r, r*pi1, 1-r, (1-r)*pi2)
  p1 = mr[1]
  p2 = mr[2]

  # variance for one sampling unit
  v0 = p1*(1-p1)/r + p2*(1-p2)/(1-r)
  v1 = pi1*(1-pi1)/r + pi2*(1-pi2)/(1-r)

  varianceRatio = ifelse(nullVariance, v0/v1, 1)

  if (!is.na(n)) { # power calculation
    if (rounding) {
      n = ceiling(n)
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
      spendingTime, varianceRatio)
  } else { # sample size calculation
    des = getDesign(
      beta, IMax = NA, theta,
      kMax, informationRates,
      efficacyStopping, futilityStopping,
      criticalValues, alpha, typeAlphaSpending,
      parameterAlphaSpending, userAlphaSpending,
      futilityBounds, typeBetaSpending,
      parameterBetaSpending, userBetaSpending,
      spendingTime, varianceRatio)

    n = des$overallResults$information*v1

    if (rounding) {
      n = ceiling(n)
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
        spendingTime, varianceRatio)
    }
  }

  des$overallResults$theta = theta
  des$overallResults$numberOfSubjects = n
  des$overallResults$expectedNumberOfSubjectsH1 =
    des$overallResults$expectedInformationH1*v1
  des$overallResults$expectedNumberOfSubjectsH0 =
    des$overallResults$expectedInformationH0*v1
  des$overallResults$riskDiffH0 = riskDiffH0
  des$overallResults$pi1 = pi1
  des$overallResults$pi2 = pi2

  if (directionUpper) {
    des$byStageResults$efficacyRiskDiff =
      des$byStageResults$efficacyTheta + riskDiffH0
    des$byStageResults$futilityRiskDiff =
      des$byStageResults$futilityTheta + riskDiffH0
  } else {
    des$byStageResults$efficacyRiskDiff =
      -des$byStageResults$efficacyTheta + riskDiffH0
    des$byStageResults$futilityRiskDiff =
      -des$byStageResults$futilityTheta + riskDiffH0
  }

  des$byStageResults$efficacyTheta = NULL
  des$byStageResults$futilityTheta = NULL
  des$byStageResults$numberOfSubjects =
    des$byStageResults$informationRates*n

  des$settings$calculationTarget = ifelse(is.na(beta), "beta", "n")
  des$settings$nullVariance = nullVariance
  des$settings$varianceRatio = varianceRatio
  des$settings$allocationRatioPlanned = allocationRatioPlanned
  des$settings$rounding = rounding

  attr(des, "class") = "designRiskDiff"

  des
}


#' @title Group sequential design for two-sample risk ratio
#' @description Obtains the power given sample size or obtains the sample
#' size given power for a group sequential design for two-sample risk
#' ratio.
#'
#' @param beta The type II error.
#' @param n The total sample size.
#' @param riskRatioH0 The risk ratio under the null hypothesis.
#'   Defaults to 1.
#' @param pi1 The assumed probability for the active treatment group.
#' @param pi2 The assumed probability for the control group.
#' @param nullVariance Whether to use the variance under the null or
#'   the empirical variance under the alternative.
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
#' @return An S3 class \code{designRiskRatio} object with three components:
#'
#' * \code{overallResults}: A data frame containing the following variables:
#'
#'    - \code{overallReject}: The overall rejection probability.
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
#'     - \code{riskRatioH0}: The risk ratio under the null hypothesis.
#'
#'     - \code{pi1}: The assumed probability for the active treatment group.
#'
#'     - \code{pi2}: The assumed probability for the control group.
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
#'     - \code{efficacyRiskRatio}: The efficacy boundaries on the risk
#'       ratio scale.
#'
#'     - \code{futilityRiskRatio}: The futility boundaries on the risk
#'       ratio scale.
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
#'     - \code{varianceRatio}: The ratio of the variance under H0 to
#'       the variance under H1.
#'
#'     - \code{calculationTarget}: The calculation target, \code{beta} or
#'       \code{n}.
#'
#'     - \code{nullVariance}: Whether to use the variance under the null or
#'       the empirical variance under the alternative.
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
#' (design1 <- getDesignRiskRatio(
#'   beta = 0.1, n = NA, pi1 = 0.5, pi2 = 0.3,
#'   alpha = 0.05))
#'
#' @export
getDesignRiskRatio <- function(
    beta = NA_real_,
    n = NA_real_,
    riskRatioH0 = 1,
    pi1 = NA_real_,
    pi2 = NA_real_,
    nullVariance = TRUE,
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
    stop("beta must lie in [0.0001, 1-alpha)");
  }

  if (!is.na(n) && n <= 0) {
    stop("n must be positive")
  }

  if (riskRatioH0 <= 0) {
    stop("riskRatioH0 must be positive")
  }

  if (pi1 <= 0 || pi1 >= 1) {
    stop("pi1 must lie between 0 and 1")
  }

  if (pi2 <= 0 || pi2 >= 1) {
    stop("pi2 must lie between 0 and 1")
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

  riskRatio = pi1/pi2

  directionUpper = riskRatio > riskRatioH0

  theta = ifelse(directionUpper, log(riskRatio) - log(riskRatioH0),
                 log(riskRatioH0) - log(riskRatio))

  # restricted maximum likelihood estimates
  mr = remlRiskRatio(riskRatioH0, r, r*pi1, 1-r, (1-r)*pi2)
  p1 = mr[1]
  p2 = mr[2]

  # variance for one sampling unit
  v0 = (1-p1)/(r*p1) + (1-p2)/((1-r)*p2)
  v1 = (1-pi1)/(r*pi1) + (1-pi2)/((1-r)*pi2)
  varianceRatio = ifelse(nullVariance, v0/v1, 1)

  if (!is.na(n)) { # power calculation
    if (rounding) {
      n = ceiling(n)
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
      spendingTime, varianceRatio)
  } else { # sample size calculation
    des = getDesign(
      beta, IMax = NA, theta,
      kMax, informationRates,
      efficacyStopping, futilityStopping,
      criticalValues, alpha, typeAlphaSpending,
      parameterAlphaSpending, userAlphaSpending,
      futilityBounds, typeBetaSpending,
      parameterBetaSpending, userBetaSpending,
      spendingTime, varianceRatio)

    n = des$overallResults$information*v1

    if (rounding) {
      n = ceiling(n)
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
        spendingTime, varianceRatio)
    }
  }

  des$overallResults$theta = theta
  des$overallResults$numberOfSubjects = n
  des$overallResults$expectedNumberOfSubjectsH1 =
    des$overallResults$expectedInformationH1*v1
  des$overallResults$expectedNumberOfSubjectsH0 =
    des$overallResults$expectedInformationH0*v1
  des$overallResults$riskRatioH0 = riskRatioH0
  des$overallResults$pi1 = pi1
  des$overallResults$pi2 = pi2

  if (directionUpper) {
    des$byStageResults$efficacyRiskRatio =
      exp(des$byStageResults$efficacyTheta)*riskRatioH0
    des$byStageResults$futilityRiskRatio =
      exp(des$byStageResults$futilityTheta)*riskRatioH0
  } else {
    des$byStageResults$efficacyRiskRatio =
      exp(-des$byStageResults$efficacyTheta)*riskRatioH0
    des$byStageResults$futilityRiskRatio =
      exp(-des$byStageResults$futilityTheta)*riskRatioH0
  }

  des$byStageResults$efficacyTheta = NULL
  des$byStageResults$futilityTheta = NULL
  des$byStageResults$numberOfSubjects =
    des$byStageResults$informationRates*n

  des$settings$calculationTarget = ifelse(is.na(beta), "beta", "n")
  des$settings$nullVariance = nullVariance
  des$settings$allocationRatioPlanned = allocationRatioPlanned
  des$settings$rounding = rounding

  attr(des, "class") = "designRiskRatio"

  des
}


#' @title Group sequential design for two-sample risk ratio based on
#' the Farrington-Manning score test
#' @description Obtains the power given sample size or obtains the sample
#' size given power for a group sequential design for two-sample risk
#' ratio based on the Farrington-Manning score test
#'
#' @param beta The type II error.
#' @param n The total sample size.
#' @param riskRatioH0 The risk ratio under the null hypothesis.
#'   Defaults to 1.
#' @param pi1 The assumed probability for the active treatment group.
#' @param pi2 The assumed probability for the control group.
#' @param nullVariance Whether to use the variance under the null or
#'   the empirical variance under the alternative.
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
#' @return An S3 class \code{designRiskRatioFM} object with three
#' components:
#'
#' * \code{overallResults}: A data frame containing the following variables:
#'
#'    - \code{overallReject}: The overall rejection probability.
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
#'     - \code{riskRatioH0}: The risk ratio under the null hypothesis.
#'
#'     - \code{pi1}: The assumed probability for the active treatment group.
#'
#'     - \code{pi2}: The assumed probability for the control group.
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
#'     - \code{efficacyRiskRatioScore}: The efficacy boundaries on the
#'       score test \code{pi1 - riskRatioH0*pi2} score.
#'
#'     - \code{futilityRiskRatioScore}: The futility boundaries on the
#'       score test \code{pi1 - riskRatioH0*pi2} scale.
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
#'     - \code{varianceRatio}: The ratio of the variance under H0 to
#'       the variance under H1.
#'
#'     - \code{calculationTarget}: The calculation target, \code{beta} or
#'       \code{n}.
#'
#'     - \code{nullVariance}: Whether to use the variance under the null or
#'       the empirical variance under the alternative.
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
#' (design1 <- getDesignRiskRatioFM(
#'   beta = 0.2, riskRatioH0 = 1.3, pi1 = 0.125, pi2 = 0.125,
#'   alpha = 0.05))
#'
#' @export
getDesignRiskRatioFM <- function(
    beta = NA_real_,
    n = NA_real_,
    riskRatioH0 = 1,
    pi1 = NA_real_,
    pi2 = NA_real_,
    nullVariance = TRUE,
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
    stop("beta must lie in [0.0001, 1-alpha)");
  }

  if (!is.na(n) && n <= 0) {
    stop("n must be positive")
  }

  if (riskRatioH0 <= 0) {
    stop("riskRatioH0 must be positive")
  }

  if (pi1 <= 0 || pi1 >= 1) {
    stop("pi1 must lie between 0 and 1")
  }

  if (pi2 <= 0 || pi2 >= 1) {
    stop("pi2 must lie between 0 and 1")
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

  riskRatio = pi1/pi2

  directionUpper = riskRatio > riskRatioH0

  theta = ifelse(directionUpper, pi1 - riskRatioH0*pi2,
                 -pi1 + riskRatioH0*pi2)

  # restricted maximum likelihood estimates
  mr = remlRiskRatio(riskRatioH0, r, r*pi1, 1-r, (1-r)*pi2)
  p1 = mr[1]
  p2 = mr[2]

  # variance for one sampling unit
  v0 = p1*(1-p1)/r + riskRatioH0^2*p2*(1-p2)/(1-r)
  v1 = pi1*(1-pi1)/r + riskRatioH0^2*pi2*(1-pi2)/(1-r)

  varianceRatio = ifelse(nullVariance, v0/v1, 1)

  if (!is.na(n)) { # power calculation
    if (rounding) {
      n = ceiling(n)
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
      spendingTime, varianceRatio)
  } else { # sample size calculation
    des = getDesign(
      beta, IMax = NA, theta,
      kMax, informationRates,
      efficacyStopping, futilityStopping,
      criticalValues, alpha, typeAlphaSpending,
      parameterAlphaSpending, userAlphaSpending,
      futilityBounds, typeBetaSpending,
      parameterBetaSpending, userBetaSpending,
      spendingTime, varianceRatio)

    n = des$overallResults$information*v1

    if (rounding) {
      n = ceiling(n)
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
        spendingTime, varianceRatio)
    }
  }

  des$overallResults$theta = theta
  des$overallResults$numberOfSubjects = n
  des$overallResults$expectedNumberOfSubjectsH1 =
    des$overallResults$expectedInformationH1*v1
  des$overallResults$expectedNumberOfSubjectsH0 =
    des$overallResults$expectedInformationH0*v1
  des$overallResults$riskRatioH0 = riskRatioH0
  des$overallResults$pi1 = pi1
  des$overallResults$pi2 = pi2

  if (directionUpper) {
    des$byStageResults$efficacyRiskRatioScore =
      des$byStageResults$efficacyTheta
    des$byStageResults$futilityRiskRatioScore =
      des$byStageResults$futilityTheta
  } else {
    des$byStageResults$efficacyRiskRatioScore =
      -des$byStageResults$efficacyTheta
    des$byStageResults$futilityRiskRatioScore =
      -des$byStageResults$futilityTheta
  }

  des$byStageResults$efficacyTheta = NULL
  des$byStageResults$futilityTheta = NULL
  des$byStageResults$numberOfSubjects =
    des$byStageResults$informationRates*n

  des$settings$calculationTarget = ifelse(is.na(beta), "beta", "n")
  des$settings$nullVariance = nullVariance
  des$settings$allocationRatioPlanned = allocationRatioPlanned
  des$settings$rounding = rounding

  attr(des, "class") = "designRiskRatioFM"

  des
}


#' @title Group sequential design for two-sample odds ratio
#' @description Obtains the power given sample size or obtains the sample
#' size given power for a group sequential design for two-sample odds
#' ratio.
#'
#' @param beta The type II error.
#' @param n The total sample size.
#' @param oddsRatioH0 The odds ratio under the null hypothesis.
#'   Defaults to 1.
#' @param pi1 The assumed probability for the active treatment group.
#' @param pi2 The assumed probability for the control group.
#' @param nullVariance Whether to use the variance under the null
#'   or the empirical variance under the alternative.
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
#' @return An S3 class \code{designOddsRatio} object with three components:
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
#'     - \code{oddsRatioH0}: The odds ratio under the null hypothesis.
#'
#'     - \code{pi1}: The assumed probability for the active treatment group.
#'
#'     - \code{pi2}: The assumed probability for the control group.
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
#'     - \code{efficacyOdddsRatio}: The efficacy boundaries on the odds
#'       ratio scale.
#'
#'     - \code{futilityOddsRatio}: The futility boundaries on the odds
#'       ratio scale.
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
#'     - \code{varianceRatio}: The ratio of the variance under H0 to
#'       the variance under H1.
#'
#'     - \code{calculationTarget}: The calculation target, \code{beta} or
#'       \code{n}.
#'
#'     - \code{nullVariance}: Whether to use the variance under the null
#'       or the empirical variance under the alternative.
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
#' (design1 <- getDesignOddsRatio(
#'   beta = 0.1, n = NA, pi1 = 0.5, pi2 = 0.3,
#'   alpha = 0.05))
#'
#' @export
getDesignOddsRatio <- function(
    beta = NA_real_,
    n = NA_real_,
    oddsRatioH0 = 1,
    pi1 = NA_real_,
    pi2 = NA_real_,
    nullVariance = FALSE,
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
    stop("beta must lie in [0.0001, 1-alpha)");
  }

  if (!is.na(n) && n <= 0) {
    stop("n must be positive")
  }

  if (oddsRatioH0 <= 0) {
    stop("oddsRatioH0 must be positive")
  }

  if (pi1 <= 0 || pi1 >= 1) {
    stop("pi1 must lie between 0 and 1")
  }

  if (pi2 <= 0 || pi2 >= 1) {
    stop("pi2 must lie between 0 and 1")
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

  oddsRatio = pi1*(1-pi2)/((1-pi1)*pi2)

  directionUpper = oddsRatio > oddsRatioH0

  theta = ifelse(directionUpper, log(oddsRatio) - log(oddsRatioH0),
                 log(oddsRatioH0) - log(oddsRatio))

  # restricted maximum likelihood estimates
  mr = remlOddsRatio(oddsRatioH0, r, r*pi1, 1-r, (1-r)*pi2)
  p1 = mr[1]
  p2 = mr[2]

  # variance for one sampling unit
  v0 = 1/(r*p1*(1-p1)) + 1/((1-r)*p2*(1-p2))
  v1 = 1/(r*pi1*(1-pi1)) + 1/((1-r)*pi2*(1-pi2))

  varianceRatio = ifelse(nullVariance, v0/v1, 1)

  if (!is.na(n)) { # power calculation
    if (rounding) {
      n = ceiling(n)
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
      spendingTime, varianceRatio)
  } else { # sample size calculation
    des = getDesign(
      beta, IMax = NA, theta,
      kMax, informationRates,
      efficacyStopping, futilityStopping,
      criticalValues, alpha, typeAlphaSpending,
      parameterAlphaSpending, userAlphaSpending,
      futilityBounds, typeBetaSpending,
      parameterBetaSpending, userBetaSpending,
      spendingTime, varianceRatio)

    n = des$overallResults$information*v1

    if (rounding) {
      n = ceiling(n)
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
        spendingTime, varianceRatio)
    }
  }

  des$overallResults$theta = theta
  des$overallResults$numberOfSubjects = n
  des$overallResults$expectedNumberOfSubjectsH1 =
    des$overallResults$expectedInformationH1*v1
  des$overallResults$expectedNumberOfSubjectsH0 =
    des$overallResults$expectedInformationH0*v1
  des$overallResults$oddsRatioH0 = oddsRatioH0
  des$overallResults$pi1 = pi1
  des$overallResults$pi2 = pi2

  if (directionUpper) {
    des$byStageResults$efficacyOddsRatio =
      exp(des$byStageResults$efficacyTheta)*oddsRatioH0
    des$byStageResults$futilityOddsRatio =
      exp(des$byStageResults$futilityTheta)*oddsRatioH0
  } else {
    des$byStageResults$efficacyOddsRatio =
      exp(-des$byStageResults$efficacyTheta)*oddsRatioH0
    des$byStageResults$futilityOddsRatio =
      exp(-des$byStageResults$futilityTheta)*oddsRatioH0
  }

  des$byStageResults$efficacyTheta = NULL
  des$byStageResults$futilityTheta = NULL
  des$byStageResults$numberOfSubjects =
    des$byStageResults$informationRates*n

  des$settings$calculationTarget = ifelse(is.na(beta), "beta", "n")
  des$settings$nullVariance = nullVariance
  des$settings$allocationRatioPlanned = allocationRatioPlanned
  des$settings$rounding = rounding

  attr(des, "class") = "designOddsRatio"

  des
}


#' @title Group sequential design for equivalence in two-sample
#' risk difference
#' @description Obtains the power given sample size or obtains the sample
#' size given power for a group sequential design for equivalence in
#' two-sample risk difference.
#'
#' @param beta The type II error.
#' @param n The total sample size.
#' @param riskDiffLower The lower equivalence limit of risk difference.
#' @param riskDiffUpper The upper equivalence limit of risk difference.
#' @param pi1 The assumed probability for the active treatment group.
#' @param pi2 The assumed probability for the control group.
#' @param nullVariance Whether to use the variance under the null or
#'   the empirical variance under the alternative.
#' @param allocationRatioPlanned Allocation ratio for the active treatment
#'   versus control. Defaults to 1 for equal randomization.
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
#' @return An S3 class \code{designRiskDiffEquiv} object with three
#' components:
#'
#' * \code{overallResults}: A data frame containing the following variables:
#'
#'     - \code{overallReject}: The overall rejection probability.
#'
#'     - \code{alpha}: The significance level for each of the two one-sided
#'       tests. Defaults to 0.05.
#'
#'     - \code{attainedAlphaH10}: The attained significance level under H10.
#'
#'     - \code{attainedAlphaH20}: The attained significance level under H20.
#'
#'     - \code{kMax}: The number of stages.
#'
#'     - \code{information}: The maximum information.
#'
#'     - \code{expectedInformationH1}: The expected information under H1.
#'
#'     - \code{expectedInformationH10}: The expected information under H10.
#'
#'     - \code{expectedInformationH20}: The expected information under H20.
#'
#'     - \code{numberOfSubjects}: The maximum number of subjects.
#'
#'     - \code{expectedNumberOfSubjectsH1}: The expected number of subjects
#'       under H1.
#'
#'     - \code{expectedNumberOfSubjectsH10}: The expected number of subjects
#'       under H10.
#'
#'     - \code{expectedNumberOfSubjectsH20}: The expected number of subjects
#'       under H20.
#'
#'     - \code{riskDiffLower}: The lower equivalence limit of risk
#'       difference.
#'
#'     - \code{riskDiffUpper}: The upper equivalence limit of risk
#'       difference.
#'
#'     - \code{pi1}: The assumed probability for the active treatment group.
#'
#'     - \code{pi2}: The assumed probability for the control group.
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
#'     - \code{cumulativeAttainedAlphaH10}: The cumulative alpha attained
#'       under H10.
#'
#'     - \code{cumulativeAttainedAlphaH20}: The cumulative alpha attained
#'       under H20.
#'
#'     - \code{efficacyP}: The efficacy bounds on the p-value scale for
#'       each of the two one-sided tests.
#'
#'     - \code{information}: The cumulative information.
#'
#'     - \code{efficacyRiskDiffLower}: The efficacy boundaries on the
#'       risk difference scale for the one-sided null hypothesis on the
#'       lower equivalence limit.
#'
#'     - \code{efficacyRiskDiffUpper}: The efficacy boundaries on the
#'       risk difference scale for the one-sided null hypothesis on the
#'       upper equivalence limit.
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
#'     - \code{calculationTarget}: The calculation target, \code{beta} or
#'       \code{n}.
#'
#'     - \code{nullVariance}: Whether to use the variance under the null or
#'       the empirical variance under the alternative.
#'
#'     - \code{varianceRatioH10}: The ratio of the variance under H10 to
#'       the variance under H1.
#'
#'     - \code{varianceRatioH20}: The ratio of the variance under H20 to
#'       the variance under H1.
#'
#'     - \code{varianceRatioH12}: The ratio of the variance under H10 to
#'       the variance under H20.
#'
#'     - \code{varianceRatioH21}: The ratio of the variance under H20 to
#'       the variance under H10.
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
#' (design1 <- getDesignRiskDiffEquiv(
#'   beta = 0.2, n = NA, riskDiffLower = -0.1,
#'   riskDiffUpper = 0.1, pi1 = 0.12, pi2 = 0.12,
#'   nullVariance = 1,
#'   kMax = 3, alpha = 0.05, typeAlphaSpending = "sfOF"))
#'
#' @export
getDesignRiskDiffEquiv <- function(
    beta = NA_real_,
    n = NA_real_,
    riskDiffLower = NA_real_,
    riskDiffUpper = NA_real_,
    pi1 = NA_real_,
    pi2 = NA_real_,
    nullVariance = FALSE,
    allocationRatioPlanned = 1,
    rounding = TRUE,
    kMax = 1L,
    informationRates = NA_real_,
    criticalValues = NA_real_,
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
    stop("beta must lie in [0.0001, 1-alpha)");
  }

  if (!is.na(n) && n <= 0) {
    stop("n must be positive")
  }

  if (is.na(riskDiffLower)) {
    stop("riskDiffLower must be provided")
  }

  if (is.na(riskDiffUpper)) {
    stop("riskDiffUpper must be provided")
  }

  if (riskDiffLower <= -1) {
    stop("riskDiffLower must be greater than -1")
  }

  if (riskDiffUpper >= 1) {
    stop("riskDiffUpper must be less than 1")
  }

  if (is.na(pi1)) {
    stop("pi1 must be provided")
  }

  if (is.na(pi2)) {
    stop("pi2 must be provided")
  }

  if (pi1 <= 0 || pi1 >= 1) {
    stop("pi1 must lie between 0 and 1")
  }

  if (pi2 <= 0 || pi2 >= 1) {
    stop("pi2 must lie between 0 and 1")
  }

  riskDiff = pi1 - pi2

  if (riskDiffLower >= riskDiff) {
    stop("riskDiffLower must be less than pi1 - pi2")
  }

  if (riskDiffUpper <= riskDiff) {
    stop("riskDiffUpper must be greater than pi1 - pi2")
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

  # variance for one sampling unit under H1
  vH1 = pi1*(1-pi1)/r + pi2*(1-pi2)/(1-r)

  if (nullVariance) {
    # variance for one sampling unit under H10
    pi1H10 = pi2 + riskDiffLower
    vH10 = pi1H10*(1-pi1H10)/r + pi2*(1-pi2)/(1-r)

    # variance for one sampling unit under H20
    pi1H20 = pi2 + riskDiffUpper
    vH20 = pi1H20*(1-pi1H20)/r + pi2*(1-pi2)/(1-r)

    # restricted maximum likelihood estimates under H10 w.r.t. H1
    mr = remlRiskDiff(riskDiffLower, r, r*pi1, 1-r, (1-r)*pi2)
    vH10H1 = mr[1]*(1-mr[1])/r + mr[2]*(1-mr[2])/(1-r)

    # restricted maximum likelihood estimates under H20 w.r.t. H1
    mr = remlRiskDiff(riskDiffUpper, r, r*pi1, 1-r, (1-r)*pi2)
    vH20H1 = mr[1]*(1-mr[1])/r + mr[2]*(1-mr[2])/(1-r)

    # restricted maximum likelihood estimates under H10 w.r.t. H20
    mr = remlRiskDiff(riskDiffLower, r, r*pi1H20, 1-r, (1-r)*pi2)
    vH10H20 = mr[1]*(1-mr[1])/r + mr[2]*(1-mr[2])/(1-r)

    # restricted maximum likelihood estimates under H20 w.r.t. H10
    mr = remlRiskDiff(riskDiffUpper, r, r*pi1H10, 1-r, (1-r)*pi2)
    vH20H10 = mr[1]*(1-mr[1])/r + mr[2]*(1-mr[2])/(1-r)

    # calculate the variance ratios
    varianceRatioH10 = vH10H1/vH1
    varianceRatioH20 = vH20H1/vH1
    varianceRatioH12 = vH10H20/vH20
    varianceRatioH21 = vH20H10/vH10
  } else {
    varianceRatioH10 = varianceRatioH20 = 1
    varianceRatioH12 = varianceRatioH21 = 1
  }


  # sample size calculation
  if (is.na(n)) {
    des = getDesignEquiv(
      beta = beta, IMax = NA, thetaLower = riskDiffLower,
      thetaUpper = riskDiffUpper, theta = riskDiff,
      kMax = kMax, informationRates = informationRates,
      alpha = alpha, typeAlphaSpending = typeAlphaSpending,
      parameterAlphaSpending = parameterAlphaSpending,
      userAlphaSpending = userAlphaSpending,
      spendingTime = spendingTime,
      varianceRatioH10 = varianceRatioH10,
      varianceRatioH20 = varianceRatioH20,
      varianceRatioH12 = varianceRatioH12,
      varianceRatioH21 = varianceRatioH21)

    n = des$overallResults$information*vH1
  }

  if (rounding) {
    n = ceiling(n)
    informationRates = round(n*informationRates)/n
  }

  if (is.na(beta) || rounding) { # calculate power
    des = getDesignEquiv(
      beta = NA, IMax = n/vH1, thetaLower = riskDiffLower,
      thetaUpper = riskDiffUpper, theta = riskDiff,
      kMax = kMax, informationRates = informationRates,
      alpha = alpha, typeAlphaSpending = typeAlphaSpending,
      parameterAlphaSpending = parameterAlphaSpending,
      userAlphaSpending = userAlphaSpending,
      spendingTime = spendingTime,
      varianceRatioH10 = varianceRatioH10,
      varianceRatioH20 = varianceRatioH20,
      varianceRatioH12 = varianceRatioH12,
      varianceRatioH21 = varianceRatioH21)
  }

  des$overallResults$numberOfSubjects = n
  des$overallResults$expectedNumberOfSubjectsH1 =
    des$overallResults$expectedInformationH1*vH1
  des$overallResults$expectedNumberOfSubjectsH10 =
    des$overallResults$expectedInformationH10*vH1
  des$overallResults$expectedNumberOfSubjectsH20 =
    des$overallResults$expectedInformationH20*vH1
  des$overallResults$riskDiffLower = riskDiffLower
  des$overallResults$riskDiffUpper = riskDiffUpper
  des$overallResults$pi1 = pi1
  des$overallResults$pi2 = pi2
  des$overallResults <-
    des$overallResults[, c("overallReject", "alpha",
                           "attainedAlphaH10", "attainedAlphaH20",
                           "kMax", "information",
                           "expectedInformationH1",
                           "expectedInformationH10",
                           "expectedInformationH20",
                           "numberOfSubjects",
                           "expectedNumberOfSubjectsH1",
                           "expectedNumberOfSubjectsH10",
                           "expectedNumberOfSubjectsH20",
                           "riskDiffLower", "riskDiffUpper",
                           "pi1", "pi2")]

  des$byStageResults$efficacyRiskDiffLower =
    des$byStageResults$efficacyThetaLower
  des$byStageResults$efficacyRiskDiffUpper =
    des$byStageResults$efficacyThetaUpper
  des$byStageResults$numberOfSubjects = n*informationRates
  des$byStageResults <-
    des$byStageResults[, c("informationRates", "efficacyBounds",
                           "rejectPerStage", "cumulativeRejection",
                           "cumulativeAlphaSpent",
                           "cumulativeAttainedAlphaH10",
                           "cumulativeAttainedAlphaH20",
                           "efficacyRiskDiffLower",
                           "efficacyRiskDiffUpper", "efficacyP",
                           "information", "numberOfSubjects")]

  des$settings$nullVariance = nullVariance
  des$settings$allocationRatioPlanned = allocationRatioPlanned
  des$settings$rounding = rounding
  des$settings$calculationTarget = ifelse(is.na(beta), "beta", "n")
  des$settings <-
    des$settings[c("typeAlphaSpending", "parameterAlphaSpending",
                   "userAlphaSpending", "spendingTime",
                   "calculationTarget", "nullVariance",
                   "varianceRatioH10", "varianceRatioH20",
                   "varianceRatioH12", "varianceRatioH21",
                   "allocationRatioPlanned", "rounding")]

  attr(des, "class") = "designRiskDiffEquiv"

  des
}


#' @title Group sequential design for equivalence in two-sample
#' risk ratio
#' @description Obtains the power given sample size or obtains the sample
#' size given power for a group sequential design for equivalence in
#' two-sample risk ratio.
#'
#' @param beta The type II error.
#' @param n The total sample size.
#' @param riskRatioLower The lower equivalence limit of risk ratio.
#' @param riskRatioUpper The upper equivalence limit of risk ratio.
#' @param pi1 The assumed probability for the active treatment group.
#' @param pi2 The assumed probability for the control group.
#' @param nullVariance Whether to use the variance under the null or
#'   the empirical variance under the alternative.
#' @param allocationRatioPlanned Allocation ratio for the active treatment
#'   versus control. Defaults to 1 for equal randomization.
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
#' @return An S3 class \code{designRiskRatioEquiv} object with three
#' components:
#'
#' * \code{overallResults}: A data frame containing the following variables:
#'
#'     - \code{overallReject}: The overall rejection probability.
#'
#'     - \code{alpha}: The significance level for each of the two one-sided
#'       tests. Defaults to 0.05.
#'
#'     - \code{attainedAlphaH10}: The attained significance level under H10.
#'
#'     - \code{attainedAlphaH20}: The attained significance level under H20.
#'
#'     - \code{kMax}: The number of stages.
#'
#'     - \code{information}: The maximum information.
#'
#'     - \code{expectedInformationH1}: The expected information under H1.
#'
#'     - \code{expectedInformationH10}: The expected information under H10.
#'
#'     - \code{expectedInformationH20}: The expected information under H20.
#'
#'     - \code{numberOfSubjects}: The maximum number of subjects.
#'
#'     - \code{expectedNumberOfSubjectsH1}: The expected number of subjects
#'       under H1.
#'
#'     - \code{expectedNumberOfSubjectsH10}: The expected number of subjects
#'       under H10.
#'
#'     - \code{expectedNumberOfSubjectsH20}: The expected number of subjects
#'       under H20.
#'
#'     - \code{riskRatioLower}: The lower equivalence limit of risk ratio.
#'
#'     - \code{riskRatioUpper}: The upper equivalence limit of risk ratio.
#'
#'     - \code{pi1}: The assumed probability for the active treatment group.
#'
#'     - \code{pi2}: The assumed probability for the control group.
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
#'     - \code{cumulativeAttainedAlphaH10}: The cumulative alpha attained
#'       under H10.
#'
#'     - \code{cumulativeAttainedAlphaH20}: The cumulative alpha attained
#'       under H20.
#'
#'     - \code{efficacyRiskRatioLower}: The efficacy boundaries on the
#'       risk ratio scale for the one-sided null hypothesis on the
#'       lower equivalence limit.
#'
#'     - \code{efficacyRiskRatioUpper}: The efficacy boundaries on the
#'       risk ratio scale for the one-sided null hypothesis on the
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
#'     - \code{calculationTarget}: The calculation target, \code{beta} or
#'       \code{n}.
#'
#'     - \code{nullVariance}: Whether to use the variance under the null or
#'       the empirical variance under the alternative.
#'
#'     - \code{varianceRatioH10}: The ratio of the variance under H10 to
#'       the variance under H1.
#'
#'     - \code{varianceRatioH20}: The ratio of the variance under H20 to
#'       the variance under H1.
#'
#'     - \code{varianceRatioH12}: The ratio of the variance under H10 to
#'       the variance under H20.
#'
#'     - \code{varianceRatioH21}: The ratio of the variance under H20 to
#'       the variance under H10.
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
#' (design1 <- getDesignRiskRatioEquiv(
#'   beta = 0.2, n = NA, riskRatioLower = 0.8,
#'   riskRatioUpper = 1.25, pi1 = 0.12, pi2 = 0.12,
#'   kMax = 3, alpha = 0.05, typeAlphaSpending = "sfOF"))
#'
#' @export
getDesignRiskRatioEquiv <- function(
    beta = NA_real_,
    n = NA_real_,
    riskRatioLower = NA_real_,
    riskRatioUpper = NA_real_,
    pi1 = NA_real_,
    pi2 = NA_real_,
    nullVariance = FALSE,
    allocationRatioPlanned = 1,
    rounding = TRUE,
    kMax = 1L,
    informationRates = NA_real_,
    criticalValues = NA_real_,
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
    stop("beta must lie in [0.0001, 1-alpha)");
  }

  if (!is.na(n) && n <= 0) {
    stop("n must be positive")
  }

  if (is.na(riskRatioLower)) {
    stop("riskRatioLower must be provided")
  }

  if (is.na(riskRatioUpper)) {
    stop("riskRatioUpper must be provided")
  }

  if (riskRatioLower <= 0) {
    stop("riskRatioLower must be positive")
  }

  if (is.na(pi1)) {
    stop("pi1 must be provided")
  }

  if (is.na(pi2)) {
    stop("pi2 must be provided")
  }

  if (pi1 <= 0 || pi1 >= 1) {
    stop("pi1 must lie between 0 and 1")
  }

  if (pi2 <= 0 || pi2 >= 1) {
    stop("pi2 must lie between 0 and 1")
  }

  riskRatio = pi1/pi2

  if (riskRatioLower >= riskRatio) {
    stop("riskRatioLower must be less than pi1/pi2")
  }

  if (riskRatioUpper <= riskRatio) {
    stop("riskRatioUpper must be greater than pi1/pi2")
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

  # variance for one sampling unit
  vH1 = (1-pi1)/(r*pi1) + (1-pi2)/((1-r)*pi2)

  if (nullVariance) {
    # variance for one sampling unit under H10
    pi1H10 = pi2*riskRatioLower
    vH10 = (1-pi1H10)/(r*pi1H10) + (1-pi2)/((1-r)*pi2)

    # variance for one sampling unit under H20
    pi1H20 = pi2*riskRatioUpper
    vH20 = (1-pi1H20)/(r*pi1H20) + (1-pi2)/((1-r)*pi2)

    # restricted maximum likelihood estimates under H10 w.r.t. H1
    mr = remlRiskRatio(riskRatioLower, r, r*pi1, 1-r, (1-r)*pi2)
    vH10H1 = (1-mr[1])/(r*mr[1]) + (1-mr[2])/((1-r)*mr[2])

    # restricted maximum likelihood estimates under H20 w.r.t. H1
    mr = remlRiskRatio(riskRatioUpper, r, r*pi1, 1-r, (1-r)*pi2)
    vH20H1 = (1-mr[1])/(r*mr[1]) + (1-mr[2])/((1-r)*mr[2])

    # restricted maximum likelihood estimates under H10 w.r.t. H20
    mr = remlRiskRatio(riskRatioLower, r, r*pi1H20, 1-r, (1-r)*pi2)
    vH10H20 = (1-mr[1])/(r*mr[1]) + (1-mr[2])/((1-r)*mr[2])

    # restricted maximum likelihood estimates under H20 w.r.t. H10
    mr = remlRiskRatio(riskRatioUpper, r, r*pi1H10, 1-r, (1-r)*pi2)
    vH20H10 = (1-mr[1])/(r*mr[1]) + (1-mr[2])/((1-r)*mr[2])

    # calculate the variance ratios
    varianceRatioH10 = vH10H1/vH1
    varianceRatioH20 = vH20H1/vH1
    varianceRatioH12 = vH10H20/vH20
    varianceRatioH21 = vH20H10/vH10
  } else {
    varianceRatioH10 = varianceRatioH20 = 1
    varianceRatioH12 = varianceRatioH21 = 1
  }

  # sample size calculation
  if (is.na(n)) {
    des = getDesignEquiv(
      beta = beta, IMax = NA, thetaLower = log(riskRatioLower),
      thetaUpper = log(riskRatioUpper), theta = log(riskRatio),
      kMax = kMax, informationRates = informationRates,
      alpha = alpha, typeAlphaSpending = typeAlphaSpending,
      parameterAlphaSpending = parameterAlphaSpending,
      userAlphaSpending = userAlphaSpending,
      spendingTime = spendingTime,
      varianceRatioH10 = varianceRatioH10,
      varianceRatioH20 = varianceRatioH20,
      varianceRatioH12 = varianceRatioH12,
      varianceRatioH21 = varianceRatioH21)

    n = des$overallResults$information*vH1
  }

  if (rounding) {
    n = ceiling(n)
    informationRates = round(n*informationRates)/n
  }

  if (is.na(beta) || rounding) { # calculate power
    des = getDesignEquiv(
      beta = NA, IMax = n/vH1, thetaLower = log(riskRatioLower),
      thetaUpper = log(riskRatioUpper), theta = log(riskRatio),
      kMax = kMax, informationRates = informationRates,
      alpha = alpha, typeAlphaSpending = typeAlphaSpending,
      parameterAlphaSpending = parameterAlphaSpending,
      userAlphaSpending = userAlphaSpending,
      spendingTime = spendingTime,
      varianceRatioH10 = varianceRatioH10,
      varianceRatioH20 = varianceRatioH20,
      varianceRatioH12 = varianceRatioH12,
      varianceRatioH21 = varianceRatioH21)
  }

  des$overallResults$numberOfSubjects = n
  des$overallResults$expectedNumberOfSubjectsH1 =
    des$overallResults$expectedInformationH1*vH1
  des$overallResults$expectedNumberOfSubjectsH10 =
    des$overallResults$expectedInformationH10*vH1
  des$overallResults$expectedNumberOfSubjectsH20 =
    des$overallResults$expectedInformationH20*vH1
  des$overallResults$riskRatioLower = riskRatioLower
  des$overallResults$riskRatioUpper = riskRatioUpper
  des$overallResults$pi1 = pi1
  des$overallResults$pi2 = pi2
  des$overallResults <-
    des$overallResults[, c("overallReject", "alpha",
                           "attainedAlphaH10", "attainedAlphaH20",
                           "kMax", "information",
                           "expectedInformationH1",
                           "expectedInformationH10",
                           "expectedInformationH20",
                           "numberOfSubjects",
                           "expectedNumberOfSubjectsH1",
                           "expectedNumberOfSubjectsH10",
                           "expectedNumberOfSubjectsH20",
                           "riskRatioLower", "riskRatioUpper",
                           "pi1", "pi2")]

  des$byStageResults$efficacyRiskRatioLower =
    exp(des$byStageResults$efficacyThetaLower)
  des$byStageResults$efficacyRiskRatioUpper =
    exp(des$byStageResults$efficacyThetaUpper)
  des$byStageResults$numberOfSubjects = n*informationRates
  des$byStageResults <-
    des$byStageResults[, c("informationRates", "efficacyBounds",
                           "rejectPerStage", "cumulativeRejection",
                           "cumulativeAlphaSpent",
                           "cumulativeAttainedAlphaH10",
                           "cumulativeAttainedAlphaH20",
                           "efficacyRiskRatioLower",
                           "efficacyRiskRatioUpper", "efficacyP",
                           "information", "numberOfSubjects")]

  des$settings$nullVariance = nullVariance
  des$settings$allocationRatioPlanned = allocationRatioPlanned
  des$settings$rounding = rounding
  des$settings$calculationTarget = ifelse(is.na(beta), "beta", "n")
  des$settings <-
    des$settings[c("typeAlphaSpending", "parameterAlphaSpending",
                   "userAlphaSpending", "spendingTime",
                   "calculationTarget", "nullVariance",
                   "varianceRatioH10", "varianceRatioH20",
                   "varianceRatioH12", "varianceRatioH21",
                   "allocationRatioPlanned", "rounding")]

  attr(des, "class") = "designRiskRatioEquiv"

  des
}


#' @title Group sequential design for equivalence in two-sample
#' odds ratio
#' @description Obtains the power given sample size or obtains the sample
#' size given power for a group sequential design for equivalence in
#' two-sample odds ratio.
#'
#' @param beta The type II error.
#' @param n The total sample size.
#' @param oddsRatioLower The lower equivalence limit of odds ratio.
#' @param oddsRatioUpper The upper equivalence limit of odds ratio.
#' @param pi1 The assumed probability for the active treatment group.
#' @param pi2 The assumed probability for the control group.
#' @param nullVariance Whether to use the variance under the null
#'   or the empirical variance under the alternative.
#' @param allocationRatioPlanned Allocation ratio for the active treatment
#'   versus control. Defaults to 1 for equal randomization.
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
#' @return An S3 class \code{designOddsRatioEquiv} object with three
#' components:
#'
#' * \code{overallResults}: A data frame containing the following variables:
#'
#'     - \code{overallReject}: The overall rejection probability.
#'
#'     - \code{alpha}: The significance level for each of the two one-sided
#'       tests. Defaults to 0.05.
#'
#'     - \code{attainedAlphaH10}: The attained significance level under H10.
#'
#'     - \code{attainedAlphaH20}: The attained significance level under H20.
#'
#'     - \code{kMax}: The number of stages.
#'
#'     - \code{information}: The maximum information.
#'
#'     - \code{expectedInformationH1}: The expected information under H1.
#'
#'     - \code{expectedInformationH10}: The expected information under H10.
#'
#'     - \code{expectedInformationH20}: The expected information under H20.
#'
#'     - \code{numberOfSubjects}: The maximum number of subjects.
#'
#'     - \code{expectedNumberOfSubjectsH1}: The expected number of subjects
#'       under H1.
#'
#'     - \code{expectedNumberOfSubjectsH10}: The expected number of subjects
#'       under H10.
#'
#'     - \code{expectedNumberOfSubjectsH20}: The expected number of subjects
#'       under H20.
#'
#'     - \code{oddsRatioLower}: The lower equivalence limit of odds ratio.
#'
#'     - \code{oddsRatioUpper}: The upper equivalence limit of odds ratio.
#'
#'     - \code{pi1}: The assumed probability for the active treatment group.
#'
#'     - \code{pi2}: The assumed probability for the control group.
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
#'     - \code{cumulativeAttainedAlphaH10}: The cumulative alpha attained
#'       under H10.
#'
#'     - \code{cumulativeAttainedAlphaH20}: The cumulative alpha attained
#'       under H20.
#'
#'     - \code{efficacyOddsRatioLower}: The efficacy boundaries on the
#'       odds ratio scale for the one-sided null hypothesis on the
#'       lower equivalence limit.
#'
#'     - \code{efficacyOddsRatioUpper}: The efficacy boundaries on the
#'       odds ratio scale for the one-sided null hypothesis on the
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
#'     - \code{calculationTarget}: The calculation target, \code{beta} or
#'       \code{n}.
#'
#'     - \code{nullVariance}: Whether to use the variance under the null
#'       or the empirical variance under the alternative.
#'
#'     - \code{varianceRatioH10}: The ratio of the variance under H10 to
#'       the variance under H1.
#'
#'     - \code{varianceRatioH20}: The ratio of the variance under H20 to
#'       the variance under H1.
#'
#'     - \code{varianceRatioH12}: The ratio of the variance under H10 to
#'       the variance under H20.
#'
#'     - \code{varianceRatioH21}: The ratio of the variance under H20 to
#'       the variance under H10.
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
#' (design1 <- getDesignOddsRatioEquiv(
#'   beta = 0.2, n = NA, oddsRatioLower = 0.8,
#'   oddsRatioUpper = 1.25, pi1 = 0.12, pi2 = 0.12,
#'   kMax = 3, alpha = 0.05, typeAlphaSpending = "sfOF"))
#'
#' @export
getDesignOddsRatioEquiv <- function(
    beta = NA_real_,
    n = NA_real_,
    oddsRatioLower = NA_real_,
    oddsRatioUpper = NA_real_,
    pi1 = NA_real_,
    pi2 = NA_real_,
    nullVariance = FALSE,
    allocationRatioPlanned = 1,
    rounding = TRUE,
    kMax = 1L,
    informationRates = NA_real_,
    criticalValues = NA_real_,
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
    stop("beta must lie in [0.0001, 1-alpha)");
  }

  if (!is.na(n) && n <= 0) {
    stop("n must be positive")
  }

  if (is.na(oddsRatioLower)) {
    stop("oddsRatioLower must be provided")
  }

  if (is.na(oddsRatioUpper)) {
    stop("oddsRatioUpper must be provided")
  }

  if (oddsRatioLower <= 0) {
    stop("oddsRatioLower must be positive")
  }

  if (pi1 <= 0 || pi1 >= 1) {
    stop("pi1 must lie between 0 and 1")
  }

  if (pi2 <= 0 || pi2 >= 1) {
    stop("pi2 must lie between 0 and 1")
  }

  oddsRatio = pi1*(1-pi2)/((1-pi1)*pi2)

  if (oddsRatioLower >= oddsRatio) {
    stop("oddsRatioLower must be less than pi1*(1-pi2)/((1-pi1)*pi2)")
  }

  if (oddsRatioUpper <= oddsRatio) {
    stop("oddsRatioUpper must be greater than pi1*(1-pi2)/((1-pi1)*pi2)")
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

  # variance for one sampling unit
  vH1 = 1/(r*pi1*(1-pi1)) + 1/((1-r)*pi2*(1-pi2))

  if (nullVariance) {
    # variance for one sampling unit under H10
    pi1H10 = pi2*oddsRatioLower/(1 + pi2*(oddsRatioLower - 1))
    vH10 = 1/(r*pi1H10*(1-pi1H10)) + 1/((1-r)*pi2*(1-pi2))

    # variance for one sampling unit under H20
    pi1H20 = pi2*oddsRatioUpper/(1 + pi2*(oddsRatioUpper - 1))
    vH20 = 1/(r*pi1H20*(1-pi1H20)) + 1/((1-r)*pi2*(1-pi2))

    # restricted maximum likelihood estimates under H10 w.r.t. H1
    mr = remlOddsRatio(oddsRatioLower, r, r*pi1, 1-r, (1-r)*pi2)
    vH10H1 = 1/(r*mr[1]*(1-mr[1])) + 1/((1-r)*mr[2]*(1-mr[2]))

    # restricted maximum likelihood estimates under H20 w.r.t. H1
    mr = remlOddsRatio(oddsRatioUpper, r, r*pi1, 1-r, (1-r)*pi2)
    vH20H1 = 1/(r*mr[1]*(1-mr[1])) + 1/((1-r)*mr[2]*(1-mr[2]))

    # restricted maximum likelihood estimates under H10 w.r.t. H20
    mr = remlOddsRatio(oddsRatioLower, r, r*pi1H20, 1-r, (1-r)*pi2)
    vH10H20 = 1/(r*mr[1]*(1-mr[1])) + 1/((1-r)*mr[2]*(1-mr[2]))

    # restricted maximum likelihood estimates under H20 w.r.t. H10
    mr = remlOddsRatio(oddsRatioUpper, r, r*pi1H10, 1-r, (1-r)*pi2)
    vH20H10 = 1/(r*mr[1]*(1-mr[1])) + 1/((1-r)*mr[2]*(1-mr[2]))

    # calculate the variance ratios
    varianceRatioH10 = vH10H1/vH1
    varianceRatioH20 = vH20H1/vH1
    varianceRatioH12 = vH10H20/vH20
    varianceRatioH21 = vH20H10/vH10
  } else {
    varianceRatioH10 = varianceRatioH20 = 1
    varianceRatioH12 = varianceRatioH21 = 1
  }

  # sample size calculation
  if (is.na(n)) {
    des = getDesignEquiv(
      beta = beta, IMax = NA, thetaLower = log(oddsRatioLower),
      thetaUpper = log(oddsRatioUpper), theta = log(oddsRatio),
      kMax = kMax, informationRates = informationRates,
      alpha = alpha, typeAlphaSpending = typeAlphaSpending,
      parameterAlphaSpending = parameterAlphaSpending,
      userAlphaSpending = userAlphaSpending,
      spendingTime = spendingTime,
      varianceRatioH10 = varianceRatioH10,
      varianceRatioH20 = varianceRatioH20,
      varianceRatioH12 = varianceRatioH12,
      varianceRatioH21 = varianceRatioH21)

    n = des$overallResults$information*vH1
  }

  if (rounding) {
    n = ceiling(n)
    informationRates = round(n*informationRates)/n
  }

  if (is.na(beta) || rounding) { # calculate power
    des = getDesignEquiv(
      beta = NA, IMax = n/vH1, thetaLower = log(oddsRatioLower),
      thetaUpper = log(oddsRatioUpper), theta = log(oddsRatio),
      kMax = kMax, informationRates = informationRates,
      alpha = alpha, typeAlphaSpending = typeAlphaSpending,
      parameterAlphaSpending = parameterAlphaSpending,
      userAlphaSpending = userAlphaSpending,
      spendingTime = spendingTime,
      varianceRatioH10 = varianceRatioH10,
      varianceRatioH20 = varianceRatioH20,
      varianceRatioH12 = varianceRatioH12,
      varianceRatioH21 = varianceRatioH21)
  }

  des$overallResults$numberOfSubjects = n
  des$overallResults$expectedNumberOfSubjectsH1 =
    des$overallResults$expectedInformationH1*vH1
  des$overallResults$expectedNumberOfSubjectsH10 =
    des$overallResults$expectedInformationH10*vH1
  des$overallResults$expectedNumberOfSubjectsH20 =
    des$overallResults$expectedInformationH20*vH1
  des$overallResults$oddsRatioLower = oddsRatioLower
  des$overallResults$oddsRatioUpper = oddsRatioUpper
  des$overallResults$pi1 = pi1
  des$overallResults$pi2 = pi2
  des$overallResults <-
    des$overallResults[, c("overallReject", "alpha",
                           "attainedAlphaH10", "attainedAlphaH20",
                           "kMax", "information",
                           "expectedInformationH1",
                           "expectedInformationH10",
                           "expectedInformationH20",
                           "numberOfSubjects",
                           "expectedNumberOfSubjectsH1",
                           "expectedNumberOfSubjectsH10",
                           "expectedNumberOfSubjectsH20",
                           "oddsRatioLower", "oddsRatioUpper",
                           "pi1", "pi2")]

  des$byStageResults$efficacyOddsRatioLower =
    exp(des$byStageResults$efficacyThetaLower)
  des$byStageResults$efficacyOddsRatioUpper =
    exp(des$byStageResults$efficacyThetaUpper)
  des$byStageResults$numberOfSubjects = n*informationRates
  des$byStageResults <-
    des$byStageResults[, c("informationRates", "efficacyBounds",
                           "rejectPerStage", "cumulativeRejection",
                           "cumulativeAlphaSpent",
                           "cumulativeAttainedAlphaH10",
                           "cumulativeAttainedAlphaH20",
                           "efficacyOddsRatioLower",
                           "efficacyOddsRatioUpper", "efficacyP",
                           "information", "numberOfSubjects")]

  des$settings$nullVariance = nullVariance
  des$settings$allocationRatioPlanned = allocationRatioPlanned
  des$settings$rounding = rounding
  des$settings$calculationTarget = ifelse(is.na(beta), "beta", "n")
  des$settings <-
    des$settings[c("typeAlphaSpending", "parameterAlphaSpending",
                   "userAlphaSpending", "spendingTime",
                   "calculationTarget", "nullVariance",
                   "varianceRatioH10", "varianceRatioH20",
                   "varianceRatioH12", "varianceRatioH21",
                   "allocationRatioPlanned", "rounding")]

  attr(des, "class") = "designOddsRatioEquiv"

  des
}


#' @title Power and sample size for Fisher's exact test for two proportions
#' @description Obtains the power given sample size or obtains the sample
#' size given power for Fisher's exact test for two proportions.
#'
#' @param beta The type II error.
#' @param n The total sample size.
#' @param pi1 The assumed probability for the active treatment group.
#' @param pi2 The assumed probability for the control group.
#' @param allocationRatioPlanned Allocation ratio for the active treatment
#'   versus control. Defaults to 1 for equal randomization.
#' @param alpha The two-sided significance level. Defaults to 0.05.
#'
#' @return A data frame with the following variables:
#'
#' * \code{alpha}: The two-sided significance level.
#'
#' * \code{power}: The power.
#'
#' * \code{n}: The sample size.
#'
#' * \code{pi1}: The assumed probability for the active treatment group.
#'
#' * \code{pi2}: The assumed probability for the control group.
#'
#' * \code{calculationTarget}: The calculation target, \code{beta} or
#'   \code{n}.
#'
#' * \code{allocationRatioPlanned}: Allocation ratio for the active
#'   treatment versus control.
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @examples
#'
#' (design1 <- getDesignFisherExact(
#'   beta = 0.2, pi1 = 0.5, pi2 = 0.2, alpha = 0.05))
#'
#' @export
getDesignFisherExact <- function(
    beta = NA_real_,
    n = NA_real_,
    pi1 = NA_real_,
    pi2 = NA_real_,
    allocationRatioPlanned = 1,
    alpha = 0.05) {

  if (is.na(beta) && is.na(n)) {
    stop("beta and n cannot be both missing")
  }

  if (!is.na(beta) && !is.na(n)) {
    stop("Only one of beta and n should be provided")
  }

  if (!is.na(beta) && (beta >= 1-alpha || beta < 0.0001)) {
    stop("beta must lie in [0.0001, 1-alpha)");
  }

  if (!is.na(n) && n <= 0) {
    stop("n must be positive")
  }

  if (pi1 <= 0 || pi1 >= 1) {
    stop("pi1 must lie between 0 and 1")
  }

  if (pi2 <= 0 || pi2 >= 1) {
    stop("pi2 must lie between 0 and 1")
  }

  if (allocationRatioPlanned <= 0) {
    stop("allocationRatioPlanned must be positive")
  }

  if (alpha < 0.00001 || alpha >= 1) {
    stop("alpha must lie in [0.00001, 1)")
  }

  if (!is.na(n)) { # power calculation
    a = powerFisherExact(n, pi1, pi2, allocationRatioPlanned, alpha)
  } else { # sample size calculation
    a = samplesizeFisherExact(beta, pi1, pi2, allocationRatioPlanned, alpha)
  }

  a$calculationTarget = ifelse(is.na(beta), "beta", "n")

  a
}


#' @title Clopper-Pearson confidence interval for one-sample proportion
#' @description Obtains the Clopper-Pearson exact confidence interval for
#' a one-sample proportion.
#'
#' @param n The sample size.
#' @param y The number of responses.
#' @param cilevel The confidence interval level.
#'
#' @return A data frame with the following variables:
#'
#' * \code{n}: The sample size.
#'
#' * \code{y}: The number of responses.
#'
#' * \code{phat}: The observed proportion of responses.
#'
#' * \code{lower}: The lower limit of the confidence interval.
#'
#' * \code{upper}: The upper limit of the confidence interval.
#'
#' * \code{cilevel}: The confidence interval level.
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @examples
#'
#' ClopperPearsonCI(20, 3)
#'
#' @export
#'
ClopperPearsonCI <- function(n, y, cilevel = 0.95) {
  if (n <= 0 || n != round(n)) {
    stop("n must be a positive integer")
  }

  if (y < 0 || y > n || y != round(y)) {
    stop("y must be an integer between 0 and n, inclusive")
  }

  if (cilevel <= 0 || cilevel >= 1) {
    stop("cilevel must lie between 0 and 1")
  }

  alpha = 1 - cilevel
  lower = 0
  upper = 1
  if (y > 0) lower = qbeta(alpha/2, y, n-y+1)
  if (y < n) upper = qbeta(1-alpha/2, y+1, n-y)

  data.frame(n = n, y = y, phat = y/n,
             lower = lower, upper = upper, cilevel = cilevel)
}


#' @title mTPI-2 decision table
#' @description Obtains the decision table for the modified toxicity
#' probability interval-2 (mTPI-2) design.
#'
#' @param nMax The maximum number of subjects in a dose cohort.
#' @param pT The target toxicity probability. Defaults to 0.3.
#' @param epsilon1 The lower equivalence margin from the target.
#'   Defaults to 0.05.
#' @param epsilon2 The upper equivalence margin from the target.
#'   Defaults to 0.05.
#' @param a The prior toxicity parameter for the beta prior.
#' @param b The prior non-toxicity parameter for the beta prior.
#' @param pExcessTox The threshold for excessive toxicity, i.e.,
#'   if Prob(p > pT | Data) > pExcessTox, then the current and
#'   all higher doses will be excluded and never be used again
#'   in the remainder of the trial to avoid any other subjects
#'   receiving treatment at those doses. Defaults to 0.95.
#'
#' @return An S3 class \code{mTPI2Table} object with the following
#' components:
#'
#' * \code{settings}: The input settings data frame with the following
#'   variables:
#'
#'     - \code{nMax}: The maximum number of subjects in a dose cohort.
#'
#'     - \code{pT}: The target toxicity probability.
#'
#'     - \code{epsilon1}: The lower equivalence margin from the target.
#'
#'     - \code{epsilon2}: The upper equivalence margin from the target.
#'
#'     - \code{a}: The prior toxicity parameter for the beta prior.
#'
#'     - \code{b}: The prior non-toxicity parameter for the beta prior.
#'
#'     - \code{pExcessTox}: The threshold for excessive toxicity.
#'
#' * \code{subintervals}: The subintervals of equal length in the mTPI-2
#'   design. It includes the following variables:
#'
#'     - \code{lower}: The lower bound of the subinterval.
#'
#'     - \code{upper}: The upper bound of the subinterval.
#'
#'     - \code{decision}: The dosing decision for the subinterval.
#'
#' * \code{decisionDataFrame}: The decision data frame for the mTPI-2 design.
#'   It includes the following variables:
#'
#'     - \code{n}: The sample size.
#'
#'     - \code{y}: The number of toxicities.
#'
#'     - \code{decision}: The dosing decision.
#'
#' * \code{decisionMatrix}: The decision matrix corresponding to the
#'   decision data frame.
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @examples
#'
#' mTPI2Table(nMax = 18, pT = 0.3, epsilon1 = 0.05, epsilon2 = 0.05)
#'
#' @export
#'
mTPI2Table <- function(
    nMax = NA_integer_, pT = 0.3,
    epsilon1 = 0.05, epsilon2 = 0.05,
    a = 1, b = 1, pExcessTox = 0.95) {

  if (nMax <= 0 || nMax != round(nMax)) {
    stop("nMax must be a positive integer")
  }

  if (pT <= 0 || pT >= 1) {
    stop("pT must lie between 0 and 1")
  }

  if (epsilon1 <= 0) {
    stop("epsilon1 must be positive")
  }

  if (pT - epsilon1 <= 0) {
    stop("epsilon1 must be less than pT")
  }

  if (epsilon2 <= 0) {
    stop("epsilon2 must be positive")
  }

  if (1 - pT - epsilon2 <= 0) {
    stop("epsilon2 must be less than 1 - pT")
  }

  if (a <= 0) {
    stop("a must be positive")
  }

  if (b <= 0) {
    stop("b must be positive")
  }

  if (pExcessTox <= 0 || pExcessTox >= 1) {
    stop("pExcessTox must lie between 0 and 1")
  }


  settings = data.frame(nMax, pT, epsilon1, epsilon2, a, b, pExcessTox)

  # construct the data frame of subintervals with equal length
  width = epsilon1 + epsilon2
  l1 = floor((pT - epsilon1)/width)
  delta1 = pT - epsilon1 - l1*width
  l2 = floor((1 - pT - epsilon2)/width)
  delta2 = 1 - pT - epsilon2 - l2*width

  df1 = data.frame(
    lower = c(0, delta1 + (0:l1)*width, 1 - delta2 - (l2:0)*width),
    upper = c(delta1 + (0:l1)*width, 1 - delta2 - (l2:0)*width, 1),
    decision = c(rep("E", l1+1), "S", rep("D", l2+1)))

  # construct the decision table
  df2 = data.frame()
  for (n in 1:nMax) {
    for (y in 0:n) {
      # check the exclusion rule
      if (1 - pbeta(pT, a+y, b+n-y) > pExcessTox & y != 1) {
        decision = "DU"
      } else { # unit probability mass
        num = pbeta(df1$upper, a+y, b+n-y) - pbeta(df1$lower, a+y, b+n-y)
        den = df1$upper - df1$lower
        postprob = num/den
        i = which.max(postprob)
        decision = df1$decision[i]
      }
      df2 = rbind(df2, data.frame(n = n, y = y, decision = decision))
    }
  }


  # generate the plot for the decision table
  df3 = matrix("", nrow = nMax+1, ncol = nMax,
               dimnames = list(as.character(0:nMax),
                               as.character(1:nMax)))
  for (i in 1:nrow(df2)) {
    df3[df2$y[i]+1, df2$n[i]] = df2$decision[i]
  }

  des = list(settings = settings,
             subintervals = df1,
             decisionDataFrame = df2,
             decisionMatrix = as.data.frame(df3))

  attr(des, "class") = "mTPI2Table"

  des
}


#' @title BOIN decision table
#' @description Obtains the decision table for the Bayesian optimal
#' interval (BOIN) design.
#'
#' @param nMax The maximum number of subjects in a dose cohort.
#' @param pT The target toxicity probability. Defaults to 0.3.
#' @param phi1 The lower equivalence limit for target toxicity
#'   probability.
#' @param phi2 The upper equivalence limit for target toxicity
#'   probability.
#' @param a The prior toxicity parameter for the beta prior.
#' @param b The prior non-toxicity parameter for the beta prior.
#' @param pExcessTox The threshold for excessive toxicity, i.e.,
#'   if Prob(p > pT | Data) > pExcessTox, then the current and
#'   all higher doses will be excluded and never be used again
#'   in the remainder of the trial to avoid any other subjects
#'   receiving treatment at those doses. Defaults to 0.95.
#'
#' @return An S3 class \code{BOINTable} object with the following
#' components:
#'
#' * \code{settings}: The input settings data frame with the following
#'   variables:
#'
#'     - \code{nMax}: The maximum number of subjects in a dose cohort.
#'
#'     - \code{pT}: The target toxicity probability.
#'
#'     - \code{phi1}: The lower equivalence limit for target toxicity
#'       probability.
#'
#'     - \code{phi2}: The upper equivalence limit for target toxicity
#'       probability.
#'
#'     - \code{lambda1}: The lower decision boundary for observed toxicity
#'       probability.
#'
#'     - \code{lambda2}: The upper decision boundary for observed toxicity
#'       probability.
#'
#'     - \code{a}: The prior toxicity parameter for the beta prior.
#'
#'     - \code{b}: The prior non-toxicity parameter for the beta prior.
#'
#'     - \code{pExcessTox}: The threshold for excessive toxicity.
#'
#' * \code{decisionDataFrame}: The decision data frame for the BOIN design.
#'   It includes the following variables:
#'
#'     - \code{n}: The sample size.
#'
#'     - \code{y}: The number of toxicities.
#'
#'     - \code{decision}: The dosing decision.
#'
#' * \code{decisionMatrix}: The decision matrix corresponding to the
#'   decision data frame.
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @examples
#'
#' BOINTable(nMax = 18, pT = 0.3, phi = 0.6*0.3, phi2 = 1.4*0.3)
#'
#' @export
#'
BOINTable <- function(
    nMax = NA_integer_, pT = 0.3,
    phi1 = 0.6*pT, phi2 = 1.4*pT,
    a = 1, b = 1, pExcessTox = 0.95) {

  if (nMax <= 0 || nMax != round(nMax)) {
    stop("nMax must be a positive integer")
  }

  if (pT <= 0 || pT >= 1) {
    stop("pT must lie between 0 and 1")
  }

  if (phi1 <= 0 || phi1 >= pT) {
    stop("phi1 must lie between 0 and", pT)
  }

  if (phi2 <= pT || phi2 >= 1) {
    stop("phi2 must lies between", pT, "and 1")
  }

  if (a <= 0) {
    stop("a must be positive")
  }

  if (b <= 0) {
    stop("b must be positive")
  }

  if (pExcessTox <= 0 || pExcessTox >= 1) {
    stop("pExcessTox must lie between 0 and 1")
  }

  lambda1 = log((1-phi1)/(1-pT))/log(pT*(1-phi1)/(phi1*(1-pT)))
  lambda2 = log((1-pT)/(1-phi2))/log(phi2*(1-pT)/(pT*(1-phi2)))

  settings = data.frame(nMax, pT, phi1 = phi1, phi2 = phi2,
                        lambda1 = lambda1, lambda2 = lambda2,
                        a, b, pExcessTox)

  # construct the decision table
  df2 = data.frame()
  for (n in 1:nMax) {
    for (y in 0:n) {
      # check the exclusion rule
      if (1 - pbeta(pT, a+y, b+n-y) > pExcessTox & n > 2) {
        decision = "DU"
      } else {
        p = y/n
        if (p <= lambda1) {
          decision = "E"
        } else if (p < lambda2) {
          decision = "S"
        } else {
          decision = "D"
        }
      }
      df2 = rbind(df2, data.frame(n = n, y = y, decision = decision))
    }
  }


  # generate the plot for the decision table
  df3 = matrix("", nrow = nMax+1, ncol = nMax,
               dimnames = list(as.character(0:nMax),
                               as.character(1:nMax)))
  for (i in 1:nrow(df2)) {
    df3[df2$y[i]+1, df2$n[i]] = df2$decision[i]
  }

  des = list(settings = settings,
             decisionDataFrame = df2,
             decisionMatrix = as.data.frame(df3))

  attr(des, "class") = "BOINTable"

  des
}


#' @title Power and sample for one-sample multinomial response
#' @description Obtains the power given sample size or obtains the sample
#' size given power for one-sample multinomial response.
#'
#' @param beta The type II error.
#' @param n The total sample size.
#' @param ncats The number of categories of the multinomial response.
#' @param piH0 The prevalence of each category under the null hypothesis.
#'   Only need to provide the values for the first \code{ncats-1}
#'   categories.
#' @param pi The prevalence of each category. Only need to provide the
#'   values for the first \code{ncats-1} categories.
#' @param rounding Whether to round up sample size. Defaults to 1 for
#'   sample size rounding.
#' @param alpha The two-sided significance level. Defaults to 0.05.
#'
#' @return An S3 class \code{designOneMultinom} object with the following
#' components:
#'
#' * \code{power}: The power to reject the null hypothesis.
#'
#' * \code{alpha}: The two-sided significance level.
#'
#' * \code{n}: The maximum number of subjects.
#'
#' * \code{ncats}: The number of categories of the multinomial response.
#'
#' * \code{piH0}: The prevalence of each category under the null hypothesis.
#'
#' * \code{pi}: The prevalence of each category.
#'
#' * \code{effectsize}: The effect size for the chi-square test.
#'
#' * \code{calculationTarget}: The calculation target, \code{beta} or
#'   \code{n}.
#'
#' * \code{rounding}: Whether to round up sample size.
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @examples
#'
#' (design1 <- getDesignOneMultinom(
#'   beta = 0.1, ncats = 3, piH0 = c(0.25, 0.25),
#'   pi = c(0.3, 0.4), alpha = 0.05))
#'
#' @export
#'
getDesignOneMultinom <- function(
    beta = NA_real_,
    n = NA_real_,
    ncats = NA_integer_,
    piH0 = NA_real_,
    pi = NA_real_,
    rounding = TRUE,
    alpha = 0.05) {

  if (is.na(beta) && is.na(n)) {
    stop("beta and n cannot be both missing")
  }

  if (!is.na(beta) && !is.na(n)) {
    stop("Only one of beta and n should be provided")
  }

  if (!is.na(beta) && (beta >= 1-alpha || beta < 0.0001)) {
    stop("beta must lie in [0.0001, 1-alpha)");
  }

  if (!is.na(n) && (n <= 0 || n != round(n))) {
    stop("n must be a positive integer")
  }

  if (is.na(ncats)) {
    stop("ncats must be provided")
  }

  if (ncats < 2 || ncats != round(ncats)) {
    stop("ncats must be a positive integer >= 2")
  }

  if (any(is.na(piH0))) {
    stop("piH0 must be provided")
  }

  if (any(piH0 <= 0 | piH0 >= 1)) {
    stop(paste("Elements of piH0 must lie between 0 and 1"))
  }

  if (length(piH0) != ncats-1 && length(piH0) != ncats) {
    stop(paste("piH0 must have length", ncats-1, "or", ncats))
  }

  if (length(piH0) == ncats-1) {
    if (sum(piH0) >= 1) {
      stop("piH0 must sum up to < 1 if it has only", ncats-1, "elements")
    }

    piH0 = c(piH0, 1 - sum(piH0))
  } else {
    if (sum(piH0) != 1) {
      stop(paste("piH0 must sum to 1 across", ncats, "categories"))
    }
  }

  if (any(is.na(pi))) {
    stop("pi must be provided")
  }

  if (any(pi <= 0 | pi >= 1)) {
    stop(paste("Elements of pi must lie between 0 and 1"))
  }

  if (length(pi) != ncats-1 && length(pi) != ncats) {
    stop(paste("pi must have length", ncats-1, "or", ncats))
  }

  if (length(pi) == ncats-1) {
    if (sum(pi) >= 1) {
      stop("pi must sum up to < 1 if it has only", ncats-1, "elements")
    }

    pi = c(pi, 1 - sum(pi))
  } else {
    if (sum(pi) != 1) {
      stop(paste("pi must sum to 1 across", ncats, "categories"))
    }
  }

  if (alpha < 0.00001 || alpha >= 1) {
    stop("alpha must lie in [0.00001, 1)")
  }

  tau = sum((pi - piH0)^2/piH0)

  nu = ncats - 1
  b = qchisq(1-alpha, nu)

  if (is.na(n)) {
    n0 = (b - nu + qnorm(1-beta)*sqrt(2*nu))/tau
    while (pchisq(b, nu, n0*tau) > beta) n0 <- 2*n0
    n = uniroot(function(n) pchisq(b, nu, n*tau) - beta, c(0.5*n0, n0))$root
  }

  if (rounding) {
    n = ceiling(n)
  }

  power = pchisq(b, nu, n*tau, lower.tail = FALSE)

  des = list(power = power, alpha = alpha, n = n,
             ncats = ncats, piH0 = piH0, pi = pi,
             effectsize = tau,
             calculationTarget = ifelse(is.na(beta), "beta", "n"),
             rounding = rounding)

  attr(des, "class") = "designOneMultinom"

  des
}


#' @title Power and sample for difference in two-sample multinomial response
#' @description Obtains the power given sample size or obtains the sample
#' size given power for difference in two-sample multinomial response.
#'
#' @param beta The type II error.
#' @param n The total sample size.
#' @param ncats The number of categories of the multinomial response.
#' @param pi1 The prevalence of each category for the treatment group.
#'   Only need to specify the valued for the first \code{ncats-1} categories.
#' @param pi2 The prevalence of each category for the control group.
#'   Only need to specify the valued for the first \code{ncats-1} categories.
#' @param allocationRatioPlanned Allocation ratio for the active treatment
#'   versus control. Defaults to 1 for equal randomization.
#' @param rounding Whether to round up sample size. Defaults to 1 for
#'   sample size rounding.
#' @param alpha The two-sided significance level. Defaults to 0.05.
#'
#' @return An S3 class \code{designTwoMultinom} object with the following
#' components:
#'
#' * \code{power}: The power to reject the null hypothesis.
#'
#' * \code{alpha}: The two-sided significance level.
#'
#' * \code{n}: The maximum number of subjects.
#'
#' * \code{ncats}: The number of categories of the multinomial response.
#'
#' * \code{pi1}: The prevalence of each category for the treatment group.
#'
#' * \code{pi2}: The prevalence of each category for the control group.
#'
#' * \code{effectsize}: The effect size for the chi-square test.
#'
#' * \code{calculationTarget}: The calculation target, \code{beta} or
#'   \code{n}.
#'
#' * \code{allocationRatioPlanned}: Allocation ratio for the active treatment
#'   versus control.
#'
#' * \code{rounding}: Whether to round up sample size.
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @examples
#'
#' (design1 <- getDesignTwoMultinom(
#'   beta = 0.1, ncats = 3, pi1 = c(0.3, 0.35),
#'   pi2 = c(0.2, 0.3), alpha = 0.05))
#'
#' @export
#'
getDesignTwoMultinom <- function(
    beta = NA_real_,
    n = NA_real_,
    ncats = NA_integer_,
    pi1 = NA_real_,
    pi2 = NA_real_,
    allocationRatioPlanned = 1,
    rounding = TRUE,
    alpha = 0.05) {

  if (is.na(beta) && is.na(n)) {
    stop("beta and n cannot be both missing")
  }

  if (!is.na(beta) && !is.na(n)) {
    stop("Only one of beta and n should be provided")
  }

  if (!is.na(beta) && (beta >= 1-alpha || beta < 0.0001)) {
    stop("beta must lie in [0.0001, 1-alpha)");
  }

  if (!is.na(n) && (n <= 0 || n != round(n))) {
    stop("n must be a positive integer")
  }

  if (is.na(ncats)) {
    stop("ncats must be provided")
  }

  if (ncats < 2 || ncats != round(ncats)) {
    stop("ncats must be a positive integer >= 2")
  }

  if (any(is.na(pi1))) {
    stop("pi1 must be provided")
  }

  if (any(pi1 <= 0 | pi1 >= 1)) {
    stop(paste("Elements of pi1 must lie between 0 and 1"))
  }

  if (length(pi1) != ncats-1 && length(pi1) != ncats) {
    stop(paste("pi1 must have length", ncats-1, "or", ncats))
  }

  if (length(pi1) == ncats-1) {
    if (sum(pi1) >= 1) {
      stop("pi1 must sum up to < 1 if it has only", ncats-1, "elements")
    }

    pi1 = c(pi1, 1 - sum(pi1))
  } else {
    if (sum(pi1) != 1) {
      stop(paste("pi1 must sum to 1 across", ncats, "categories"))
    }
  }

  if (any(is.na(pi2))) {
    stop("pi2 must be provided")
  }

  if (any(pi2 <= 0 | pi2 >= 1)) {
    stop(paste("Elements of pi2 must lie between 0 and 1"))
  }

  if (length(pi2) != ncats-1 && length(pi2) != ncats) {
    stop(paste("pi2 must have length", ncats-1, "or", ncats))
  }

  if (length(pi2) == ncats-1) {
    if (sum(pi2) >= 1) {
      stop("pi2 must sum up to < 1 if it has only", ncats-1, "elements")
    }

    pi2 = c(pi2, 1 - sum(pi2))
  } else {
    if (sum(pi2) != 1) {
      stop(paste("pi2 must sum to 1 across", ncats, "categories"))
    }
  }

  if (allocationRatioPlanned <= 0) {
    stop("allocationRatioPlanned must be positive")
  }

  if (alpha < 0.00001 || alpha >= 1) {
    stop("alpha must lie in [0.00001, 1)")
  }


  r = allocationRatioPlanned/(1 + allocationRatioPlanned)

  tau = r*(1-r)*sum((pi1 - pi2)^2/(r*pi1 + (1-r)*pi2))

  nu = ncats - 1
  b = qchisq(1-alpha, nu)

  if (is.na(n)) {
    n0 = (b - nu + qnorm(1-beta)*sqrt(2*nu))/tau
    while (pchisq(b, nu, n0*tau) > beta) n0 <- 2*n0
    n = uniroot(function(n) pchisq(b, nu, n*tau) - beta, c(0.5*n0, n0))$root
  }

  if (rounding) {
    n = ceiling(n)
  }

  power = pchisq(b, nu, n*tau, lower.tail = FALSE)

  des = list(power = power, alpha = alpha, n = n,
             ncats = ncats, pi1 = pi1, pi2 = pi2,
             effectsize = tau,
             calculationTarget = ifelse(is.na(beta), "beta", "n"),
             allocationRatioPlanned = allocationRatioPlanned,
             rounding = rounding)

  attr(des, "class") = "designTwoMultinom"

  des
}


#' @title Power and sample size for the Wilcoxon test for two-sample
#' ordinal response
#' @description Obtains the power given sample size or obtains the sample
#' size given power for the Wilcoxon test for two-sample ordinal response.
#'
#' @param beta The type II error.
#' @param n The total sample size.
#' @param ncats The number of categories of the ordinal response.
#' @param pi1 The prevalence of each category for the treatment group.
#'   Only need to specify the valued for the first \code{ncats-1} categories.
#' @param pi2 The prevalence of each category for the control group.
#'   Only need to specify the valued for the first \code{ncats-1} categories.
#' @param allocationRatioPlanned Allocation ratio for the active treatment
#'   versus control. Defaults to 1 for equal randomization.
#' @param rounding Whether to round up sample size. Defaults to 1 for
#'   sample size rounding.
#' @inheritParams param_alpha
#'
#' @return An S3 class \code{designTwoOrdinal} object with the following
#' components:
#'
#' * \code{power}: The power to reject the null hypothesis.
#'
#' * \code{alpha}: The two-sided significance level.
#'
#' * \code{n}: The maximum number of subjects.
#'
#' * \code{ncats}: The number of categories of the ordinal response.
#'
#' * \code{pi1}: The prevalence of each category for the treatment group.
#'
#' * \code{pi2}: The prevalence of each category for the control group.
#'
#' * \code{meanscore1}: The mean midrank score for the treatment group.
#'
#' * \code{meanscore2}: The mean midrank score for the control group.
#'
#' * \code{calculationTarget}: The calculation target, \code{beta} or
#'   \code{n}.
#'
#' * \code{allocationRatioPlanned}: Allocation ratio for the active treatment
#'   versus control.
#'
#' * \code{rounding}: Whether to round up sample size.
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @examples
#'
#' (design1 <- getDesignTwoOrdinal(
#'   beta = 0.1, ncats = 4, pi1 = c(0.55, 0.3, 0.1),
#'   pi2 = c(0.214, 0.344, 0.251), alpha = 0.025))
#'
#' @export
getDesignTwoOrdinal <- function(
    beta = NA_real_,
    n = NA_real_,
    ncats = NA_integer_,
    pi1 = NA_real_,
    pi2 = NA_real_,
    allocationRatioPlanned = 1,
    rounding = TRUE,
    alpha = 0.05) {

  if (is.na(beta) && is.na(n)) {
    stop("beta and n cannot be both missing")
  }

  if (!is.na(beta) && !is.na(n)) {
    stop("Only one of beta and n should be provided")
  }

  if (!is.na(beta) && (beta >= 1-alpha || beta < 0.0001)) {
    stop("beta must lie in [0.0001, 1-alpha)");
  }

  if (!is.na(n) && (n <= 0 || n != round(n))) {
    stop("n must be a positive integer")
  }

  if (is.na(ncats)) {
    stop("ncats must be provided")
  }

  if (ncats < 2 || ncats != round(ncats)) {
    stop("ncats must be a positive integer >= 2")
  }

  if (any(is.na(pi1))) {
    stop("pi1 must be provided")
  }

  if (any(pi1 <= 0 | pi1 >= 1)) {
    stop(paste("Elements of pi1 must lie between 0 and 1"))
  }

  if (length(pi1) != ncats-1 && length(pi1) != ncats) {
    stop(paste("pi1 must have length", ncats-1, "or", ncats))
  }

  if (length(pi1) == ncats-1) {
    if (sum(pi1) >= 1) {
      stop("pi1 must sum up to < 1 if it has only", ncats-1, "elements")
    }

    pi1 = c(pi1, 1 - sum(pi1))
  } else {
    if (sum(pi1) != 1) {
      stop(paste("pi1 must sum to 1 across", ncats, "categories"))
    }
  }

  if (any(is.na(pi2))) {
    stop("pi2 must be provided")
  }

  if (any(pi2 <= 0 | pi2 >= 1)) {
    stop(paste("Elements of pi2 must lie between 0 and 1"))
  }

  if (length(pi2) != ncats-1 && length(pi2) != ncats) {
    stop(paste("pi2 must have length", ncats-1, "or", ncats))
  }

  if (length(pi2) == ncats-1) {
    if (sum(pi2) >= 1) {
      stop("pi2 must sum up to < 1 if it has only", ncats-1, "elements")
    }

    pi2 = c(pi2, 1 - sum(pi2))
  } else {
    if (sum(pi2) != 1) {
      stop(paste("pi2 must sum to 1 across", ncats, "categories"))
    }
  }

  if (allocationRatioPlanned <= 0) {
    stop("allocationRatioPlanned must be positive")
  }

  if (alpha < 0.00001 || alpha >= 1) {
    stop("alpha must lie in [0.00001, 1)")
  }


  r = allocationRatioPlanned/(1 + allocationRatioPlanned)

  pi = r*pi1 + (1-r)*pi2
  p = cumsum(pi)
  w = p - 0.5*pi

  riskDiff = sum(w*(pi1 - pi2))

  directionUpper = riskDiff > 0

  theta = ifelse(directionUpper, riskDiff, -riskDiff)

  varH0 = (sum(w^2*pi) - sum(w*pi)^2)/(r*(1-r))
  varH1 = (sum(w^2*pi1) - sum(w*pi1)^2)/r +
    (sum(w^2*pi2) - sum(w*pi2)^2)/(1-r)


  if (is.na(n)) {
    n = (qnorm(1-alpha)*sqrt(varH0) + qnorm(1-beta)*sqrt(varH1))^2/theta^2
  }

  if (rounding) {
    n = ceiling(n)
  }

  power = pnorm((sqrt(n)*theta - qnorm(1-alpha)*sqrt(varH0))/sqrt(varH1))

  des = list(power = power, alpha = alpha, n = n,
             ncats = ncats, pi1 = pi1, pi2 = pi2,
             meanscore1 = n*sum(w*pi1), meanscore2 = n*sum(w*pi2),
             calculationTarget = ifelse(is.na(beta), "beta", "n"),
             allocationRatioPlanned = allocationRatioPlanned,
             rounding = rounding)

  attr(des, "class") = "designTwoOrdinal"

  des
}


#' @title Power and sample size for Cochran-Armitage trend test for
#' ordered multi-sample binomial response
#' @description Obtains the power given sample size or obtains the sample
#' size given power for the Cochran-Armitage trend test for ordered
#' multi-sample binomial response.
#'
#' @param beta The type II error.
#' @param n The total sample size.
#' @param ngroups The number of treatment groups.
#' @param pi The response probabilities for the treatment groups.
#' @param w The scores assigned to the treatment groups. This should
#'   reflect the ordinal nature of the treatment groups, e.g. dose levels.
#'   Defaults to equally spaced scores.
#' @param allocationRatioPlanned Allocation ratio for the treatment groups.
#' @param rounding Whether to round up sample size. Defaults to 1 for
#'   sample size rounding.
#' @param alpha The two-sided significance level. Defaults to 0.05.
#'
#' @return An S3 class \code{designOrderedBinom} object with the following
#' components:
#'
#' * \code{power}: The power to reject the null hypothesis.
#'
#' * \code{alpha}: The two-sided significance level.
#'
#' * \code{n}: The maximum number of subjects.
#'
#' * \code{ngroups}: The number of treatment groups.
#'
#' * \code{pi}: The response probabilities for the treatment groups.
#'
#' * \code{w}: The scores assigned to the treatment groups.
#'
#' * \code{trendstat}: The Cochran-Armitage trend test statistic.
#'
#' * \code{calculationTarget}: The calculation target, \code{beta} or
#'   \code{n}.
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
#' (design1 <- getDesignOrderedBinom(
#'   beta = 0.1, ngroups = 3, pi = c(0.1, 0.25, 0.5), alpha = 0.05))
#'
#' @export
getDesignOrderedBinom <- function(
    beta = NA_real_,
    n = NA_real_,
    ngroups = NA_integer_,
    pi = NA_real_,
    w = NA_real_,
    allocationRatioPlanned = NA_integer_,
    rounding = TRUE,
    alpha = 0.05) {

  if (is.na(beta) && is.na(n)) {
    stop("beta and n cannot be both missing")
  }

  if (!is.na(beta) && !is.na(n)) {
    stop("Only one of beta and n should be provided")
  }

  if (!is.na(beta) && (beta >= 1-alpha || beta < 0.0001)) {
    stop("beta must lie in [0.0001, 1-alpha)");
  }

  if (!is.na(n) && (n <= 0 || n != round(n))) {
    stop("n must be a positive integer")
  }

  if (is.na(ngroups)) {
    stop("ncats must be provided")
  }

  if (ngroups < 2 || ngroups != round(ngroups)) {
    stop("ngroups must be a positive integer >= 2")
  }

  if (any(is.na(pi))) {
    stop("pi must be provided")
  }

  if (any(pi <= 0 | pi >= 1)) {
    stop(paste("Elements of pi must lie between 0 and 1"))
  }

  if (length(pi) != ngroups) {
    stop(paste("pi must have length", ngroups))
  }

  if (any(is.na(w))) {
    w = seq(1, ngroups)
  }

  if (length(w) != ngroups) {
    stop(paste("w must have length", ngroups))
  }

  if (any(is.na(allocationRatioPlanned))) {
    allocationRatioPlanned = rep(1, ngroups)
  }

  if (length(allocationRatioPlanned) != ngroups-1 &&
      length(allocationRatioPlanned) != ngroups) {
    stop(paste("allocationRatioPlanned should have", ngroups-1,
               "or", ngroups, "elements"))
  }

  if (length(allocationRatioPlanned) == ngroups-1) {
    allocationRatioPlanned = c(allocationRatioPlanned, 1)
  }

  if (any(allocationRatioPlanned <= 0)) {
    stop("Elements of allocationRatioPlanned must be positive")
  }

  if (alpha < 0.00001 || alpha >= 1) {
    stop("alpha must lie in [0.00001, 1)")
  }


  # randomization probabilities
  r = allocationRatioPlanned/sum(allocationRatioPlanned)

  wbar = sum(r*w)
  pibar = sum(r*pi)

  theta = sum(r*(w-wbar)*(pi-pibar))

  varH0 = pibar*(1-pibar)*sum(r*(w-wbar)^2)
  varH1 = sum(r*(w-wbar)^2*pi*(1-pi))

  if (is.na(n)) {
    n = (qnorm(1-alpha/2)*sqrt(varH0) + qnorm(1-beta)*sqrt(varH1))^2/theta^2
  }

  if (rounding) {
    n = ceiling(n)
  }

  power = pnorm((sqrt(n)*abs(theta) - qnorm(1-alpha/2)*sqrt(varH0))/
                  sqrt(varH1))

  des = list(power = power, alpha = alpha, n = n,
             ngroups = ngroups, pi = pi, w = w,
             trendstat = n*sum(r*(w-wbar)*pi),
             calculationTarget = ifelse(is.na(beta), "beta", "n"),
             allocationRatioPlanned = allocationRatioPlanned,
             rounding = rounding)

  attr(des, "class") = "designOrderedBinom"

  des
}


#' @title Power and sample size for unordered multi-sample binomial response
#' @description Obtains the power given sample size or obtains the sample
#' size given power for the chi-square test for unordered multi-sample
#' binomial response.
#'
#' @param beta The type II error.
#' @param n The total sample size.
#' @param ngroups The number of treatment groups.
#' @param pi The response probabilities for the treatment groups.
#' @param allocationRatioPlanned Allocation ratio for the treatment groups.
#' @param rounding Whether to round up sample size. Defaults to 1 for
#'   sample size rounding.
#' @param alpha The two-sided significance level. Defaults to 0.05.
#'
#' @return An S3 class \code{designUnorderedBinom} object with the following
#' components:
#'
#' * \code{power}: The power to reject the null hypothesis.
#'
#' * \code{alpha}: The two-sided significance level.
#'
#' * \code{n}: The maximum number of subjects.
#'
#' * \code{ngroups}: The number of treatment groups.
#'
#' * \code{pi}: The response probabilities for the treatment groups.
#'
#' * \code{effectsize}: The effect size for the chi-square test.
#'
#' * \code{calculationTarget}: The calculation target, \code{beta} or
#'   \code{n}.
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
#' (design1 <- getDesignUnorderedBinom(
#'   beta = 0.1, ngroups = 3, pi = c(0.1, 0.25, 0.5), alpha = 0.05))
#'
#' @export
getDesignUnorderedBinom <- function(
    beta = NA_real_,
    n = NA_real_,
    ngroups = NA_integer_,
    pi = NA_real_,
    allocationRatioPlanned = NA_integer_,
    rounding = TRUE,
    alpha = 0.05) {

  if (is.na(beta) && is.na(n)) {
    stop("beta and n cannot be both missing")
  }

  if (!is.na(beta) && !is.na(n)) {
    stop("Only one of beta and n should be provided")
  }

  if (!is.na(beta) && (beta >= 1-alpha || beta < 0.0001)) {
    stop("beta must lie in [0.0001, 1-alpha)");
  }

  if (!is.na(n) && (n <= 0 || n != round(n))) {
    stop("n must be a positive integer")
  }

  if (is.na(ngroups)) {
    stop("ngroups must be provided")
  }

  if (ngroups < 2 || ngroups != round(ngroups)) {
    stop("ngroups must be a positive integer >= 2")
  }

  if (any(is.na(pi))) {
    stop("pi must be provided")
  }

  if (any(pi <= 0 | pi >= 1)) {
    stop(paste("Elements of pi must lie between 0 and 1"))
  }

  if (length(pi) != ngroups) {
    stop(paste("pi must have length", ngroups))
  }

  if (any(is.na(allocationRatioPlanned))) {
    allocationRatioPlanned = rep(1, ngroups)
  }

  if (length(allocationRatioPlanned) != ngroups-1 &&
      length(allocationRatioPlanned) != ngroups) {
    stop(paste("allocationRatioPlanned should have", ngroups-1,
               "or", ngroups, "elements"))
  }

  if (length(allocationRatioPlanned) == ngroups-1) {
    allocationRatioPlanned = c(allocationRatioPlanned, 1)
  }

  if (any(allocationRatioPlanned <= 0)) {
    stop("Elements of allocationRatioPlanned must be positive")
  }

  if (alpha < 0.00001 || alpha >= 1) {
    stop("alpha must lie in [0.00001, 1)")
  }


  # randomization probabilities
  r = allocationRatioPlanned/sum(allocationRatioPlanned)

  pibar = sum(r*pi)

  tau = sum(r*(pi - pibar)^2)/(pibar*(1-pibar))

  nu = ngroups - 1
  b = qchisq(1-alpha, nu)

  if (is.na(n)) {
    n0 = (qchisq(1-alpha, nu) - nu + qnorm(1-beta)*sqrt(2*nu))/tau
    while (pchisq(b, nu, n0*tau) > beta) n0 <- 2*n0
    n = uniroot(function(n) pchisq(b, nu, n*tau) - beta, c(0.5*n0, n0))$root
  }

  if (rounding) {
    n = ceiling(n)
  }

  power = pchisq(b, nu, n*tau, lower.tail = FALSE)

  des = list(power = power, alpha = alpha, n = n,
             ngroups = ngroups, pi = pi,
             effectsize = tau,
             calculationTarget = ifelse(is.na(beta), "beta", "n"),
             allocationRatioPlanned = allocationRatioPlanned,
             rounding = rounding)

  attr(des, "class") = "designUnorderedBinom"

  des
}


#' @title Power and sample size for unordered multi-sample
#' multinomial response
#' @description Obtains the power given sample size or obtains the sample
#' size given power for the chi-square test for unordered multi-sample
#' multinomial response.
#'
#' @param beta The type II error.
#' @param n The total sample size.
#' @param ngroups The number of treatment groups.
#' @param ncats The number of categories of the multinomial response.
#' @param pi The matrix of response probabilities for the treatment groups.
#'   It should have \code{ngroups} rows and \code{ncats-1} or \code{ncats}
#'   columns.
#' @param allocationRatioPlanned Allocation ratio for the treatment groups.
#' @param rounding Whether to round up sample size. Defaults to 1 for
#'   sample size rounding.
#' @param alpha The two-sided significance level. Defaults to 0.05.
#'
#' @return An S3 class \code{designUnorderedMultinom} object with the
#' following components:
#'
#' * \code{power}: The power to reject the null hypothesis.
#'
#' * \code{alpha}: The two-sided significance level.
#'
#' * \code{n}: The maximum number of subjects.
#'
#' * \code{ngroups}: The number of treatment groups.
#'
#' * \code{ncats}: The number of categories of the multinomial response.
#'
#' * \code{pi}: The response probabilities for the treatment groups.
#'
#' * \code{effectsize}: The effect size for the chi-square test.
#'
#' * \code{calculationTarget}: The calculation target, \code{beta} or
#'   \code{n}.
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
#' (design1 <- getDesignUnorderedMultinom(
#'   beta = 0.1, ngroups = 3, ncats = 4,
#'   pi = matrix(c(0.230, 0.320, 0.272,
#'                 0.358, 0.442, 0.154,
#'                 0.142, 0.036, 0.039),
#'               3, 3, byrow = TRUE),
#'   allocationRatioPlanned = c(2, 2, 1),
#'   alpha = 0.05))
#'
#' @export
getDesignUnorderedMultinom <- function(
    beta = NA_real_,
    n = NA_real_,
    ngroups = NA_integer_,
    ncats = NA_integer_,
    pi = NA_real_,
    allocationRatioPlanned = NA_integer_,
    rounding = TRUE,
    alpha = 0.05) {

  if (is.na(beta) && is.na(n)) {
    stop("beta and n cannot be both missing")
  }

  if (!is.na(beta) && !is.na(n)) {
    stop("Only one of beta and n should be provided")
  }

  if (!is.na(beta) && (beta >= 1-alpha || beta < 0.0001)) {
    stop("beta must lie in [0.0001, 1-alpha)");
  }

  if (!is.na(n) && (n <= 0 || n != round(n))) {
    stop("n must be a positive integer")
  }

  if (is.na(ngroups)) {
    stop("ngroups must be provided")
  }

  if (ngroups < 2 || ngroups != round(ngroups)) {
    stop("ngroups must be a positive integer >= 2")
  }

  if (is.na(ncats)) {
    stop("ncats must be provided")
  }

  if (ncats < 2 || ncats != round(ncats)) {
    stop("ncats must be a positive integer >= 2")
  }

  if (any(is.na(pi))) {
    stop("pi must be provided")
  }

  if (any(pi <= 0 | pi >= 1)) {
    stop(paste("Elements of pi must lie between 0 and 1"))
  }

  if (nrow(pi) != ngroups) {
    stop(paste("pi must have", ngroups, "rows"))
  }

  if (ncol(pi) != ncats-1 && ncol(pi) != ncats) {
    stop(paste("pi must have", ncats-1, "or", ncats, "columns"))
  }

  if (ncol(pi) == ncats-1) {
    if (any(rowSums(pi) >= 1)) {
      stop("pi should sum up to <1 for any treatment group")
    }

    pi = matrix(c(pi, 1-rowSums(pi)), ngroups, ncats)
  } else {
    if (any(rowSums(pi) != 1)) {
      stop("pi should sum up to 1 for any treatment group")
    }
  }

  if (any(is.na(allocationRatioPlanned))) {
    allocationRatioPlanned = rep(1, ngroups)
  }

  if (length(allocationRatioPlanned) != ngroups-1 &&
      length(allocationRatioPlanned) != ngroups) {
    stop(paste("allocationRatioPlanned should have", ngroups-1,
               "or", ngroups, "elements"))
  }

  if (length(allocationRatioPlanned) == ngroups-1) {
    allocationRatioPlanned = c(allocationRatioPlanned, 1)
  }

  if (any(allocationRatioPlanned <= 0)) {
    stop("Elements of allocationRatioPlanned must be positive")
  }

  if (alpha < 0.00001 || alpha >= 1) {
    stop("alpha must lie in [0.00001, 1)")
  }


  # randomization probabilities
  r = allocationRatioPlanned/sum(allocationRatioPlanned)
  mr = matrix(r, ngroups, ncats)

  pibar = colSums(r*pi)
  mpibar = matrix(pibar, ngroups, ncats, byrow = TRUE)

  tau = sum(mr*(pi - mpibar)^2/mpibar)

  nu = (ngroups-1)*(ncats-1)
  b = qchisq(1-alpha, nu)

  if (is.na(n)) {
    n0 = (qchisq(1-alpha, nu) - nu + qnorm(1-beta)*sqrt(2*nu))/tau
    while (pchisq(b, nu, n0*tau) > beta) n0 <- 2*n0
    n = uniroot(function(n) pchisq(b, nu, n*tau) - beta, c(0.5*n0, n0))$root
  }

  if (rounding) {
    n = ceiling(n)
  }

  power = pchisq(b, nu, n*tau, lower.tail = FALSE)

  des = list(power = power, alpha = alpha, n = n,
             ngroups = ngroups, ncats = ncats,
             pi = pi, effectsize = tau,
             calculationTarget = ifelse(is.na(beta), "beta", "n"),
             allocationRatioPlanned = allocationRatioPlanned,
             rounding = rounding)

  attr(des, "class") = "designUnorderedMultinom"

  des
}


#' @title Power and sample size for logistic regression
#' @description Obtains the power given sample size or obtains the sample
#' size given power for logistic regression of a binary response given
#' the covariate of interest and other covariates.
#'
#' @param beta The type II error.
#' @param n The total sample size.
#' @param ncovariates The number of covariates.
#' @param nconfigs The number of configurations of discretized covariate
#'   values.
#' @param x The matrix of covariate values.
#' @param pconfigs The vector of probabilities for the configurations.
#' @param corr The multiple correlation between the predictor and other
#'   covariates. Defaults to 0.
#' @param oddsratios The odds ratios for one unit increase in the
#'   covariates.
#' @param responseprob The response probability in the full model when
#'   all predictor variables are equal to their means.
#' @param rounding Whether to round up sample size. Defaults to 1 for
#'   sample size rounding.
#' @param alpha The two-sided significance level. Defaults to 0.05.
#'
#' @details
#' We consider the logistic regression of a binary response variable
#' \eqn{Y} on a set of predictor variables \eqn{x = (x_1,\ldots,x_K)^T}
#' with \eqn{x_1} being the covariate of interest:
#' \eqn{\log \frac{P(Y_i=1)}{1 - P(Y_i = 1)} = \psi_0 + x_i^T \psi,}
#' where \eqn{\psi = (\psi_1,\ldots,\psi_K)^T}.
#' Similar to Self et al (1992), we assume that all covariates are
#' either inherently discrete or discretized from continuous
#' distributions (e.g. using the quantiles). Let \eqn{m} denote the total
#' number of configurations of the covariate values. Let
#' \deqn{\pi_i = P(x = x_i), i = 1,\ldots, m}
#' denote the probabilities for the configurations of the covariates
#' under independence. The likelihood ratio test statistic for testing
#' \eqn{H_0: \psi_1 = 0} can be approximated by a noncentral chi-square
#' distribution with one degree of freedom and noncentrality parameter
#' \deqn{\Delta = 2 \sum_{i=1}^m \pi_i [b'(\theta_i)(\theta_i - \theta_i^*)
#' - \{b(\theta_i) - b(\theta_i^*)\}],} where
#' \deqn{\theta_i = \psi_0 + \sum_{j=1}^{k} \psi_j x_{ij},}
#' \deqn{\theta_i^* = \psi_0^* + \sum_{j=2}^{k} \psi_j^* x_{ij},}
#' for \eqn{\psi_0^* = \psi_0 + \psi_1 \mu_1}, and \eqn{\psi_j^* = \psi_j}
#' for \eqn{j=2,\ldots,K}. Here \eqn{\mu_1} is the mean of \eqn{x_1},
#' e.g., \eqn{\mu_1 = \sum_i \pi_i x_{i1}.} In addition, by
#' formulating the logistic regression in the framework of generalized
#' linear models, \deqn{b(\theta) = \log(1 + \exp(\theta)),} and
#' \deqn{b'(\theta) = \frac{\exp(\theta)}{1 + \exp(\theta)}.}
#'
#' The regression coefficients \eqn{\psi} can be obtained by taking the
#' log of the odds ratios for the covariates. The intercept \eqn{\psi_0}
#' can be derived as \deqn{\psi_0 = \log(\bar{\mu}/(1- \bar{\mu})) -
#' \sum_{j=1}^{K} \psi_j \mu_j,} where \eqn{\bar{\mu}} denotes the
#' response probability when all predictor variables are equal to their
#' means.
#'
#' Finally, let \eqn{\rho} denote the multiple correlation between
#' the predictor and other covariates. The noncentrality parameter
#' of the chi-square test is adjusted downward by multiplying by
#' \eqn{1-\rho^2}.
#'
#' @return An S3 class \code{designLogistic} object with the
#' following components:
#'
#' * \code{power}: The power to reject the null hypothesis.
#'
#' * \code{alpha}: The two-sided significance level.
#'
#' * \code{n}: The total sample size.
#'
#' * \code{ncovariates}: The number of covariates.
#'
#' * \code{nconfigs}: The number of configurations of discretized covariate
#'   values.
#'
#' * \code{x}: The matrix of covariate values.
#'
#' * \code{pconfigs}: The vector of probabilities for the configurations.
#'
#' * \code{corr}: The multiple correlation between the predictor and other
#'   covariates.
#'
#' * \code{oddsratios}: The odds ratios for one unit increase in the
#'   covariates.
#'
#' * \code{responseprob}: The response probability in the full model when
#'   all predictor variables are equal to their means.
#'
#' * \code{effectsize}: The effect size for the chi-square test.
#'
#' * \code{calculationTarget}: The calculation target, \code{beta} or
#'   \code{n}.
#'
#' * \code{rounding}: Whether to round up sample size.
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @references
#'
#' Steven G. Self, Robert H. Mauritsen and Jill Ohara.
#' Power calculations for likelihood ratio tests in generalized linear
#' models. Biometrics 1992; 48:31-39.
#'
#' @examples
#'
#' # two ordinal covariates
#' x1 = c(5, 10, 15, 20)
#' px1 = c(0.2, 0.3, 0.3, 0.2)
#'
#' x2 = c(2, 4, 6)
#' px2 = c(0.4, 0.4, 0.2)
#'
#' # discretizing a normal distribution with mean 4 and standard deviation 2
#' nbins = 10
#' x3 = qnorm(((1:nbins) - 0.5)/nbins)*2 + 4
#' px3 = rep(1/nbins, nbins)
#'
#' # combination of covariate values
#' nconfigs = length(x1)*length(x2)*length(x3)
#' x = expand.grid(x3 = x3, x2 = x2, x1 = x1)
#' x = as.matrix(x[, ncol(x):1])
#'
#' # probabilities for the covariate configurations under independence
#' pconfigs = as.numeric(px1 %x% px2 %x% px3)
#'
#' # convert the odds ratio for the predictor variable in 5-unit change
#' # to the odds ratio in 1-unit change
#' (design1 <- getDesignLogistic(
#'   beta = 0.1, ncovariates = 3,
#'   nconfigs = nconfigs,
#'   x = x,
#'   pconfigs = pconfigs,
#'   oddsratios = c(1.2^(1/5), 1.4, 1.3),
#'   responseprob = 0.25,
#'   alpha = 0.1))
#'
#' @export
#'
getDesignLogistic <- function(
    beta = NA_real_,
    n = NA_real_,
    ncovariates = NA_integer_,
    nconfigs = NA_integer_,
    x = NA_real_,
    pconfigs = NA_real_,
    corr = 0,
    oddsratios = NA_real_,
    responseprob = NA_real_,
    rounding = TRUE,
    alpha = 0.05) {

  if (is.na(beta) && is.na(n)) {
    stop("beta and n cannot be both missing")
  }

  if (!is.na(beta) && !is.na(n)) {
    stop("Only one of beta and n should be provided")
  }

  if (!is.na(beta) && (beta >= 1-alpha || beta < 0.0001)) {
    stop("beta must lie in [0.0001, 1-alpha)");
  }

  if (!is.na(n) && (n <= 0 || n != round(n))) {
    stop("n must be a positive integer")
  }

  if (is.na(ncovariates)) {
    stop("ncovariates must be provided")
  }

  if (ncovariates < 1 || ncovariates != round(ncovariates)) {
    stop("ncovariates must be a positive integer")
  }

  if (is.na(nconfigs)) {
    stop("nconfigs must be provided")
  }

  if (nconfigs < 2 || nconfigs != round(nconfigs)) {
    stop("nconfigs must be a positive integer >= 2")
  }

  if (any(is.na(x))) {
    stop("x must be provided")
  }

  if (ncovariates != 1 && !is.matrix(x)) {
    stop("x must be a matrix")
  }

  if (ncovariates == 1 && !is.matrix(x)) {
    x = as.matrix(x, ncol = 1)
  }

  if (is.matrix(x) && nrow(x) != nconfigs) {
    stop("x must have", nconfigs, "rows")
  }

  if (is.matrix(x) && ncol(x) != ncovariates) {
    stop("x must have", ncovariates, "columns")
  }

  if (any(is.na(pconfigs))) {
    stop("pconfigs must be provided")
  }

  if (any(pconfigs <= 0 | pconfigs >= 1)) {
    stop(paste("Elements of pconfigs must lie between 0 and 1"))
  }

  if (length(pconfigs) != nconfigs-1 && length(pconfigs) != nconfigs) {
    stop(paste("pconfigs must have", nconfigs-1, "or", nconfigs, "elements"))
  }

  if (length(pconfigs) == nconfigs-1) {
    if (sum(pconfigs) >= 1) {
      stop("pconfigs should sum up to <1")
    }

    pconfigs = c(pconfigs, 1-sum(pconfigs))
  } else {
    if (sum(pconfigs) != 1) {
      stop("pconfigs should sum up to 1")
    }
  }

  if (corr <= -1 || corr >= 1) {
    stop("corr must lie between 0 and 1")
  }

  if (any(is.na(oddsratios))) {
    stop("oddsratios must be provided")
  }

  if (length(oddsratios) != ncovariates) {
    stop("oddsratios must have", ncovariates, "elements")
  }

  if (any(oddsratios <= 0)) {
    stop("oddsratios must be positive")
  }

  if (is.na(responseprob)) {
    stop("responseprob must be provided")
  }

  if (responseprob <= 0 || responseprob >= 1) {
    stop("responseprob must lie between 0 and 1")
  }

  if (alpha < 0.00001 || alpha >= 1) {
    stop("alpha must lie in [0.00001, 1)")
  }


  xmean = colSums(pconfigs*x)
  names(xmean) = NULL
  psi = log(oddsratios)
  psi0 = qlogis(responseprob) - sum(xmean*psi)
  psi0star = psi0 + psi[1]*xmean[1]
  psistar = psi
  psistar[1] = 0

  theta = psi0 + as.numeric(x %*% psi)
  thetastar = psi0star + as.numeric(x %*% psistar)


  tau = 2*sum(pconfigs*(plogis(theta)*(theta - thetastar) -
                          log(1 + exp(theta)) +
                          log(1 + exp(thetastar))))

  tau = tau*(1-corr^2)


  nu = 1
  b = qchisq(1-alpha, nu)

  if (is.na(n)) {
    n0 = (qchisq(1-alpha, nu) - nu + qnorm(1-beta)*sqrt(2*nu))/tau
    while (pchisq(b, nu, n0*tau) > beta) n0 <- 2*n0
    n = uniroot(function(n) pchisq(b, nu, n*tau) - beta, c(0.5*n0, n0))$root
  }

  if (rounding) {
    n = ceiling(n)
  }

  power = pchisq(b, nu, n*tau, lower.tail = FALSE)

  des = list(power = power, alpha = alpha, n = n,
             ncovariates = ncovariates, nconfigs = nconfigs,
             x = x, pconfigs = pconfigs, corr = corr,
             oddsratios = oddsratios,
             responseprob = responseprob, effectsize = tau,
             calculationTarget = ifelse(is.na(beta), "beta", "n"),
             rounding = rounding)

  attr(des, "class") = "designLogistic"

  des
}


#' @title Power and sample size for Cohen's kappa
#' @description Obtains the power given sample size or obtains the sample
#' size given power for Cohen's kappa.
#'
#' @param beta The type II error.
#' @param n The total sample size.
#' @param ncats The number of categories.
#' @param kappaH0 The kappa coefficient under the null hypothesis.
#' @param kappa The kappa coefficient under the alternative hypothesis.
#' @param p1 The marginal probabilities for the first rater.
#' @param p2 The marginal probabilities for the second rater.
#'   Defaults to be equal to the marginal probabilities for the first
#'   rater if not provided.
#' @param rounding Whether to round up sample size. Defaults to 1 for
#'   sample size rounding.
#' @param alpha The one-sided significance level. Defaults to 0.025.
#'
#' @details
#' The kappa coefficient is defined as
#' \deqn{\kappa = \frac{\pi_o - \pi_e}{1 - \pi_e},} where
#' \eqn{\pi_o = \sum_i \pi_{ii}} is the observed agreement, and
#' \eqn{\pi_e = \sum_i \pi_{i.} \pi_{.i}} is the expected agreement
#' by chance.
#'
#' By Fleiss et al. (1969), the variance of \eqn{\hat{\kappa}} is given by
#' \deqn{Var(\hat{\kappa}) = \frac{v_1}{n},} where
#' \deqn{v_1 = \frac{Q_1 + Q_2 - Q3 - Q4}{(1-\pi_e)^4},}
#' \deqn{Q_1 = \pi_o(1-\pi_e)^2,}
#' \deqn{Q_2 = (1-\pi_o)^2 \sum_i \sum_j \pi_{ij}(\pi_{i.} + \pi_{.j})^2,}
#' \deqn{Q_3 = 2(1-\pi_o)(1-\pi_e) \sum_i \pi_{ii}(\pi_{i.} + \pi_{.i}),}
#' \deqn{Q_4 = (\pi_o \pi_e - 2\pi_e + \pi_o)^2.}
#'
#' Given \eqn{\kappa} and marginals
#' \eqn{\{(\pi_{i.}, \pi_{.i}): i=1,\ldots,k\}},
#' we obtain \eqn{\pi_o}. The only unknowns are the double summation
#' in \eqn{Q_2} and the single summation in \eqn{Q_3}.
#'
#' We find the optimal configuration of cell probabilities that yield the
#' maximum variance of \eqn{\hat{\kappa}} by treating the problem as a
#' linear programming problem with constraints to match the given
#' marginal probabilities and the observed agreement and ensure that
#' the cell probabilities are nonnegative. This is an extension of
#' Flack et al. (1988) by allowing unequal marginal probabilities
#' of the two raters.
#'
#' We perform the optimization under both the null and alternative
#' hypotheses to obtain \eqn{\max Var(\hat{\kappa} | \kappa = \kappa_0)}
#' and \eqn{\max Var(\hat{\kappa} | \kappa = \kappa_1)} for a
#' single subject, and then calculate the sample size or power
#' according to the following equation:
#' \deqn{\sqrt{n} |\kappa - \kappa_0| = z_{1-\alpha}
#' \sqrt{\max Var(\hat{\kappa} | \kappa = \kappa_0)} +
#' z_{1-\beta} \sqrt{\max Var(\hat{\kappa} | \kappa = \kappa_1)}.}
#'
#' @return An S3 class \code{designAgreement} object with the
#' following components:
#'
#' * \code{power}: The power to reject the null hypothesis.
#'
#' * \code{alpha}: The one-sided significance level.
#'
#' * \code{n}: The total sample size.
#'
#' * \code{ncats}: The number of categories.
#'
#' * \code{kappaH0}: The kappa coefficient under the null hypothesis.
#'
#' * \code{kappa}: The kappa coefficient under the alternative hypothesis.
#'
#' * \code{p1}: The marginal probabilities for the first rater.
#'
#' * \code{p2}: The marginal probabilities for the second rater.
#'
#' * \code{piH0}: The cell probabilities that maximize the
#'   variance of estimated kappa under H0.
#'
#' * \code{pi}: The cell probabilities that maximize the
#'   variance of estimated kappa under H1.
#'
#' * \code{calculationTarget}: The calculation target, \code{beta} or
#'   \code{n}.
#'
#' * \code{rounding}: Whether to round up sample size.
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @references
#'
#' V. F. Flack, A. A. Afifi, and P. A. Lachenbruch.
#' Sample size determinations for the two rater kappa statistic.
#' Psychometrika 1988; 53:321-325.
#'
#' @examples
#'
#' (design1 <- getDesignAgreement(
#'   beta = 0.2, n = NA, ncats = 4, kappaH0 = 0.4, kappa = 0.6,
#'   p1 = c(0.1, 0.2, 0.3, 0.4), p2 = c(0.15, 0.2, 0.24, 0.41),
#'   rounding = TRUE, alpha = 0.05))
#'
#' @export
#'
getDesignAgreement <- function(
    beta = NA_real_,
    n = NA_real_,
    ncats = NA_integer_,
    kappaH0 = NA_real_,
    kappa = NA_real_,
    p1 = NA_real_,
    p2 = NA_real_,
    rounding = TRUE,
    alpha = 0.025) {

  if (is.na(beta) && is.na(n)) {
    stop("beta and n cannot be both missing")
  }

  if (!is.na(beta) && !is.na(n)) {
    stop("Only one of beta and n should be provided")
  }

  if (!is.na(beta) && (beta >= 1-alpha || beta < 0.0001)) {
    stop("beta must lie in [0.0001, 1-alpha)");
  }

  if (!is.na(n) && (n <= 0 || n != round(n))) {
    stop("n must be a positive integer")
  }

  if (is.na(ncats)) {
    stop("ncats must be provided")
  }

  if (ncats < 2 || ncats != round(ncats)) {
    stop("ncats must be a positive integer >= 2")
  }

  if (is.na(kappaH0)) {
    stop("kappaH0 must be provided")
  }

  if (kappaH0 <= -1 || kappaH0 >= 1) {
    stop("kappaH0 must lie between -1 and 1")
  }

  if (is.na(kappa)) {
    stop("kappa must be provided")
  }

  if (kappa <= -1 || kappa >= 1) {
    stop("kappa must lie between -1 and 1")
  }

  if (any(is.na(p1))) {
    stop("p1 must be provided")
  }

  if (any(p1 <= 0 | p1 >= 1)) {
    stop(paste("Elements of p1 must lie between 0 and 1"))
  }

  if (length(p1) != ncats-1 && length(p1) != ncats) {
    stop(paste("p1 must have length", ncats-1, "or", ncats))
  }

  if (length(p1) == ncats-1) {
    if (sum(p1) >= 1) {
      stop("p1 must sum up to < 1 if it has only", ncats-1, "elements")
    }

    p1 = c(p1, 1 - sum(p1))
  } else {
    if (sum(p1) != 1) {
      stop(paste("p1 must sum to 1 across", ncats, "categories"))
    }
  }

  if (any(is.na(p2))) {
    p2 = p1
  }

  if (any(p2 <= 0 | p2 >= 1)) {
    stop(paste("Elements of p2 must lie between 0 and 1"))
  }

  if (length(p2) != ncats-1 && length(p2) != ncats) {
    stop(paste("p2 must have length", ncats-1, "or", ncats))
  }

  if (length(p2) == ncats-1) {
    if (sum(p2) >= 1) {
      stop("p2 must sum up to < 1 if it has only", ncats-1, "elements")
    }

    p2 = c(p2, 1 - sum(p2))
  } else {
    if (sum(p2) != 1) {
      stop(paste("p2 must sum to 1 across", ncats, "categories"))
    }
  }

  if (alpha < 0.00001 || alpha >= 1) {
    stop("alpha must lie in [0.00001, 1)")
  }


  k = ncats

  # obtain the maximum variance for given kappa and marginal probabilities
  fpi <- function(kappa, p1, p2) {
    pe = sum(p1*p2)
    po = pe + (1-pe)*kappa

    # equation (3) in Flack et al (1988) as a function of cell probabilities
    f.obj = rep(0, k^2)
    for (i in 1:k) {
      for (j in 1:k) {
        t = (i-1)*k + j
        if (i != j) {
          f.obj[t] = (1-po)^2*(p1[i]+p2[j])^2
        } else {
          f.obj[t] = -(1-po)*(p1[i]+p2[i])*(2*(1-pe) - (1-po)*(p1[i]+p2[i]))
        }
      }
    }

    # row margin, column margin, and observed agreement
    f.con = matrix(0, 2*k+1+k^2, k^2)
    for (i in 1:k) {
      for (j in 1:k) {
        t = (i-1)*k + j
        f.con[i,t] = 1
        f.con[k+j,t] = 1
        if (i==j) f.con[2*k+1,t] = 1
      }
    }

    # nonnegative cell probabilities
    f.con[(2*k+2):(2*k+1+k^2),] = diag(k^2)

    f.dir = c(rep("==", 2*k+1), rep(">=", k^2))
    f.rhs = c(p1, p2, po, rep(0, k^2))

    # linear programming
    opt = lpSolve::lp("max", f.obj, f.con, f.dir, f.rhs)
    matrix(opt$solution, k, k, byrow = TRUE)
  }


  # variance of kappa estimate for given cell probabilities
  fv1 <- function(pi) {
    p1 = rowSums(pi)
    p2 = colSums(pi)
    pd = diag(pi)
    po = sum(pd)
    pe = sum(p1*p2)

    mp1 = matrix(p1, k, k, byrow = TRUE)
    mp2 = matrix(p2, k, k)

    q1 = po*(1-pe)^2
    q2 = (1-po)^2*sum(pi*(mp1+mp2)^2)
    q3 = 2*(1-po)*(1-pe)*sum(pd*(p1+p2))
    q4 = (po*pe-2*pe+po)^2
    v1 = (q1+q2-q3-q4)/(1-pe)^4
    v1
  }


  piH0 = fpi(kappaH0, p1, p2)
  v1H0 = fv1(piH0)
  pi = fpi(kappa, p1, p2)
  v1 = fv1(pi)

  theta = abs(kappa - kappaH0)

  if (is.na(n)) {
    n = (qnorm(1-alpha)*sqrt(v1H0) + qnorm(1-beta)*sqrt(v1))^2/theta^2
  }

  if (rounding) {
    n = ceiling(n)
  }

  power = pnorm((theta*sqrt(n) - qnorm(1-alpha)*sqrt(v1H0))/sqrt(v1))

  des = list(power = power, alpha = alpha, n = n,
             ncats = ncats, kappaH0 = kappaH0, kappa = kappa,
             p1 = p1, p2 = p2, piH0 = piH0, pi = pi,
             calculationTarget = ifelse(is.na(beta), "beta", "n"),
             rounding = rounding)

  attr(des, "class") = "designAgreement"

  des
}


#' @title Power and sample size for one-sample Poisson rate exact test
#' @description Obtains the power given sample size or obtains the sample
#' size given power for one-sample Poisson rate.
#'
#' @param beta The type II error.
#' @param n The total sample size.
#' @param lambdaH0 The Poisson rate under the null hypothesis.
#' @param lambda The Poisson rate under the alternative hypothesis.
#' @param D The average exposure per subject.
#' @param alpha The one-sided significance level. Defaults to 0.025.
#'
#' @return A data frame containing the following variables:
#'
#' * \code{alpha}: The specified significance level.
#'
#' * \code{attainedAlpha}: The attained type I error of the exact test.
#'
#' * \code{power}: The actual power of the exact test.
#'
#' * \code{n}: The sample size.
#'
#' * \code{lambdaH0}: The Poisson rate under the null hypothesis.
#'
#' * \code{lambda}: The Poisson rate under the alternative hypothesis.
#'
#' * \code{D}: The average exposure per subject.
#'
#' * \code{r}: The critical value of the number of events for rejecting
#'   the null hypothesis. Reject H0 if \code{Y >= r} for upper-tailed
#'   test, and reject H0 if \code{Y <= r} for lower-tailed test.
#'
#' * \code{calculationTarget}: The calculation target, \code{beta} or
#'   \code{n}.
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @examples
#'
#' # Example 1: power calculation
#' (design1 <- getDesignOneRateExact(
#'   n = 525, lambdaH0 = 0.049, lambda = 0.012,
#'   D = 0.5, alpha = 0.025))
#'
#' # Example 2: sample size calculation
#' (design2 <- getDesignOneRateExact(
#'   beta = 0.2, lambdaH0 = 0.2, lambda = 0.3,
#'   D = 1, alpha = 0.05))
#'
#' @export
getDesignOneRateExact <- function(
    beta = NA_real_,
    n = NA_real_,
    lambdaH0 = NA_real_,
    lambda = NA_real_,
    D = 1,
    alpha = 0.025) {

  if (is.na(beta) && is.na(n)) {
    stop("beta and n cannot be both missing")
  }

  if (!is.na(beta) && !is.na(n)) {
    stop("Only one of beta and n should be provided")
  }

  if (!is.na(beta) && (beta >= 1-alpha || beta < 0.0001)) {
    stop("beta must lie in [0.0001, 1-alpha)");
  }

  if (!is.na(n) && (n <= 0 || n != round(n))) {
    stop("n must be a positive integer")
  }

  if (lambdaH0 <= 0) {
    stop("lambdaH0 must be positive");
  }

  if (lambda <= 0) {
    stop("lambda must be positive");
  }

  if (D <= 0) {
    stop("D must be positive")
  }

  if (alpha < 0.00001 || alpha >= 1) {
    stop("alpha must lie in [0.00001, 1)")
  }


  if (!is.na(n)) { # power calculation
    des = powerOneRateExact(n, lambdaH0, lambda, D, alpha)
  } else { # sample size calculation
    des = samplesizeOneRateExact(beta, lambdaH0, lambda, D, alpha)
  }

  des$calculationTarget = ifelse(is.na(beta), "beta", "n")

  des
}


#' @title Power and sample size for exact unconditional test for risk
#' difference
#' @description Obtains the power given sample size or obtains the sample
#' size given power for exact unconditional test of risk difference.
#'
#' @param beta The type II error.
#' @param n The total sample size.
#' @param riskDiffH0 The risk difference under the null hypothesis.
#'   Defaults to 0.
#' @param pi1 The assumed probability for the active treatment group.
#' @param pi2 The assumed probability for the control group.
#' @param allocationRatioPlanned Allocation ratio for the active treatment
#'   versus control. Defaults to 1 for equal randomization.
#' @param alpha The one-sided significance level. Defaults to 0.025.
#'
#' @return A data frame with the following variables:
#'
#' * \code{alpha}: The specified one-sided significance level.
#'
#' * \code{attainedAlpha}: The attained one-sided significance level.
#'
#' * \code{power}: The power.
#'
#' * \code{n}: The sample size.
#'
#' * \code{riskDiffH0}: The risk difference under the null hypothesis.
#'
#' * \code{pi1}: The assumed probability for the active treatment group.
#'
#' * \code{pi2}: The assumed probability for the control group.
#'
#' * \code{allocationRatioPlanned}: Allocation ratio for the active
#'   treatment versus control.
#'
#' * \code{zstatRiskDiffBound}: The critical value on the scale of
#'   score test statistic for risk difference.
#'
#' * \code{pi2star}: The response probability in the control group
#'   at which the critical value of the test statistic is attained.
#'
#' * \code{calculationTarget}: The calculation target, \code{beta} or
#'   \code{n}.
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @examples
#'
#' # Superiority test
#'
#' getDesignRiskDiffExact(n = 50, pi1 = 0.6, pi2 = 0.25, alpha = 0.025)
#'
#'
#' # Non-inferiority test
#'
#' getDesignRiskDiffExact(beta = 0.1, riskDiffH0 = -0.2,
#'                        pi1 = 0.8, pi2 = 0.8, alpha = 0.025)
#'
#'
#' @export
getDesignRiskDiffExact <- function(
    beta = NA_real_,
    n = NA_real_,
    riskDiffH0 = 0,
    pi1 = NA_real_,
    pi2 = NA_real_,
    allocationRatioPlanned = 1,
    alpha = 0.025) {

  if (is.na(beta) && is.na(n)) {
    stop("beta and n cannot be both missing")
  }

  if (!is.na(beta) && !is.na(n)) {
    stop("Only one of beta and n should be provided")
  }

  if (!is.na(beta) && (beta >= 1-alpha || beta < 0.0001)) {
    stop("beta must lie in [0.0001, 1-alpha)");
  }

  if (!is.na(n) && (n <= 0 || n != round(n))) {
    stop("n must be a positive integer")
  }

  if (riskDiffH0 <= -1 || riskDiffH0 >= 1) {
    stop("riskDiffH0 must lie between -1 and 1")
  }

  if (pi1 <= 0 || pi1 >= 1) {
    stop("pi1 must lie between 0 and 1")
  }

  if (pi2 <= 0 || pi2 >= 1) {
    stop("pi2 must lie between 0 and 1")
  }

  if (allocationRatioPlanned <= 0) {
    stop("allocationRatioPlanned must be positive")
  }

  if (alpha < 0.00001 || alpha >= 1) {
    stop("alpha must lie in [0.00001, 1)")
  }


  if (!is.na(n)) { # power calculation
    a = powerRiskDiffExact(n, riskDiffH0, pi1, pi2,
                           allocationRatioPlanned, alpha)
  } else { # sample size calculation
    a = samplesizeRiskDiffExact(beta, riskDiffH0, pi1, pi2,
                                allocationRatioPlanned, alpha)
  }

  a$calculationTarget = ifelse(is.na(beta), "beta", "n")

  a
}


#' @title Power and sample size for exact unconditional test for risk
#' ratio
#' @description Obtains the power given sample size or obtains the sample
#' size given power for exact unconditional test of risk ratio.
#'
#' @param beta The type II error.
#' @param n The total sample size.
#' @param riskRatioH0 The risk ratio under the null hypothesis.
#'   Defaults to 0.
#' @param pi1 The assumed probability for the active treatment group.
#' @param pi2 The assumed probability for the control group.
#' @param allocationRatioPlanned Allocation ratio for the active treatment
#'   versus control. Defaults to 1 for equal randomization.
#' @param alpha The one-sided significance level. Defaults to 0.025.
#'
#' @return A data frame with the following variables:
#'
#' * \code{alpha}: The specified one-sided significance level.
#'
#' * \code{attainedAlpha}: The attained one-sided significance level.
#'
#' * \code{power}: The power.
#'
#' * \code{n}: The sample size.
#'
#' * \code{riskRatioH0}: The risk ratio under the null hypothesis.
#'
#' * \code{pi1}: The assumed probability for the active treatment group.
#'
#' * \code{pi2}: The assumed probability for the control group.
#'
#' * \code{allocationRatioPlanned}: Allocation ratio for the active
#'   treatment versus control.
#'
#' * \code{zstatRiskRatioBound}: The critical value on the scale of
#'   score test statistic for risk ratio.
#'
#' * \code{pi2star}: The response probability in the control group
#'   at which the critical value of the test statistic is attained.
#'
#' * \code{calculationTarget}: The calculation target, \code{beta} or
#'   \code{n}.
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @examples
#'
#' # Non-inferiority test
#'
#' getDesignRiskRatioExact(beta = 0.1, riskRatioH0 = 0.833,
#'                         pi1 = 0.9, pi2 = 0.9, alpha = 0.025)
#'
#'
#' @export
getDesignRiskRatioExact <- function(
    beta = NA_real_,
    n = NA_real_,
    riskRatioH0 = 1,
    pi1 = NA_real_,
    pi2 = NA_real_,
    allocationRatioPlanned = 1,
    alpha = 0.025) {

  if (is.na(beta) && is.na(n)) {
    stop("beta and n cannot be both missing")
  }

  if (!is.na(beta) && !is.na(n)) {
    stop("Only one of beta and n should be provided")
  }

  if (!is.na(beta) && (beta >= 1-alpha || beta < 0.0001)) {
    stop("beta must lie in [0.0001, 1-alpha)");
  }

  if (!is.na(n) && (n <= 0 || n != round(n))) {
    stop("n must be a positive integer")
  }

  if (riskRatioH0 <= 0) {
    stop("riskRatioH0 must be positive")
  }

  if (pi1 <= 0 || pi1 >= 1) {
    stop("pi1 must lie between 0 and 1")
  }

  if (pi2 <= 0 || pi2 >= 1) {
    stop("pi2 must lie between 0 and 1")
  }

  if (allocationRatioPlanned <= 0) {
    stop("allocationRatioPlanned must be positive")
  }

  if (alpha < 0.00001 || alpha >= 1) {
    stop("alpha must lie in [0.00001, 1)")
  }


  if (!is.na(n)) { # power calculation
    a = powerRiskRatioExact(n, riskRatioH0, pi1, pi2,
                            allocationRatioPlanned, alpha)
  } else { # sample size calculation
    a = samplesizeRiskRatioExact(beta, riskRatioH0, pi1, pi2,
                                 allocationRatioPlanned, alpha)
  }

  a$calculationTarget = ifelse(is.na(beta), "beta", "n")

  a
}


#' @title Power and sample size for exact unconditional test for
#' equivalence in risk difference
#' @description Obtains the power given sample size or obtains the sample
#' size given power for exact unconditional test of equivalence in risk
#' difference.
#'
#' @param beta The type II error.
#' @param n The total sample size.
#' @param riskDiffLower The lower equivalence limit of risk difference.
#' @param riskDiffUpper The upper equivalence limit of risk difference.
#' @param pi1 The assumed probability for the active treatment group.
#' @param pi2 The assumed probability for the control group.
#' @param allocationRatioPlanned Allocation ratio for the active treatment
#'   versus control. Defaults to 1 for equal randomization.
#' @param alpha The significance level for each of the two one-sided
#'   tests. Defaults to 0.05.
#' @return A data frame with the following variables:
#'
#' * \code{alpha}: The specified significance level for each of the two
#'   one-sided tests.
#'
#' * \code{attainedAlpha}: The attained significance level.
#'
#' * \code{power}: The power.
#'
#' * \code{n}: The sample size.
#'
#' * \code{riskDiffLower}: The lower equivalence limit of risk difference.
#'
#' * \code{riskDiffUpper}: The upper equivalence limit of risk difference.
#'
#' * \code{pi1}: The assumed probability for the active treatment group.
#'
#' * \code{pi2}: The assumed probability for the control group.
#'
#' * \code{allocationRatioPlanned}: Allocation ratio for the active
#'   treatment versus control.
#'
#' * \code{zstatRiskDiffLower}: The efficacy boundaries on the
#'   z-test statistic scale for the one-sided null hypothesis on the
#'   lower equivalence limit.
#'
#' * \code{zstatRiskDiffUpper}: The efficacy boundaries on the
#'   z-test statistic scale for the one-sided null hypothesis on the
#'   upper equivalence limit.
#'
#' * \code{calculationTarget}: The calculation target, \code{beta} or
#'   \code{n}.
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @examples
#'
#' getDesignRiskDiffExactEquiv(
#'   n = 200, riskDiffLower = -0.2, riskDiffUpper = 0.2,
#'   pi1 = 0.775, pi2 = 0.775, alpha = 0.05)
#'
#' @export
getDesignRiskDiffExactEquiv <- function(
    beta = NA_real_,
    n = NA_real_,
    riskDiffLower = NA_real_,
    riskDiffUpper = NA_real_,
    pi1 = NA_real_,
    pi2 = NA_real_,
    allocationRatioPlanned = 1,
    alpha = 0.05) {

  if (is.na(beta) && is.na(n)) {
    stop("beta and n cannot be both missing")
  }

  if (!is.na(beta) && !is.na(n)) {
    stop("Only one of beta and n should be provided")
  }

  if (!is.na(beta) && (beta >= 1-alpha || beta < 0.0001)) {
    stop("beta must lie in [0.0001, 1-alpha)");
  }

  if (!is.na(n) && (n <= 0 || n != round(n))) {
    stop("n must be a positive integer")
  }

  if (is.na(riskDiffLower)) {
    stop("riskDiffLower must be provided")
  }

  if (is.na(riskDiffUpper)) {
    stop("riskDiffUpper must be provided")
  }

  if (riskDiffLower <= -1) {
    stop("riskDiffLower must be greater than -1")
  }

  if (riskDiffUpper >= 1) {
    stop("riskDiffUpper must be less than 1")
  }

  if (riskDiffLower >= riskDiffUpper) {
    stop("riskDiffLower must be less than riskDiffUpper")
  }

  if (pi1 <= 0 || pi1 >= 1) {
    stop("pi1 must lie between 0 and 1")
  }

  if (pi2 <= 0 || pi2 >= 1) {
    stop("pi2 must lie between 0 and 1")
  }

  if (allocationRatioPlanned <= 0) {
    stop("allocationRatioPlanned must be positive")
  }

  if (alpha < 0.00001 || alpha >= 1) {
    stop("alpha must lie in [0.00001, 1)")
  }


  if (!is.na(n)) { # power calculation
    a = powerRiskDiffExactEquiv(n, riskDiffLower, riskDiffUpper,
                                pi1, pi2, allocationRatioPlanned, alpha)
  } else { # sample size calculation
    a = samplesizeRiskDiffExactEquiv(beta, riskDiffLower, riskDiffUpper,
                                     pi1, pi2, allocationRatioPlanned, alpha)
  }

  a$calculationTarget = ifelse(is.na(beta), "beta", "n")

  a
}


#' @title Power and sample size for exact unconditional test for
#' equivalence in risk ratio
#' @description Obtains the power given sample size or obtains the sample
#' size given power for exact unconditional test of equivalence in risk
#' ratio.
#'
#' @param beta The type II error.
#' @param n The total sample size.
#' @param riskRatioLower The lower equivalence limit of risk ratio.
#' @param riskRatioUpper The upper equivalence limit of risk ratio.
#' @param pi1 The assumed probability for the active treatment group.
#' @param pi2 The assumed probability for the control group.
#' @param allocationRatioPlanned Allocation ratio for the active treatment
#'   versus control. Defaults to 1 for equal randomization.
#' @param alpha The significance level for each of the two one-sided
#'   tests. Defaults to 0.05.
#' @return A data frame with the following variables:
#'
#' * \code{alpha}: The specified significance level for each of the two
#'   one-sided tests.
#'
#' * \code{attainedAlpha}: The attained significance level.
#'
#' * \code{power}: The power.
#'
#' * \code{n}: The sample size.
#'
#' * \code{riskRatioLower}: The lower equivalence limit of risk ratio.
#'
#' * \code{riskRatioUpper}: The upper equivalence limit of risk ratio.
#'
#' * \code{pi1}: The assumed probability for the active treatment group.
#'
#' * \code{pi2}: The assumed probability for the control group.
#'
#' * \code{allocationRatioPlanned}: Allocation ratio for the active
#'   treatment versus control.
#'
#' * \code{zstatRiskRatioLower}: The efficacy boundaries on the
#'   z-test statistic scale for the one-sided null hypothesis on the
#'   lower equivalence limit.
#'
#' * \code{zstatRiskRatioUpper}: The efficacy boundaries on the
#'   z-test statistic scale for the one-sided null hypothesis on the
#'   upper equivalence limit.
#'
#' * \code{calculationTarget}: The calculation target, \code{beta} or
#'   \code{n}.
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @examples
#'
#' getDesignRiskRatioExactEquiv(
#'   n = 200, riskRatioLower = 0.8, riskRatioUpper = 1.25,
#'   pi1 = 0.775, pi2 = 0.775, alpha = 0.05)
#'
#' @export
getDesignRiskRatioExactEquiv <- function(
    beta = NA_real_,
    n = NA_real_,
    riskRatioLower = NA_real_,
    riskRatioUpper = NA_real_,
    pi1 = NA_real_,
    pi2 = NA_real_,
    allocationRatioPlanned = 1,
    alpha = 0.05) {

  if (is.na(beta) && is.na(n)) {
    stop("beta and n cannot be both missing")
  }

  if (!is.na(beta) && !is.na(n)) {
    stop("Only one of beta and n should be provided")
  }

  if (!is.na(beta) && (beta >= 1-alpha || beta < 0.0001)) {
    stop("beta must lie in [0.0001, 1-alpha)");
  }

  if (!is.na(n) && (n <= 0 || n != round(n))) {
    stop("n must be a positive integer")
  }

  if (is.na(riskRatioLower)) {
    stop("riskRatioLower must be provided")
  }

  if (is.na(riskRatioUpper)) {
    stop("riskRatioUpper must be provided")
  }

  if (riskRatioLower <= 0) {
    stop("riskRatioLower must be positive")
  }

  if (riskRatioLower >= riskRatioUpper) {
    stop("riskRatioLower must be less than riskRatioUpper")
  }

  if (pi1 <= 0 || pi1 >= 1) {
    stop("pi1 must lie between 0 and 1")
  }

  if (pi2 <= 0 || pi2 >= 1) {
    stop("pi2 must lie between 0 and 1")
  }

  if (allocationRatioPlanned <= 0) {
    stop("allocationRatioPlanned must be positive")
  }

  if (alpha < 0.00001 || alpha >= 1) {
    stop("alpha must lie in [0.00001, 1)")
  }


  if (!is.na(n)) { # power calculation
    a = powerRiskRatioExactEquiv(n, riskRatioLower, riskRatioUpper,
                                 pi1, pi2, allocationRatioPlanned, alpha)
  } else { # sample size calculation
    a = samplesizeRiskRatioExactEquiv(beta, riskRatioLower, riskRatioUpper,
                                      pi1, pi2, allocationRatioPlanned,
                                      alpha)
  }

  a$calculationTarget = ifelse(is.na(beta), "beta", "n")

  a
}


