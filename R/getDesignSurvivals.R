#' @title The Quantiles of a Survival Distribution
#' @description Obtains the quantiles of a survival distribution.
#'
#' @param S The survival function of a univariate survival time.
#' @param probs The numeric vector of probabilities.
#' @param ... Additional arguments to be passed to S.
#'
#' @return A vector of \code{length(probs)} for the quantiles.
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @examples
#'
#' fquantile(pweibull, probs = c(0.25, 0.5, 0.75),
#'           shape = 1.37, scale = 1/0.818, lower.tail = FALSE)
#'
#' @export
fquantile <- function(S, probs, ...) {
  Sinf = S(Inf, ...)

  if (any(probs > 1 - Sinf)) {
    stop(paste("probs must be less than or equal to", 1 - Sinf))
  }

  sapply(probs, function(p) {
    q = 1 - p

    lower = 0
    upper = 1
    while (S(upper, ...) > q) {
      lower = upper
      upper = 2*upper
    }

    uniroot(f = function(t) S(t, ...) - q,
            interval = c(lower, upper))$root

  })
}


#' @title Profile Log-Likelihood Function for the Change Points in
#' Piecewise Exponential Approximation
#' @description Obtains the profile log-likelihood function for the
#' change points in the piecewise exponential approximation to
#' a survival function.
#'
#' @param tau The numeric vector of change points.
#' @param S The survival function of a univariate survival time.
#' @param ... Additional arguments to be passed to S.
#'
#' @return A list with the following three components:
#'
#' * \code{piecewiseSurvivalTime}: A vector that specifies the starting
#'   time of piecewise exponential survival time intervals.
#'
#' * \code{lambda}: A vector of hazard rates for the event. One for
#'   each analysis time interval.
#'
#' * \code{loglik}: The value of the profile log-likelihood.
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @examples
#'
#' pwexploglik(tau = c(0.5, 1.2, 2.8), pweibull,
#'             shape = 1.37, scale = 1/0.818, lower.tail = FALSE)
#'
#' @export
pwexploglik <- function(tau, S, ...) {
  J = length(tau) + 1
  t = c(0, tau, Inf)

  surv = S(t, ...)
  d = -diff(surv)

  ex = rep(0,J)
  for (j in 1:J) {
    ex[j] = integrate(f = function(x) S(x, ...), t[j], t[j+1])$value
  }

  lambda = d/ex
  loglik = sum(d*log(lambda)) - 1

  list(piecewiseSurvivalTime = t[1:J], lambda = lambda, loglik = loglik)
}


#' @title Piecewise Exponential Approximation to a Survival Distribution
#' @description Obtains the piecewise exponential distribution that
#' approximates a survival distribution.
#'
#' @param S The survival function of a univariate survival time.
#' @param ... Additional arguments to be passed to S.
#'
#' @return A list with three components:
#'
#' * \code{piecewiseSurvivalTime}: A vector that specifies the starting
#'   time of piecewise exponential survival time intervals.
#'   Must start with 0, e.g., c(0, 6) breaks the time axis into 2 event
#'   intervals: [0, 6) and [6, Inf).
#'
#' * \code{lambda}: A vector of hazard rates for the event. One for
#'   each analysis time interval.
#'
#' * \code{loglik}: The sequence of the asymptotic limit of the
#'   piecewise exponential log-likelihood for an increasing number
#'   of change points.
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @examples
#'
#' # Example 1: Piecewise exponential
#' pwexpcuts(ptpwexp, piecewiseSurvivalTime = c(0, 3.4, 5.5),
#'           lambda = c(0.0168, 0.0833, 0.0431), lowerBound = 0,
#'           lower.tail = FALSE)
#'
#' # Example 2: Weibull
#' pwexpcuts(pweibull, shape = 1.37, scale = 1/0.818, lower.tail = FALSE)
#'
#' @export
pwexpcuts <- function(S, ...) {
  Sinf = S(Inf, ...)
  tmax = fquantile(S, 0.999*(1 - Sinf), ...)

  Jmax = 100
  tau = rep(0, Jmax-1)
  loglik = rep(0, Jmax-1)
  for (J in 2:Jmax) {
    if (J == 2) {
      eps = tmax*1.0e-6
      lower = eps
      upper = tmax - eps

      x = optimize(f = function(x) pwexploglik(x, S, ...)$loglik,
                   interval = c(lower, upper), maximum = TRUE)$maximum

      out = pwexploglik(x, S, ...)
      tau[J-1] = out$piecewiseSurvivalTime[J]
      loglik[J-1] = out$loglik
    } else {
      eps = (tmax - tau[J-2])*1.0e-6
      lower = tau[J-2] + eps
      upper = tmax - eps

      x = optimize(f = function(x) {
        pwexploglik(c(tau[1:(J-2)], x), S, ...)$loglik
      }, interval = c(lower, upper), maximum = TRUE)$maximum

      out = pwexploglik(c(tau[1:(J-2)], x), S, ...)
      tau[J-1] = out$piecewiseSurvivalTime[J]
      loglik[J-1] = out$loglik

      if (abs((loglik[J-1] - loglik[J-2])/loglik[J-2]) < 0.0001) break
    }
  }

  list(piecewiseSurvivalTime = out$piecewiseSurvivalTime,
       lambda = out$lambda, loglik = loglik[1:(J-1)])
}


#' @title Schoenfeld Method for Log-Rank Test Sample Size Calculation
#' @description Obtains the sample size and study duration by calibrating
#' the number of events calculated using the Schoenfeld formula
#' under the proportional hazards assumption.
#'
#' @param beta Type II error. Defaults to 0.2.
#' @inheritParams param_kMax
#' @param informationRates The information rates in terms of number
#'   of events for the conventional log-rank test and in terms of
#'   the actual information for weighted log-rank tests.
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
#' @inheritParams param_hazardRatioH0
#' @inheritParams param_allocationRatioPlanned
#' @inheritParams param_accrualTime
#' @inheritParams param_accrualIntensity
#' @inheritParams param_piecewiseSurvivalTime
#' @inheritParams param_stratumFraction
#' @param hazardRatio Hazard ratio under the alternative hypothesis for
#'   the active treatment versus control.
#' @inheritParams param_lambda2_stratified
#' @inheritParams param_gamma1_stratified
#' @inheritParams param_gamma2_stratified
#' @inheritParams param_followupTime
#' @inheritParams param_fixedFollowup
#' @param interval The interval to search for the solution of
#'   followupTime. Defaults to \code{c(0.001, 240)}.
#' @param spendingTime A vector of length \code{kMax} for the error spending
#'   time at each analysis. Defaults to missing, in which case, it is the
#'   same as \code{informationRates}.
#' @param rounding Whether to round up sample size and events.
#'   Defaults to 1 for sample size rounding.
#' @param calibrate Whether to use simulations to calibrate the number of
#'   events calculated using the Schoenfeld formula.
#' @param maxNumberOfIterations The number of simulation iterations.
#'   Defaults to 10000.
#' @param maxNumberOfRawDatasetsPerStage The number of raw datasets per
#'   stage to extract.
#' @param seed The seed to reproduce the simulation results.
#'   The seed from the environment will be used if left unspecified.
#'
#' @return A list of two components:
#'
#' * \code{resultsUnderH1}: An S3 class \code{lrpower} object under the
#' alternative hypothesis.
#'
#' * \code{resultsUnderH0}: An S3 class \code{lrpower} object under the
#' null hypothesis.
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @examples
#'
#' (lr1 <- lrschoenfeld(
#'   beta = 0.1, kMax = 2, alpha = 0.025,
#'   hazardRatioH0 = 1, allocationRatioPlanned = 1,
#'   accrualIntensity = 20, hazardRatio = 0.3,
#'   lambda2 = 1.9/12,
#'   gamma1 = -log(1-0.1)/24, gamma2 = -log(1-0.1)/24,
#'   fixedFollowup = 0, rounding = 1,
#'   calibrate = 0, maxNumberOfIterations = 1000,
#'   seed = 12345))
#'
#' (lr2 <- lrschoenfeld(
#'   beta = 0.1, kMax = 2, alpha = 0.025,
#'   hazardRatioH0 = 1, allocationRatioPlanned = 1,
#'   accrualIntensity = 20, hazardRatio = 0.3,
#'   lambda2 = 1.9/12,
#'   gamma1 = -log(1-0.1)/24, gamma2 = -log(1-0.1)/24,
#'   fixedFollowup = 0, rounding = 1,
#'   calibrate = 1, maxNumberOfIterations = 1000,
#'   seed = 12345))
#'
#' @export
#'
lrschoenfeld <- function(
    beta = 0.2,
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
    hazardRatioH0 = 1,
    allocationRatioPlanned = 1,
    accrualTime = 0L,
    accrualIntensity = NA_real_,
    piecewiseSurvivalTime = 0L,
    stratumFraction = 1L,
    hazardRatio = NA_real_,
    lambda2 = NA_real_,
    gamma1 = 0L,
    gamma2 = 0L,
    followupTime = NA_real_,
    fixedFollowup = 0L,
    interval = as.numeric(c(0.001, 240)),
    spendingTime = NA_real_,
    rounding = 1L,
    calibrate = 1L,
    maxNumberOfIterations = 10000L,
    maxNumberOfRawDatasetsPerStage = 0L,
    seed = NA_integer_) {

  lambda1 <- lambda2*hazardRatio

  if (!fixedFollowup) { # obtain accrual duration and then find followup time
    d <- getNeventsFromHazardRatio(
      beta, kMax, informationRates,
      efficacyStopping, futilityStopping,
      criticalValues, alpha, typeAlphaSpending,
      parameterAlphaSpending, userAlphaSpending,
      futilityBounds, typeBetaSpending,
      parameterBetaSpending, userBetaSpending,
      spendingTime, hazardRatioH0, hazardRatio,
      allocationRatioPlanned, rounding)

    durations <- getDurationFromNevents(
      d, allocationRatioPlanned,
      accrualTime, accrualIntensity,
      piecewiseSurvivalTime, stratumFraction,
      lambda1, lambda2, gamma1, gamma2,
      followupTime, fixedFollowup, 2)

    nsubjects <- ceiling(mean(durations$subjects))

    accrualDuration <- getAccrualDurationFromN(
      nsubjects, accrualTime, accrualIntensity)

    lrc1 <- lrsamplesize(
      beta, kMax, informationRates,
      efficacyStopping, futilityStopping,
      criticalValues, alpha, typeAlphaSpending,
      parameterAlphaSpending, userAlphaSpending,
      futilityBounds, typeBetaSpending,
      parameterBetaSpending, userBetaSpending,
      hazardRatioH0, allocationRatioPlanned,
      accrualTime, accrualIntensity,
      piecewiseSurvivalTime, stratumFraction,
      lambda1, lambda2, gamma1, gamma2,
      accrualDuration, followupTime,
      fixedFollowup, 0, 0, 0, "schoenfeld",
      interval, spendingTime, rounding)

    lrp1 <- lrc1$resultsUnderH1
  } else { # look for accrual duration directly
    lrc1 <- lrsamplesize(
      beta, kMax, informationRates,
      efficacyStopping, futilityStopping,
      criticalValues, alpha, typeAlphaSpending,
      parameterAlphaSpending, userAlphaSpending,
      futilityBounds, typeBetaSpending,
      parameterBetaSpending, userBetaSpending,
      hazardRatioH0, allocationRatioPlanned,
      accrualTime, accrualIntensity,
      piecewiseSurvivalTime, stratumFraction,
      lambda1, lambda2, gamma1, gamma2,
      NA, followupTime,
      fixedFollowup, 0, 0, 0, "schoenfeld",
      interval, spendingTime, rounding)

    lrp1 <- lrc1$resultsUnderH1
    accrualDuration <- lrp1$overallResults$accrualDuration
  }

  # run a simulation to verify the analytic results
  d <- lrp1$overallResults$numberOfEvents
  plannedEvents <- floor(lrp1$byStageResults$numberOfEvents + 0.5)
  informationRates <- lrp1$byStageResults$informationRates
  allocation <- float_to_fraction(allocationRatioPlanned)

  lrs1 <- lrsim(
    kMax, informationRates,
    lrp1$byStageResults$efficacyBounds,
    lrp1$byStageResults$futilityBounds,
    hazardRatioH0, allocation[1], allocation[2],
    accrualTime, accrualIntensity,
    piecewiseSurvivalTime, stratumFraction,
    lambda1, lambda2, gamma1, gamma2,
    accrualDuration, followupTime,
    fixedFollowup, 0, 0,
    plannedEvents, NA,
    maxNumberOfIterations,
    maxNumberOfRawDatasetsPerStage, seed)

  if (calibrate) {
    p <- lrs1$overview$overallReject

    # calibrate number of events
    d <- d*((qnorm(1-alpha) + qnorm(1-beta))/(qnorm(1-alpha) + qnorm(p)))^2
    if (rounding) {
      d <- ceiling(d - 1.0e-12)
      nevents <- floor(d*informationRates + 0.5)
      informationRates <- nevents/d
    } else {
      nevents <- d*informationRates
    }

    if (!fixedFollowup) { # fix accrual duration and find followup time
      studyDuration <- caltime(
        d, allocationRatioPlanned,
        accrualTime, accrualIntensity,
        piecewiseSurvivalTime, stratumFraction,
        lambda1, lambda2, gamma1, gamma2,
        accrualDuration, 1000, fixedFollowup)

      followupTime <- studyDuration - accrualDuration
    } else { # update accrual duration directly
      durations <- getDurationFromNevents(
        d, allocationRatioPlanned,
        accrualTime, accrualIntensity,
        piecewiseSurvivalTime, stratumFraction,
        lambda1, lambda2, gamma1, gamma2,
        followupTime, fixedFollowup, 2)

      nsubjects <- ceiling(durations$subjects[1])

      accrualDuration <- getAccrualDurationFromN(
        nsubjects, accrualTime, accrualIntensity)

      studyDuration <- caltime(
        d, allocationRatioPlanned,
        accrualTime, accrualIntensity,
        piecewiseSurvivalTime, stratumFraction,
        lambda1, lambda2, gamma1, gamma2,
        accrualDuration, followupTime, fixedFollowup)
    }

    lrp1 <- lrpower(
      kMax, informationRates,
      efficacyStopping, futilityStopping,
      criticalValues, alpha, typeAlphaSpending,
      parameterAlphaSpending, userAlphaSpending,
      futilityBounds, typeBetaSpending,
      parameterBetaSpending,
      hazardRatioH0, allocationRatioPlanned,
      accrualTime, accrualIntensity,
      piecewiseSurvivalTime, stratumFraction,
      lambda1, lambda2, gamma1, gamma2,
      accrualDuration, followupTime,
      fixedFollowup, 0, 0, 0, "schoenfeld",
      spendingTime, studyDuration)

    lrs1 <- lrsim(
      kMax, lrp1$byStageResults$informationRates,
      lrp1$byStageResults$efficacyBounds,
      lrp1$byStageResults$futilityBounds,
      hazardRatioH0, allocation[1], allocation[2],
      accrualTime, accrualIntensity,
      piecewiseSurvivalTime, stratumFraction,
      lambda1, lambda2, gamma1, gamma2,
      accrualDuration, followupTime,
      fixedFollowup, 0, 0, nevents, NA,
      maxNumberOfIterations,
      maxNumberOfRawDatasetsPerStage, seed)
  }

  list(analyticalResults = lrp1,
       simulationResults = lrs1)
}
