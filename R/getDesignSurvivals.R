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


#' @title Profile Log-Likelihood Function for Change Points in
#' Piecewise Exponential Approximation
#' @description Obtains the profile log-likelihood function for
#' change points in the piecewise exponential approximation to
#' a survival function.
#'
#' @param tau The numeric vector of change points.
#' @param S The survival function of a univariate survival time.
#' @param ... Additional arguments to be passed to S.
#'
#' @details
#' This function computes the profile log-likelihood for
#' change points in a piecewise exponential survival model.
#'
#' Let \eqn{S(t)} denote the survival function of a univariate
#' survival time, and \eqn{\tau} be a vector of \eqn{J-1} change points.
#' The piecewise exponential survival model divides the time
#' axis into \eqn{J} intervals defined by the change points \eqn{\tau},
#' where each interval \eqn{[t_j, t_{j+1})} has a constant
#' hazard rate \eqn{\lambda_j}. The time intervals are specified as:
#' \deqn{[t_1, t_2), [t_2, t_3), \ldots, [t_{J}, t_{J+1})}
#' where \eqn{t_1 = 0}, \eqn{t_{J+1} = \infty}, and
#' \eqn{t_j = \tau_{j-1}} for \eqn{j = 2, \ldots, J}.
#'
#' For each subject, the expected number of events occurring
#' in the \eqn{j}-th interval is
#' \deqn{d_j = E\{I(t_j < Y \leq t_{j+1})\} = S(t_j) - S(t_{j+1})}
#' The expected exposure in the \eqn{j}-th interval is:
#' \deqn{e_j = E\{(Y-t_j)I(t_j < Y \leq t_{j+1}) +
#' (t_{j+1} - t_j)I(Y > t_{j+1})\}}
#' which can be shown to be equivalent to
#' \deqn{e_j = \int_{t_j}^{t_{j+1}} S(t) dt}
#'
#' The log-likelihood for the piecewise exponential model is:
#' \deqn{\ell(\tau,\lambda) =
#' \sum_{j=1}^J \{d_j \log(\lambda_j) - e_j \lambda_j\}}
#' The profile log-likelihood for \eqn{\tau} is obtained by maximizing
#' \eqn{\ell(\tau,\lambda)} with respect to
#' \eqn{\lambda} for fixed \eqn{\tau}.
#' The maximum likelihood estimate of the hazard rate in the
#' \eqn{j}-th interval is
#' \deqn{\lambda_j = \frac{d_j}{e_j}}
#' Substituting back, the profile log-likelihood is
#' \deqn{\ell(\tau) = \sum_{j=1}^J d_j \log(d_j/e_j) - 1}
#' where we use the fact that
#' \eqn{\sum_{j=1}^J d_j = 1}.
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
#' @param tol The tolerance for convergence of the profile log-likelihood.
#'   Defaults to 0.0001.
#'
#' @details
#' This function computes the piecewise exponential approximation
#' to a survival distribution.
#' The piecewise exponential model divides the time axis into
#' \eqn{J} intervals defined by the change points, where each
#' interval \eqn{[t_j, t_{j+1})} has a constant hazard rate
#' \eqn{\lambda_j}. The time intervals are specified as:
#' \deqn{[t_1, t_2), [t_2, t_3), \ldots, [t_{J}, t_{J+1})}
#' where \eqn{t_1 = 0}, \eqn{t_{J+1} = \infty}, and
#' \eqn{t_j = \tau_{j-1}} for \eqn{j = 2, \ldots, J}.
#' The function starts with \eqn{J = 2} (1 change point) and
#' gradually increases \eqn{J} by adding one change point at a time
#' until the maximized profile log-likelihood for change points
#' stabilizes, i.e., the relative increase in the maximum of the
#' profile log-likelihood function is less than \code{tol}.
#' If the relative change in the hazard rate is also less than
#' \code{tol}, the function stops and returns the results.
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
pwexpcuts <- function(S, ..., tol = 0.0001) {
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

      if (abs((loglik[J-1] - loglik[J-2])/loglik[J-2]) < tol) break
    }
  }

  if (abs((out$lambda[J] - out$lambda[J-1])/out$lambda[J-1]) < tol) {
    J = J - 1
  }

  list(piecewiseSurvivalTime = out$piecewiseSurvivalTime[1:J],
       lambda = out$lambda[1:J], loglik = loglik[1:(J-1)])
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
#' @details
#' This function calculates the sample size and study duration
#' by calibrating the number of events estimated using the
#' Schoenfeld formula under the proportional hazards assumption,
#' particularly when the hazard ratio is far away from one and/or
#' the allocation between groups is unequal.
#'
#' For a fixed design, the Schoenfeld formula for the required
#' number of events is
#' \deqn{D = \frac{(\Phi^{-1}(1-\alpha) + \Phi^{-1}(1-\beta))^2}
#' {(\theta - \theta_0)^2 r(1-r)}}
#' where \eqn{D} is the total number of events required,
#' \eqn{\alpha} is the type I error rate,
#' \eqn{\beta} is the type II error rate,
#' \eqn{r} is the randomization probability for the active treatment group,
#' \eqn{\theta_0} and \eqn{\theta} are the log hazard ratios under
#' the null and alternative hypotheses, respectively.
#'
#' The function first computes the number of events using the
#' Schoenfeld formula. If \code{calibrate} is set to 1, the
#' function uses simulations to calibrate the number of
#' events, accounting for scenarios where the Schoenfeld formula
#' may be inaccurate (e.g., when allocation is unequal or the hazard
#' ratio is extreme).
#'
#' Let \eqn{D_{schoenfeld}} be the number of events calculated
#' by the Schoenfeld formula, and \eqn{D_{calibrated}}
#' be the calibrated number of events. The calibrated number of
#' events is calculated as
#' #' \deqn{D_{\text{calibrated}} =
#' \frac{\left\{\Phi^{-1}(1-\alpha) + \Phi^{-1}(1-\beta)\right\}^2}
#' {\left\{\Phi^{-1}(1-\alpha) +
#' \Phi^{-1}(1-\beta_{\text{schoenfeld}})\right\}^2}
#' D_{\text{schoenfeld}}}
#' where \eqn{\beta_{schoenfeld}} is the empirical type II error
#' estimated via simulation.
#'
#' A second round of simulation is performed to obtain the
#' empirical power using the calibrated number of events.
#'
#' @return A list of two components:
#'
#' * \code{analyticalResults}: An S3 class \code{lrpower} object for
#'   the asymptotic power.
#'
#' * \code{simulationResults}: An S3 class \code{lrsim} object for
#'   the empirical power.
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
