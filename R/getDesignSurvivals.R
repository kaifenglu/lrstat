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


#' @title Brookmeyer-Crowley Confidence Interval for Quantiles of
#' Right-Censored Time-to-Event Data
#' @description Obtains the Brookmeyer-Crowley confidence
#' interval for quantiles of right-censored time-to-event data.
#'
#' @param time The vector of possibly right-censored survival times.
#' @param event The vector of event indicators.
#' @param cilevel The confidence interval level. Defaults to 0.95.
#' @param transform The transformation of the survival function to use
#'   to construct the confidence interval. Options include "linear",
#'   "loglog", "log", "asinsqrt", and "logit". Defaults to "loglog".
#' @param probs The vector of probabilities to calculate the quantiles.
#'   Defaults to c(0.25, 0.5, 0.75).
#'
#' @return A data frame containing the estimated quantile and
#' confidence interval corresponding to each specified probability.
#' It includes the following variables:
#'
#' * \code{prob}: The probability to calculate the quantile.
#'
#' * \code{quantile}: The estimated quantile.
#'
#' * \code{lower}: The lower limit of the confidence interval.
#'
#' * \code{upper}: The upper limit of the confidence interval.
#'
#' * \code{cilevel}: The confidence interval level.
#'
#' * \code{transform}: The transformation of the survival function to use
#'   to construct the confidence interval.
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @examples
#'
#' survQuantile(
#'   time = c(33.7, 3.9, 10.5, 5.4, 19.5, 23.8, 7.9, 16.9, 16.6,
#'            33.7, 17.1, 7.9, 10.5, 38),
#'   event = c(0, 1, 1, 1, 1, 0, 1, 0, 0, 0, 0, 0, 1, 1),
#'   probs = c(0.25, 0.5, 0.75))
#'
#' @export
#'
survQuantile <- function(
    time = NA_real_, event = NA_real_,
    cilevel = 0.95, transform = "loglog",
    probs = c(0.25, 0.5, 0.75)) {

  if (any(is.na(time))) {
    stop("time must be provided")
  }

  if (any(is.na(event))) {
    stop("event must be provided")
  }

  if (any(time <= 0)) {
    stop("time must be positive for each subject")
  }

  if (any(event != 1 & event != 0)) {
    stop("event must be 1 or 0 for each subject")
  }

  if (cilevel <= 0 || cilevel >= 1) {
    stop("cilevel must lie between 0 and 1")
  }

  transform = tolower(transform)
  if (!(transform %in% c("linear", "loglog", "log", "asinsqrt", "logit"))) {
    stop(paste("transform must be one of the options:",
               "linear, loglog, log, asinsqrt, or logit"))
  }

  if (any(probs <= 0 | probs >= 1)) {
    stop("Elements of probs must lie between 0 and 1")
  }


  data = data.frame(time = time, event = event)

  # sort the data by time with event appearing before censoring for ties
  df1 = data[order(data$time, -data$event),]
  n = nrow(df1)

  # construct the data for # at risk and # events at distinct event times
  df2 = data.frame()
  cache = 0 # buffer for the current event time
  for (i in 1:n) {
    if ((i == 1 && df1$event[i] == 1) ||
        (i >= 2 && df1$event[i] == 1 &&
         df1$time[i] > df1$time[i-1])) { # new event
      # add the info for the previous event
      if (cache) {
        df2 = rbind(df2, data.frame(time = t, nrisk = nrisk,
                                    nevent = nevent))
      }

      # update the buffer for the current event time
      t = df1$time[i]
      nrisk = n-i+1
      nevent = 1
      cache = 1
    } else if (i >= 2 && df1$event[i] == 1 &&
               df1$event[i-1] == 1 &&
               df1$time[i] == df1$time[i-1]) { # tied event
      nevent = nevent + 1
    } else if (i >= 2 && df1$event[i] == 0 && df1$event[i-1] == 1) {
      # new censoring
      # add the info for the previous event
      df2 = rbind(df2, data.frame(time = t, nrisk = nrisk,
                                  nevent = nevent))

      # empty the buffer for the current event time
      cache = 0
    }
  }

  # add the info for the last event
  if (cache) {
    df2 = rbind(df2, data.frame(time = t, nrisk = nrisk,
                                nevent = nevent))
  }


  # construct the Kaplan-Meier estimates of survival probabilities
  df2$surv = cumprod(1 - df2$nevent/df2$nrisk)

  # obtain the Greenwood variance estimate of survival probabilities
  vcumhaz = cumsum(df2$nevent/(df2$nrisk*(df2$nrisk-df2$nevent)))
  df2$sesurv = df2$surv*sqrt(vcumhaz)

  # obtain the quantile estimate and confidence interval
  a = lapply(probs, function(p) {
    # Brookmeyer & Crowley confidence interval for quantile
    if (transform == "linear") {
      z = (df2$surv - (1-p))/df2$sesurv
    } else if (transform == "loglog") {
      grad = 1/(df2$surv*log(df2$surv))
      z = (log(-log(df2$surv)) - log(-log(1-p)))/(grad*df2$sesurv)
    } else if (transform == "log") {
      grad = 1/df2$surv
      z = (log(df2$surv) - log(1-p))/(grad*df2$sesurv)
    } else if (transform == "asinsqrt") {
      grad = 1/(2*sqrt(df2$surv*(1-df2$surv)))
      z = (asin(sqrt(df2$surv)) - asin(sqrt(1-p)))/(grad*df2$sesurv)
    } else if (transform == "logit") {
      grad = 1/(df2$surv*(1-df2$surv))
      z = (qlogis(df2$surv) - qlogis(1-p))/(grad*df2$sesurv)
    }

    i = which(abs(z[!is.nan(z)]) <= qnorm((1+cilevel)/2))
    if (length(i) == 0) {
      lower = NA
      upper = NA
    } else {
      lower = df2$time[min(i)]

      if (max(i) <= nrow(df2)) {
        upper = df2$time[max(i)+1]
      } else {
        upper = NA
      }
    }

    if (any(df2$surv < 1 - p)) {
      q = df2$time[min(which(df2$surv < 1 - p))]
      data.frame(prob = p, quantile = q, lower = lower, upper = upper)
    } else {
      data.frame(prob = p, quantile = NA, lower = lower, upper = upper)
    }
  })

  b = do.call("rbind", a)
  b$cilevel = cilevel
  b$transform = transform
  b
}

