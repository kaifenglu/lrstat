#' @title The quantiles of a survival distribution
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


#' @title Profile log-likelihood function for the change points in 
#' piecewise exponential approximation
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


#' @title Piecewise exponential approximation to a survival distribution
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
