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
errorSpent <- function(t, error, sf = "sfOF", sfpar = NA) {
  sapply(t, errorSpentcpp, error = error, sf = sf, sfpar = sfpar)
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

  rtpwexpcpp(n = n, piecewiseSurvivalTime = piecewiseSurvivalTime,
             lambda = lambda, lowerBound = lowerBound)
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


#' @title Adjusted p-Values for Bonferroni-Based Graphical Approaches
#' @description Obtains the adjusted p-values for graphical approaches
#' using weighted Bonferroni tests.
#'
#' @param w The vector of initial weights for elementary hypotheses.
#' @param G The initial transition matrix.
#' @param p The raw p-values for elementary hypotheses.
#'
#' @return A matrix of adjusted p-values.
#'
#' @references
#' Frank Bretz, Willi Maurer, Werner Brannath and Martin Posch. A
#' graphical approach to sequentially rejective multiple test
#' procedures. Statistics in Medicine. 2009; 28:586-604.
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @examples
#'
#' pvalues <- matrix(c(0.01,0.005,0.015,0.022, 0.02,0.015,0.010,0.023),
#'                   nrow=2, ncol=4, byrow=TRUE)
#' w <- c(0.5,0.5,0,0)
#' g <- matrix(c(0,0,1,0,0,0,0,1,0,1,0,0,1,0,0,0),
#'             nrow=4, ncol=4, byrow=TRUE)
#' fadjpbon(w, g, pvalues)
#'
#' @export
fadjpbon <- function(w, G, p) {
  m = length(w)

  if (!is.matrix(p)) {
    p = matrix(p, ncol=m)
  }

  x = fadjpboncpp(w = w, G = G, p = p)
  if (nrow(x) == 1) {
    x = as.vector(x)
  }
  x
}


#' @title Adjusted p-Values for Dunnett-Based Graphical Approaches
#' @description Obtains the adjusted p-values for graphical approaches
#' using weighted Dunnett tests.
#'
#' @param wgtmat The weight matrix for intersection hypotheses.
#' @param p The raw p-values for elementary hypotheses.
#' @param family The matrix of family indicators for elementary hypotheses.
#' @param corr The correlation matrix that should be used for the parametric
#'   test. Can contain NAs for unknown correlations between families.
#'
#' @return A matrix of adjusted p-values.
#'
#' @references
#' Frank Bretz, Martin Posch, Ekkehard Glimm, Florian Klinglmueller,
#' Willi Maurer, and Kornelius Rohmeyer. Graphical approach for multiple
#' comparison procedures using weighted Bonferroni, Simes, or
#' parameter tests. Biometrical Journal. 2011; 53:894-913.
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @examples
#'
#' pvalues <- matrix(c(0.01,0.005,0.015,0.022, 0.02,0.015,0.010,0.023),
#'                   nrow=2, ncol=4, byrow=TRUE)
#' w <- c(0.5,0.5,0,0)
#' g <- matrix(c(0,0,1,0,0,0,0,1,0,1,0,0,1,0,0,0),
#'             nrow=4, ncol=4, byrow=TRUE)
#' wgtmat = fwgtmat(w,g)
#'
#' family = matrix(c(1,1,0,0,0,0,1,1), nrow=2, ncol=4, byrow=TRUE)
#' corr = matrix(c(1,0.5,NA,NA, 0.5,1,NA,NA,
#'                 NA,NA,1,0.5, NA,NA,0.5,1),
#'               nrow = 4, byrow = TRUE)
#' fadjpdun(wgtmat, pvalues, family, corr)
#'
#' @export
fadjpdun <- function(wgtmat, p, family = NULL, corr = NULL) {
  ntests = nrow(wgtmat)
  m = ncol(wgtmat)

  if (!is.matrix(p)) {
    p = matrix(p, ncol=m)
  }

  r = nrow(p)

  if (is.null(family)) {
    family = matrix(1, 1, m)
  } else if (!is.matrix(family)) {
    family = matrix(family, ncol = m)
  }

  if (is.null(corr)) {
    corr = 0.5*diag(m) + 0.5
  }

  pinter = matrix(0, r, ntests)
  incid = matrix(0, ntests, m)
  for (i in 1:ntests) {
    number = ntests - i + 1
    cc = floor(number/2^(m - (1:m))) %% 2
    w = wgtmat[i,]

    J = which(cc == 1)
    J1 = intersect(J, which(w > 0))
    l = nrow(family)

    if (length(J1) > 1) {
      if (r > 1) {
        q = apply(p[,J1]/w[J1], 1, min)
      } else {
        q = min(p[,J1]/w[J1])
      }
    } else {
      q = p[,J1]/w[J1]
    }

    for (k in 1:r) {
      aval = 0
      for (h in 1:l) {
        I_h = which(family[h,] == 1)
        J_h = intersect(J1, I_h)
        if (length(J_h) > 0) {
          sigma = corr[J_h, J_h]
          upper = qnorm(1 - w[J_h]*q[k])
          v = pmvnorm(upper = upper, sigma = sigma, algorithm = "Miwa")
          aval = aval + (1 - v)
        }
      }
      pinter[k,i] = aval
    }

    incid[i,] = cc
  }


  x = matrix(0, r, m)
  for (j in 1:m) {
    ind = matrix(rep(incid[,j], each=r), nrow=r)
    x[,j] = apply(pinter*ind, 1, max)
  }
  x[x > 1] = 1
  x

  if (nrow(x) == 1) {
    x = as.vector(x)
  }
  x
}


#' @title Adjusted p-Values for Simes-Based Graphical Approaches
#' @description Obtains the adjusted p-values for graphical approaches
#' using weighted Simes tests.
#'
#' @param wgtmat The weight matrix for intersection hypotheses.
#' @param p The raw p-values for elementary hypotheses.
#' @param family The matrix of family indicators for elementary hypotheses.
#'
#' @return A matrix of adjusted p-values.
#'
#' @references
#' Frank Bretz, Martin Posch, Ekkehard Glimm, Florian Klinglmueller,
#' Willi Maurer, and Kornelius Rohmeyer. Graphical approach for multiple
#' comparison procedures using weighted Bonferroni, Simes, or
#' parameter tests. Biometrical Journal. 2011; 53:894-913.
#'
#' Kaifeng Lu. Graphical approaches using a Bonferroni mixture of weighted
#' Simes tests. Statistics in Medicine. 2016; 35:4041-4055.
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @examples
#'
#' pvalues <- matrix(c(0.01,0.005,0.015,0.022, 0.02,0.015,0.010,0.023),
#'                   nrow=2, ncol=4, byrow=TRUE)
#' w <- c(0.5,0.5,0,0)
#' g <- matrix(c(0,0,1,0,0,0,0,1,0,1,0,0,1,0,0,0),
#'             nrow=4, ncol=4, byrow=TRUE)
#' wgtmat = fwgtmat(w,g)
#'
#' family = matrix(c(1,1,0,0,0,0,1,1), nrow=2, ncol=4, byrow=TRUE)
#' fadjpsim(wgtmat, pvalues, family)
#'
#' @export
fadjpsim <- function(wgtmat, p, family = NULL) {
  m = ncol(wgtmat)

  if (!is.matrix(p)) {
    p = matrix(p, ncol=m)
  }

  if (is.null(family)) {
    family = matrix(1, 1, m)
  } else if (!is.matrix(family)) {
    family = matrix(family, ncol = m)
  }

  x = fadjpsimcpp(wgtmat = wgtmat, p = p, family = family)
  if (nrow(x) == 1) {
    x = as.vector(x)
  }
  x
}


#' @title Adjusted p-Values for Holm, Hochberg, and Hommel Procedures
#' @description Obtains the adjusted p-values for possibly truncated
#' Holm, Hochberg, and Hommel procedures.
#'
#' @param p The raw p-values for elementary hypotheses.
#' @param test The test to use, e.g., "holm", "hochberg", or
#'   "hommel" (default).
#' @param gamma The value of the truncation parameter. Defaults to 1
#'   for the regular Holm, Hochberg, or Hommel procedure.
#'
#' @return A matrix of adjusted p-values.
#'
#' @references
#' Alex Dmitrienko, Ajit C. Tamhane, and Brian L. Wiens.
#' General multistage gatekeeping procedures.
#' Biometrical Journal. 2008; 5:667-677.
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @examples
#'
#' pvalues <- matrix(c(0.01,0.005,0.015,0.022, 0.02,0.015,0.010,0.023),
#'                   nrow=2, ncol=4, byrow=TRUE)
#' ftrunc(pvalues, "hochberg")
#'
#' @export
ftrunc <- function(p, test = "hommel", gamma = 1) {
  if (!(test == "holm" || test == "hochberg" || test == "hommel")) {
    stop("test must be Holm, Hochberg, or Hommel");
  }

  if (gamma < 0 || gamma > 1) {
    stop("gamma must be within [0, 1]");
  }

  if (!is.matrix(p)) {
    p = matrix(p, nrow=1)
  }

  x = ftrunccpp(p = p, test = test, gamma = gamma)
  if (nrow(x) == 1) {
    x = as.vector(x)
  }
  x
}


#' @title Repeated p-Values for Group Sequential Design
#' @description Obtains the repeated p-values for a group sequential design.
#'
#' @inheritParams param_kMax
#' @inheritParams param_typeAlphaSpending
#' @inheritParams param_parameterAlphaSpending
#' @param maxInformation The target maximum information. Defaults to 1,
#'   in which case, \code{information} represents \code{informationRates}.
#' @param p The raw p-values at look 1 to look \code{k}. It can be a matrix
#'   with \code{k} columns for \code{k <= kMax}.
#' @param information The observed information by look. It can be a matrix
#'   with \code{k} columns.
#' @param spendingTime The error spending time at each analysis, must be
#'   increasing and less than or equal to 1. Defaults to \code{NULL},
#'   in which case, it is the same as \code{informationRates} derived from
#'   \code{information} and \code{maxInformation}. It can be a matrix with
#'   \code{k} columns.
#'
#' @return The repeated p-values at look 1 to look \code{k}.
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @examples
#'
#' # Example 1: informationRates different from spendingTime
#' repeatedPValue(kMax = 3, typeAlphaSpending = "sfOF",
#'                maxInformation = 800,
#'                p = c(0.2, 0.15, 0.1),
#'                information = c(529, 700, 800),
#'                spendingTime = c(0.6271186, 0.8305085, 1))
#'
#' # Example 2: Maurer & Bretz (2013), current look is not the last look
#' repeatedPValue(kMax = 3, typeAlphaSpending = "sfOF",
#'                p = matrix(c(0.0062, 0.017,
#'                             0.009, 0.13,
#'                             0.0002, 0.0035,
#'                             0.002, 0.06),
#'                           nrow=4, ncol=2),
#'                information = c(1/3, 2/3))
#'
#' @export
repeatedPValue <- function(kMax,
                           typeAlphaSpending = "sfOF",
                           parameterAlphaSpending = NA,
                           maxInformation = 1,
                           p,
                           information,
                           spendingTime = NULL) {

  if (is.matrix(p)) {
    p1 = p
  } else {
    p1 = matrix(p, nrow=1)
  }

  if (is.matrix(information)) {
    information1 = information
  } else {
    information1 = matrix(information, nrow=1)
  }

  if (is.null(spendingTime)) {
    spendingTime1 = matrix(0, 1, 1)
  } else if (is.matrix(spendingTime)) {
    spendingTime1 = spendingTime
  } else {
    spendingTime1 = matrix(spendingTime, nrow=1)
  }

  repp1 = repeatedPValuecpp(
    kMax = kMax, typeAlphaSpending = typeAlphaSpending,
    parameterAlphaSpending = parameterAlphaSpending,
    maxInformation = maxInformation, p = p1,
    information = information1, spendingTime = spendingTime1)

  if (is.vector(p)) { # convert the result to a vector
    repp = c(repp1)
  } else {
    repp = repp1
  }

  repp
}


#' @title Group Sequential Trials Using Bonferroni-Based Graphical
#' Approaches
#'
#' @description Obtains the test results for group sequential trials using
#' graphical approaches based on weighted Bonferroni tests.
#'
#' @param w The vector of initial weights for elementary hypotheses.
#' @param G The initial transition matrix.
#' @param alpha The significance level. Defaults to 0.025.
#' @inheritParams param_kMax
#' @param typeAlphaSpending The vector of alpha spending functions.
#'   Each element is one of the following:
#'   "OF" for O'Brien-Fleming boundaries,
#'   "P" for Pocock boundaries, "WT" for Wang & Tsiatis boundaries,
#'   "sfOF" for O'Brien-Fleming type spending function,
#'   "sfP" for Pocock type spending function,
#'   "sfKD" for Kim & DeMets spending function,
#'   "sfHSD" for Hwang, Shi & DeCani spending function,
#'   and "none" for no early efficacy stopping.
#'   Defaults to "sfOF" if not provided.
#' @param parameterAlphaSpending The vector of parameter values for the
#'   alpha spending functions. Each element corresponds to the value of
#'   Delta for "WT", rho for "sfKD", or gamma for "sfHSD".
#'   Defaults to missing if not provided.
#' @param incidenceMatrix The incidence matrix indicating whether the
#'   specific hypothesis will be tested at the given look. The number of
#'   columns of incidenceMatrix must be equal to the maximum number of
#'   study looks (\code{kMax}). If not provided, defaults to testing each
#'   hypothesis at all study looks.
#' @param maxInformation The vector of target maximum information for each
#'   hypothesis. Defaults to a vector of 1s if not provided.
#' @param information The matrix of observed information for each hypothesis
#'   by study look.
#' @param p The matrix of raw p-values for each hypothesis by study look.
#' @param spendingTime The spending time for alpha spending by study look.
#'   If not provided, it is the same as \code{informationRates} calculated
#'   from \code{information} and \code{maxInformation}.
#'
#' @return A vector to indicate the first look the specific hypothesis is
#'   rejected (0 if the hypothesis is not rejected).
#'
#' @references
#' Willi Maurer and Frank Bretz. Multiple testing in group sequential
#' trials using graphical approaches. Statistics in Biopharmaceutical
#' Research. 2013; 5:311-320.
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @examples
#'
#' # Case study from Maurer & Bretz (2013)
#'
#' fseqbon(
#'   w = c(0.5, 0.5, 0, 0),
#'   G = matrix(c(0, 0.5, 0.5, 0,  0.5, 0, 0, 0.5,
#'                0, 1, 0, 0,  1, 0, 0, 0),
#'              nrow=4, ncol=4, byrow=TRUE),
#'   alpha = 0.025,
#'   kMax = 3,
#'   typeAlphaSpending = rep("sfOF", 4),
#'   maxInformation = rep(1, 4),
#'   p = matrix(c(0.0062, 0.017, 0.009, 0.13,
#'                0.0002, 0.0035, 0.002, 0.06),
#'              nrow=4, ncol=2),
#'   information = matrix(c(rep(1/3, 4), rep(2/3, 4)),
#'                        nrow=4, ncol=2))
#'
#'
#' @export
fseqbon <- function(w, G, alpha = 0.025, kMax,
                    typeAlphaSpending = NULL,
                    parameterAlphaSpending = NULL,
                    incidenceMatrix = NULL,
                    maxInformation = NULL,
                    p, information, spendingTime = NULL) {
  m = length(w)

  if (is.null(typeAlphaSpending)) {
    typeAlphaSpending = rep("sfOF", m)
  }

  if (is.null(parameterAlphaSpending)) {
    parameterAlphaSpending = rep(NA, m)
  }

  if (is.null(incidenceMatrix)) {
    incidenceMatrix = matrix(1, nrow=m, ncol=kMax)
  }

  if (is.null(maxInformation)) {
    maxInformation = rep(1, m)
  }

  if (is.null(spendingTime)) {
    spendingTime = matrix(0, 1, 1)
  }

  fseqboncpp(w = w, G = G, alpha = alpha, kMax = kMax,
             typeAlphaSpending = typeAlphaSpending,
             parameterAlphaSpending = parameterAlphaSpending,
             incidenceMatrix = incidenceMatrix,
             maxInformation = maxInformation, p = p,
             information = information, spendingTime = spendingTime)
}


#' @title Adjusted p-Values for Stepwise Testing Procedures for Two Sequences
#' @description Obtains the adjusted p-values for the stepwise gatekeeping
#' procedures for multiplicity problems involving two sequences of
#' hypotheses.
#'
#' @param p The raw p-values for elementary hypotheses.
#' @param gamma The truncation parameters for each family. The truncation
#'   parameter for the last family is automatically set to 1.
#' @param test The component multiple testing procedure. It is either "Holm"
#'   or "Hochberg", and it defaults to "Hochberg".
#' @param retest Whether to allow retesting. It defaults to \code{TRUE}.
#'
#' @return A matrix of adjusted p-values.
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @examples
#'
#' p = c(0.0194, 0.0068, 0.0271, 0.0088, 0.0370, 0.0018, 0.0814, 0.0066)
#' gamma = c(0.6, 0.6, 0.6, 1)
#' fstp2seq(p, gamma, test="hochberg", retest=1)
#'
#' @export
fstp2seq <- function(p, gamma, test="hochberg", retest=TRUE) {
  if (!is.matrix(p)) {
    p = matrix(p, nrow=1)
  }
  m = ncol(p)

  if (m %% 2 != 0) {
    stop("p must have an even number of columns")
  }

  if (!all(gamma >= 0 & gamma <= 1)) {
    stop("gamma must lie between 0 and 1");
  }

  if (m == 2) {
    gamma = 1
  } else if (length(gamma) == m/2 - 1) {
    gamma = c(gamma, 1)
  } else if (length(gamma) != m/2) {
    stop("The number of families must be half of the number of hypotheses")
  } else {
    gamma[m/2] = 1
  }

  if (!(tolower(test) %in% c("holm", "hochberg"))) {
    stop("test must be either Holm or Hochberg")
  }

  x = fstp2seqcpp(p = p, gamma = gamma, test = test, retest = retest)
  if (nrow(x) == 1) {
    x = as.vector(x)
  }
  x
}

#' @title Adjusted p-Values for Standard Mixture Gatekeeping Procedures
#' @description Obtains the adjusted p-values for the standard gatekeeping
#' procedures for multiplicity problems involving serial and parallel
#' logical restrictions.
#'
#' @param p The raw p-values for elementary hypotheses.
#' @param family The matrix of family indicators for the hypotheses.
#' @param serial The matrix of serial rejection set for the hypotheses.
#' @param parallel The matrix of parallel rejection set for the hypotheses.
#' @param gamma The truncation parameters for each family. The truncation
#'   parameter for the last family is automatically set to 1.
#' @param test The component multiple testing procedure. Options include
#'   "holm", "hochberg", or "hommel". Defaults to "hommel".
#' @param exhaust Whether to use alpha-exhausting component testing procedure
#'   for the last family with active hypotheses. It defaults to \code{TRUE}.
#'
#' @return A matrix of adjusted p-values.
#'
#' @references
#' Alex Dmitrienko and Ajit C Tamhane. Mixtures of multiple testing
#' procedures for gatekeeping applications in clinical trials.
#' Statistics in Medicine. 2011; 30(13):1473–1488.
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @examples
#'
#' p = c(0.0194, 0.0068, 0.0271, 0.0088, 0.0370, 0.0018, 0.0814, 0.0066)
#' family = matrix(c(1, 1, 0, 0, 0, 0, 0, 0,
#'                   0, 0, 1, 1, 0, 0, 0, 0,
#'                   0, 0, 0, 0, 1, 1, 0, 0,
#'                   0, 0, 0, 0, 0, 0, 1, 1),
#'                 nrow=4, byrow=TRUE)
#'
#' serial = matrix(c(0, 0, 0, 0, 0, 0, 0, 0,
#'                   0, 0, 0, 0, 0, 0, 0, 0,
#'                   1, 0, 0, 0, 0, 0, 0, 0,
#'                   0, 1, 0, 0, 0, 0, 0, 0,
#'                   0, 0, 1, 0, 0, 0, 0, 0,
#'                   0, 0, 0, 1, 0, 0, 0, 0,
#'                   0, 0, 0, 0, 1, 0, 0, 0,
#'                   0, 0, 0, 0, 0, 1, 0, 0),
#'                 nrow=8, byrow=TRUE)
#'
#' parallel = matrix(0, 8, 8)
#' gamma = c(0.6, 0.6, 0.6, 1)
#' fstdmix(p, family, serial, parallel, gamma, test = "hommel",
#'         exhaust = FALSE)
#'
#' @export
fstdmix <- function(p, family = NULL, serial, parallel,
                    gamma, test = "hommel", exhaust = TRUE) {
  if (!is.matrix(p)) {
    p = matrix(p, nrow=1)
  }
  m = ncol(p)

  if (is.null(family)) {
    family = matrix(1, 1, m)
  } else if (!is.matrix(family)) {
    family = matrix(family, ncol = m)
  }
  k = nrow(family)

  if (ncol(family) != m) {
    stop("number of columns of family must be the number of hypotheses")
  }
  if (any(family != 0 & family != 1)) {
    stop("elements of family must be 0 or 1")
  }
  if (any(colSums(family) != 1)) {
    stop("elements of family must sum to 1 for each column")
  }

  if (nrow(serial) != m || ncol(serial) != m) {
    stop(paste("serial must be a square matrix with the number of",
               "rows/columns equal to the number of hypotheses"))
  }
  if (any(serial != 0 & serial != 1)) {
    stop("elements of serial must be 0 or 1")
  }

  if (nrow(parallel) != m || ncol(parallel) != m) {
    stop(paste("parallel must be a square matrix with the number of",
               "rows/columns equal to the number of hypotheses"))
  }
  if (any(parallel != 0 & parallel != 1)) {
    stop("elements of parallel must be 0 or 1")
  }

  if (!all(gamma >= 0 & gamma <= 1)) {
    stop("gamma must lie between 0 and 1")
  }

  if (k == 1) {
    gamma = 1
  } else if (length(gamma) == k - 1) {
    gamma = c(gamma, 1)
  } else if (length(gamma) != k) {
    stop("The length of gamma must match the number of families")
  } else {
    gamma[k] = 1
  }

  if (!(tolower(test) %in% c("holm", "hochberg", "hommel"))) {
    stop("test must be either Holm, Hochberg, or Hommel")
  }

  x = fstdmixcpp(p = p, family = family, serial = serial,
                 parallel = parallel, gamma = gamma, test = test,
                 exhaust = exhaust)

  if (nrow(x) == 1) {
    x = as.vector(x)
  }
  x
}


#' @title Adjusted p-Values for Modified Mixture Gatekeeping Procedures
#' @description Obtains the adjusted p-values for the modified gatekeeping
#' procedures for multiplicity problems involving serial and parallel
#' logical restrictions.
#'
#' @param p The raw p-values for elementary hypotheses.
#' @param family The matrix of family indicators for the hypotheses.
#' @param serial The matrix of serial rejection set for the hypotheses.
#' @param parallel The matrix of parallel rejection set for the hypotheses.
#' @param gamma The truncation parameters for each family. The truncation
#'   parameter for the last family is automatically set to 1.
#' @param test The component multiple testing procedure. Options include
#'   "holm", "hochberg", or "hommel". Defaults to "hommel".
#' @param exhaust Whether to use alpha-exhausting component testing procedure
#'   for the last family with active hypotheses. It defaults to \code{TRUE}.
#'
#' @return A matrix of adjusted p-values.
#'
#' @references
#' Alex Dmitrienko, George Kordzakhia, and Thomas Brechenmacher.
#' Mixture-based gatekeeping procedures for multiplicity problems with
#' multiple sequences of hypotheses. Journal of Biopharmaceutical
#' Statistics. 2016; 26(4):758–780.
#'
#' George Kordzakhia, Thomas Brechenmacher, Eiji Ishida, Alex Dmitrienko,
#' Winston Wenxiang Zheng, and David Fuyuan Li. An enhanced mixture method
#' for constructing gatekeeping procedures in clinical trials.
#' Journal of Biopharmaceutical Statistics. 2018; 28(1):113–128.
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @examples
#'
#' p = c(0.0194, 0.0068, 0.0271, 0.0088, 0.0370, 0.0018, 0.0814, 0.0066)
#' family = matrix(c(1, 1, 0, 0, 0, 0, 0, 0,
#'                   0, 0, 1, 1, 0, 0, 0, 0,
#'                   0, 0, 0, 0, 1, 1, 0, 0,
#'                   0, 0, 0, 0, 0, 0, 1, 1),
#'                 nrow=4, byrow=TRUE)
#'
#' serial = matrix(c(0, 0, 0, 0, 0, 0, 0, 0,
#'                   0, 0, 0, 0, 0, 0, 0, 0,
#'                   1, 0, 0, 0, 0, 0, 0, 0,
#'                   0, 1, 0, 0, 0, 0, 0, 0,
#'                   0, 0, 1, 0, 0, 0, 0, 0,
#'                   0, 0, 0, 1, 0, 0, 0, 0,
#'                   0, 0, 0, 0, 1, 0, 0, 0,
#'                   0, 0, 0, 0, 0, 1, 0, 0),
#'                 nrow=8, byrow=TRUE)
#'
#' parallel = matrix(0, 8, 8)
#' gamma = c(0.6, 0.6, 0.6, 1)
#' fmodmix(p, family, serial, parallel, gamma, test = "hommel", exhaust = TRUE)
#'
#' @export
fmodmix <- function(p, family = NULL, serial, parallel,
                    gamma, test = "hommel", exhaust = TRUE) {
  if (!is.matrix(p)) {
    p = matrix(p, nrow=1)
  }
  m = ncol(p)

  if (is.null(family)) {
    family = matrix(1, 1, m)
  } else if (!is.matrix(family)) {
    family = matrix(family, ncol = m)
  }
  k = nrow(family)

  if (ncol(family) != m) {
    stop("number of columns of family must be the number of hypotheses")
  }
  if (any(family != 0 & family != 1)) {
    stop("elements of family must be 0 or 1")
  }
  if (any(colSums(family) != 1)) {
    stop("elements of family must sum to 1 for each column")
  }

  if (nrow(serial) != m || ncol(serial) != m) {
    stop(paste("serial must be a square matrix with the number of",
               "rows/columns equal to the number of hypotheses"))
  }
  if (any(serial != 0 & serial != 1)) {
    stop("elements of serial must be 0 or 1")
  }

  if (nrow(parallel) != m || ncol(parallel) != m) {
    stop(paste("parallel must be a square matrix with the number of",
               "rows/columns equal to the number of hypotheses"))
  }
  if (any(parallel != 0 & parallel != 1)) {
    stop("elements of parallel must be 0 or 1")
  }

  if (!all(gamma >=0 & gamma <= 1)) {
    stop("gamma must lie between 0 and 1");
  }

  if (k == 1) {
    gamma = 1
  } else if (length(gamma) == k - 1) {
    gamma = c(gamma, 1)
  } else if (length(gamma) != k) {
    stop("The length of gamma must match the number of families")
  } else {
    gamma[k] = 1
  }

  if (!(tolower(test) %in% c("holm", "hochberg", "hommel"))) {
    stop("test must be either Holm, Hochberg, or Hommel")
  }

  x = fmodmixcpp(p = p, family = family, serial = serial,
                 parallel = parallel, gamma = gamma, test = test,
                 exhaust = exhaust)

  if (nrow(x) == 1) {
    x = as.vector(x)
  }
  x
}



#' @title Efficacy Boundaries for Group Sequential Design
#' @description Obtains the efficacy stopping boundaries for a group
#' sequential design.
#'
#' @param k Look number for the current analysis.
#' @param informationRates Information rates up to the current look. Must be
#'   increasing and less than or equal to 1.
#' @inheritParams param_alpha
#' @inheritParams param_typeAlphaSpending
#' @inheritParams param_parameterAlphaSpending
#' @inheritParams param_userAlphaSpending
#' @param spendingTime A vector of length \code{k} for the error spending
#'   time at each analysis. Must be increasing and less than or equal to 1.
#'   Defaults to missing, in which case, it is the same as
#'   \code{informationRates}.
#' @inheritParams param_efficacyStopping
#'
#' @details
#' If \code{typeAlphaSpending} is "OF", "P", or "WT", then the boundaries
#' will be based on equally spaced looks.
#'
#' @return A numeric vector of critical values up to the current look.
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @examples
#'
#' getBound(k = 2, informationRates = c(0.5,1),
#'          alpha = 0.025, typeAlphaSpending = "sfOF")
#'
#' @export
getBound <- function(k = NA, informationRates = NA, alpha = 0.025,
                     typeAlphaSpending = "sfOF", parameterAlphaSpending = NA,
                     userAlphaSpending = NA, spendingTime = NA,
                     efficacyStopping = NA) {
  getBoundcpp(k = k, informationRates = informationRates, alpha = alpha,
              typeAlphaSpending = typeAlphaSpending,
              parameterAlphaSpending = parameterAlphaSpending,
              userAlphaSpending = userAlphaSpending,
              spendingTime = spendingTime,
              efficacyStopping = efficacyStopping)
}


#' @title Parametric Regression Models for Failure Time Data
#' @description Obtains the parameter estimates from parametric
#' regression models with uncensored, right censored, left censored, or
#' interval censored data.
#'
#' @param data The input data frame that contains the following variables:
#'
#'   * \code{rep}: The replication for by-group processing.
#'
#'   * \code{stratum}: The stratum.
#'
#'   * \code{time}: The follow-up time for right censored data, or
#'     the left end of each interval for interval censored data.
#'
#'   * \code{time2}: The right end of each interval for interval
#'     censored data.
#'
#'   * \code{event}: The event indicator, 1=event, 0=no event.
#'
#'   * \code{covariates}: The values of baseline covariates.
#'
#'   * \code{weight}: The weight for each observation.
#'
#'   * \code{offset}: The offset for each observation.
#'
#'   * \code{id}: The optional subject ID to group the score residuals
#'     in computing the robust sandwich variance.
#'
#' @param rep The name(s) of the replication variable(s) in the input data.
#' @param stratum The name(s) of the stratum variable(s) in the input data.
#' @param time The name of the time variable or the left end of each
#'   interval for interval censored data in the input data.
#' @param time2 The name of the right end of each interval for
#'   interval censored data in the input data.
#' @param event The name of the event variable in the input data
#'   for right censored data.
#' @param covariates The vector of names of baseline covariates
#'   in the input data.
#' @param weight The name of the weight variable in the input data.
#' @param offset The name of the offset variable in the input data.
#' @param id The name of the id variable in the input data.
#' @param dist The assumed distribution for time to event. Options include
#'   "exponential", "weibull", "lognormal", and "loglogistic" to be
#'   modeled on the log-scale, and "normal" and "logistic" to be modeled
#'   on the original scale.
#' @param init A vector of initial values for the model parameters,
#'   including regression coefficients and the log scale parameter.
#'   By default, initial values are derived from an intercept-only model.
#'   If this approach fails, ordinary least squares (OLS) estimates,
#'   ignoring censoring, are used instead.
#' @param robust Whether a robust sandwich variance estimate should be
#'   computed. In the presence of the id variable, the score residuals
#'   will be aggregated for each id when computing the robust sandwich
#'   variance estimate.
#' @param plci Whether to obtain profile likelihood confidence interval.
#' @param alpha The two-sided significance level.
#' @param maxiter The maximum number of iterations.
#' @param eps The tolerance to declare convergence.
#'
#' @details There are two ways to specify the model, one for right censored
#' data through the time and event variables, and the other for interval
#' censored data through the time (lower) and time2 (upper) variables.
#' For the second form, we follow the convention used in SAS PROC LIFEREG:
#'
#' * If lower is not missing, upper is not missing, and lower is equal
#'   to upper, then there is no censoring and the event occurred at
#'   time lower.
#'
#' * If lower is not missing, upper is not missing, and lower < upper,
#'   then the event time is censored within the interval (lower, upper).
#'
#' * If lower is missing, but upper is not missing, then upper will be
#'   used as the left censoring value.
#'
#' * If lower is not missing, but upper is missing, then lower will be
#'   used as the right censoring value.
#'
#' * If lower is not missing, upper is not missing, but lower > upper,
#'   or if both lower and upper are missing, then the observation will
#'   not be used.
#'
#' @return A list with the following components:
#'
#' * \code{sumstat}: The data frame of summary statistics of model fit
#'   with the following variables:
#'
#'     - \code{n}: The number of observations.
#'
#'     - \code{nevents}: The number of events.
#'
#'     - \code{loglik0}: The log-likelihood under null.
#'
#'     - \code{loglik1}: The maximum log-likelihood.
#'
#'     - \code{niter}: The number of Newton-Raphson iterations.
#'
#'     - \code{dist}: The assumed distribution.
#'
#'     - \code{p}: The number of parameters, including the intercept,
#'       regression coefficients associated with the covariates, and
#'       the log scale parameters for the strata.
#'
#'     - \code{nvar}: The number of regression coefficients associated
#'       with the covariates (excluding the intercept).
#'
#'     - \code{robust}: Whether the robust sandwich variance estimate
#'       is requested.
#'
#'     - \code{fail}: Whether the model fails to converge.
#'
#'     - \code{rep}: The replication.
#'
#' * \code{parest}: The data frame of parameter estimates with the
#'   following variables:
#'
#'     - \code{param}: The name of the covariate for the parameter estimate.
#'
#'     - \code{beta}: The parameter estimate.
#'
#'     - \code{sebeta}: The standard error of parameter estimate.
#'
#'     - \code{z}: The Wald test statistic for the parameter.
#'
#'     - \code{expbeta}: The exponentiated parameter estimate.
#'
#'     - \code{vbeta}: The covariance matrix for parameter estimates.
#'
#'     - \code{lower}: The lower limit of confidence interval.
#'
#'     - \code{upper}: The upper limit of confidence interval.
#'
#'     - \code{p}: The p-value from the chi-square test.
#'
#'     - \code{method}: The method to compute the confidence interval and
#'       p-value.
#'
#'     - \code{sebeta_naive}: The naive standard error of parameter estimate
#'       if robust variance is requested.
#'
#'     - \code{vbeta_naive}: The naive covariance matrix for parameter
#'       estimates if robust variance is requested.
#'
#'     - \code{rep}: The replication.
#'
#' * \code{linear_predictors}: The vector of linear predictors.
#'
#' * \code{p}: The number of parameters.
#'
#' * \code{nvar}: The number of columns of the design matrix excluding
#'   the intercept.
#'
#' * \code{param}: The parameter names.
#'
#' * \code{beta}: The parameter estimate.
#'
#' * \code{vbeta}: The covariance matrix for parameter estimates.
#'
#' * \code{vbeta_naive}: The naive covariance matrix for parameter estimates.
#'
#' * \code{terms}: The terms object.
#'
#' * \code{xlevels}: A record of the levels of the factors used in fitting.
#'
#' * \code{settings}: A list containing the input parameter values.
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @references
#' John D. Kalbfleisch and Ross L. Prentice.
#' The Statistical Analysis of Failure Time Data.
#' Wiley: New York, 1980.
#'
#' @examples
#'
#' library(dplyr)
#'
#' # right censored data
#' (fit1 <- liferegr(
#'   data = rawdata %>% mutate(treat = 1*(treatmentGroup == 1)),
#'   rep = "iterationNumber", stratum = "stratum",
#'   time = "timeUnderObservation", event = "event",
#'   covariates = "treat", dist = "weibull"))
#'
#' # tobit regression for left censored data
#' (fit2 <- liferegr(
#'   data = tobin %>% mutate(time = ifelse(durable>0, durable, NA)),
#'   time = "time", time2 = "durable",
#'   covariates = c("age", "quant"), dist = "normal"))
#'
#' @export
liferegr <- function(data, rep = "", stratum = "",
                     time = "time", time2 = "", event = "event",
                     covariates = "", weight = "", offset = "",
                     id = "", dist = "weibull", init = NA_real_,
                     robust = FALSE, plci = FALSE, alpha = 0.05,
                     maxiter = 50, eps = 1.0e-9) {
  rownames(data) = NULL

  elements = c(rep, stratum, covariates, weight, offset, id)
  elements = unique(elements[elements != "" & elements != "none"])
  if (!(length(elements) == 0)) {
    fml = formula(paste("~", paste(elements, collapse = "+")))
    mf = model.frame(fml, data = data, na.action = na.omit)
  } else {
    mf = model.frame(formula("~1"), data = data)
  }

  rownum = as.integer(rownames(mf))
  df = data[rownum,]

  nvar = length(covariates)
  if (missing(covariates) || is.null(covariates) || (nvar == 1 && (
    covariates[1] == "" || tolower(covariates[1]) == "none"))) {
    p = 0
    t1 = terms(formula("~1"))
  } else {
    fml1 = formula(paste("~", paste(covariates, collapse = "+")))
    p = length(rownames(attr(terms(fml1), "factors")))
    t1 = terms(fml1)
  }

  if (p >= 1) {
    mf1 <- model.frame(fml1, data = df, na.action = na.pass)
    mm <- model.matrix(fml1, mf1)
    xlevels = mf1$xlev
    param = colnames(mm)
    colnames(mm) = make.names(colnames(mm))
    varnames = colnames(mm)[-1]
    for (i in 1:length(varnames)) {
      if (!(varnames[i] %in% names(df))) {
        df[,varnames[i]] = mm[,varnames[i]]
      }
    }
  } else {
    xlevels = NULL
    param = "(Intercept)"
    varnames = ""
  }

  fit <- liferegcpp(data = df, rep = rep, stratum = stratum, time = time,
                    time2 = time2, event = event, covariates = varnames,
                    weight = weight, offset = offset, id = id, dist = dist,
                    init = init, robust = robust, plci = plci,
                    alpha = alpha, maxiter = maxiter, eps = eps)

  fit$p <- fit$sumstat$p[1]
  fit$nvar <- fit$sumstat$nvar[1]

  if (fit$p > 0) {
    par = fit$parest$param[1:fit$p]
    if (length(par) > length(param)) {
      fit$param = c(param, par[(1+length(param)):length(par)])
    } else {
      fit$param = param
    }

    fit$beta = fit$parest$beta
    names(fit$beta) = rep(fit$param, length(fit$beta)/fit$p)

    if (fit$p > 1) {
      fit$vbeta = as.matrix(fit$parest[, paste0("vbeta.", 1:fit$p)])
      if (robust) {
        fit$vbeta_naive = as.matrix(fit$parest[, paste0("vbeta_naive.",
                                                        1:fit$p)])
      }
    } else {
      fit$vbeta = as.matrix(fit$parest[, "vbeta"])
      if (robust) {
        fit$vbeta_naive = as.matrix(fit$parest[, "vbeta_naive"])
      }
    }

    dimnames(fit$vbeta) = list(names(fit$beta), fit$param)
    if (robust) {
      dimnames(fit$vbeta_naive) = list(names(fit$beta), fit$param)
    }
  }

  fit$terms = t1
  if (fit$p > 0) fit$xlevels = xlevels

  fit$settings <- list(
    data = data, rep = rep, stratum = stratum, time = time, time2 = time2,
    event = event, covariates = covariates, weight = weight, offset = offset,
    id = id, dist = dist, init = init, robust = robust, plci = plci,
    alpha = alpha, maxiter = maxiter, eps = eps
  )

  class(fit) = "liferegr"
  fit
}



#' @title Proportional Hazards Regression Models
#' @description Obtains the hazard ratio estimates from the proportional
#' hazards regression model with right censored or counting process data.
#'
#' @param data The input data frame that contains the following variables:
#'
#'   * \code{rep}: The replication for by-group processing.
#'
#'   * \code{stratum}: The stratum.
#'
#'   * \code{time}: The follow-up time for right censored data, or
#'     the left end of each interval for counting process data.
#'
#'   * \code{time2}: The right end of each interval for counting process
#'     data. Intervals are assumed to be open on the left
#'     and closed on the right, and event indicates whether an event
#'     occurred at the right end of each interval.
#'
#'   * \code{event}: The event indicator, 1=event, 0=no event.
#'
#'   * \code{covariates}: The values of baseline covariates (and
#'     time-dependent covariates in each interval for counting
#'     process data).
#'
#'   * \code{weight}: The weight for each observation.
#'
#'   * \code{offset}: The offset for each observation.
#'
#'   * \code{id}: The optional subject ID for counting process data
#'     with time-dependent covariates.
#'
#' @param rep The name(s) of the replication variable(s) in the input data.
#' @param stratum The name(s) of the stratum variable(s) in the input data.
#' @param time The name of the time variable or the left end of each
#'   interval for counting process data in the input data.
#' @param time2 The name of the right end of each interval for counting
#'   process data in the input data.
#' @param event The name of the event variable in the input data.
#' @param covariates The vector of names of baseline and time-dependent
#'   covariates in the input data.
#' @param weight The name of the weight variable in the input data.
#' @param offset The name of the offset variable in the input data.
#' @param id The name of the id variable in the input data.
#' @param ties The method for handling ties, either "breslow" or
#'   "efron" (default).
#' @param init The vector of initial values. Defaults to zero for all
#'   variables.
#' @param robust Whether a robust sandwich variance estimate should be
#'   computed. In the presence of the id variable, the score residuals
#'   will be aggregated for each id when computing the robust sandwich
#'   variance estimate.
#' @param est_basehaz Whether to estimate the baseline hazards.
#'   Defaults to \code{TRUE}.
#' @param est_resid Whether to estimate the martingale residuals.
#'   Defaults to \code{TRUE}.
#' @param firth Whether to use Firth’s penalized likelihood method.
#'   Defaults to \code{FALSE}.
#' @param plci Whether to obtain profile likelihood confidence interval.
#' @param alpha The two-sided significance level.
#' @param maxiter The maximum number of iterations.
#' @param eps The tolerance to declare convergence.
#'
#' @return A list with the following components:
#'
#' * \code{sumstat}: The data frame of summary statistics of model fit
#'   with the following variables:
#'
#'     - \code{n}: The number of observations.
#'
#'     - \code{nevents}: The number of events.
#'
#'     - \code{loglik0}: The (penalized) log-likelihood under null.
#'
#'     - \code{loglik1}: The maximum (penalized) log-likelihood.
#'
#'     - \code{scoretest}: The score test statistic.
#'
#'     - \code{niter}: The number of Newton-Raphson iterations.
#'
#'     - \code{ties}: The method for handling ties, either "breslow" or
#'       "efron".
#'
#'     - \code{p}: The number of columns of the Cox model design matrix.
#'
#'     - \code{robust}: Whether to use the robust variance estimate.
#'
#'     - \code{firth}: Whether to use Firth's penalized likelihood method.
#'
#'     - \code{fail}: Whether the model fails to converge.
#'
#'     - \code{loglik0_unpenalized}: The unpenalized log-likelihood under null.
#'
#'     - \code{loglik1_unpenalized}: The maximum unpenalized log-likelihood.
#'
#'     - \code{rep}: The replication.
#'
#' * \code{parest}: The data frame of parameter estimates with the
#'   following variables:
#'
#'     - \code{param}: The name of the covariate for the parameter estimate.
#'
#'     - \code{beta}: The log hazard ratio estimate.
#'
#'     - \code{sebeta}: The standard error of log hazard ratio estimate.
#'
#'     - \code{z}: The Wald test statistic for log hazard ratio.
#'
#'     - \code{expbeta}: The hazard ratio estimate.
#'
#'     - \code{vbeta}: The covariance matrix for parameter estimates.
#'
#'     - \code{lower}: The lower limit of confidence interval.
#'
#'     - \code{upper}: The upper limit of confidence interval.
#'
#'     - \code{p}: The p-value from the chi-square test.
#'
#'     - \code{method}: The method to compute the confidence interval and
#'       p-value.
#'
#'     - \code{sebeta_naive}: The naive standard error of log hazard ratio
#'       estimate if robust variance is requested.
#'
#'     - \code{vbeta_naive}: The naive covariance matrix for parameter
#'       estimates if robust variance is requested.
#'
#'     - \code{rep}: The replication.
#'
#' * \code{basehaz}: The data frame of baseline hazards with the following
#'   variables (if est_basehaz is TRUE):
#'
#'     - \code{time}: The observed event time.
#'
#'     - \code{nrisk}: The number of patients at risk at the time point.
#'
#'     - \code{nevent}: The number of events at the time point.
#'
#'     - \code{haz}: The baseline hazard at the time point.
#'
#'     - \code{varhaz}: The variance of the baseline hazard at the time point
#'       assuming the parameter beta is known.
#'
#'     - \code{gradhaz}: The gradient of the baseline hazard with respect to
#'       beta at the time point (in the presence of covariates).
#'
#'     - \code{stratum}: The stratum.
#'
#'     - \code{rep}: The replication.
#'
#' * \code{residuals}: The martingale residuals.
#'
#' * \code{linear_predictors}: The vector of linear predictors.
#'
#' * \code{p}: The number of parameters.
#'
#' * \code{param}: The parameter names.
#'
#' * \code{beta}: The parameter estimate.
#'
#' * \code{vbeta}: The covariance matrix for parameter estimates.
#'
#' * \code{vbeta_naive}: The naive covariance matrix for parameter estimates.
#'
#' * \code{terms}: The terms object.
#'
#' * \code{xlevels}: A record of the levels of the factors used in fitting.
#'
#' * \code{settings}: A list containing the input parameter values.
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @references
#' Per K. Anderson and Richard D. Gill.
#' Cox's regression model for counting processes, a large sample study.
#' Annals of Statistics 1982; 10:1100-1120.
#'
#' Terry M. Therneau and Patricia M. Grambsch.
#' Modeling Survival Data: Extending the Cox Model.
#' Springer-Verlag, 2000.
#'
#' @examples
#'
#' library(dplyr)
#'
#' # Example 1 with right-censored data
#' (fit1 <- phregr(
#'   data = rawdata %>% mutate(treat = 1*(treatmentGroup == 1)),
#'   rep = "iterationNumber", stratum = "stratum",
#'   time = "timeUnderObservation", event = "event",
#'   covariates = "treat", est_basehaz = FALSE, est_resid = FALSE))
#'
#' # Example 2 with counting process data and robust variance estimate
#' (fit2 <- phregr(
#'   data = heart %>% mutate(rx = as.numeric(transplant) - 1),
#'   time = "start", time2 = "stop", event = "event",
#'   covariates = c("rx", "age"), id = "id",
#'   robust = TRUE, est_basehaz = TRUE, est_resid = TRUE))
#'
#' @export
phregr <- function(data, rep = "", stratum = "",
                   time = "time", time2 = "", event = "event",
                   covariates = "", weight = "", offset = "",
                   id = "", ties = "efron",
                   init = NA_real_,  robust = FALSE,
                   est_basehaz = TRUE, est_resid = TRUE, firth = FALSE,
                   plci = FALSE, alpha = 0.05,
                   maxiter = 50, eps = 1.0e-9) {

  rownames(data) = NULL

  elements = c(rep, stratum, time, event, covariates, weight, offset, id)
  elements = unique(elements[elements != "" & elements != "none"])
  fml = formula(paste("~", paste(elements, collapse = "+")))
  mf = model.frame(fml, data = data, na.action = na.omit)

  rownum = as.integer(rownames(mf))
  df = data[rownum,]

  nvar = length(covariates)
  if (missing(covariates) || is.null(covariates) || (nvar == 1 && (
    covariates[1] == "" || tolower(covariates[1]) == "none"))) {
    p = 0
    t1 = terms(formula("~1"))
  } else {
    fml1 = formula(paste("~", paste(covariates, collapse = "+")))
    p = length(rownames(attr(terms(fml1), "factors")))
    t1 = terms(fml1)
  }

  if (p >= 1) {
    mf1 <- model.frame(fml1, data = df, na.action = na.pass)
    mm <- model.matrix(fml1, mf1)
    xlevels = mf1$xlev
    param = colnames(mm)
    colnames(mm) = make.names(colnames(mm))
    varnames = colnames(mm)[-1]
    for (i in 1:length(varnames)) {
      if (!(varnames[i] %in% names(df))) {
        df[,varnames[i]] = mm[,varnames[i]]
      }
    }
  } else {
    xlevels = NULL
    param = "(Intercept)"
    varnames = ""
  }

  fit <- phregcpp(data = df, rep = rep, stratum = stratum, time = time,
                  time2 = time2, event = event, covariates = varnames,
                  weight = weight, offset = offset, id = id,
                  ties = ties, init = init,
                  robust = robust, est_basehaz = est_basehaz,
                  est_resid = est_resid, firth = firth,
                  plci = plci, alpha = alpha,
                  maxiter = maxiter, eps = eps)

  fit$p = fit$sumstat$p[1]

  if (fit$p > 0) {
    fit$param = param[-1]
    fit$beta = fit$parest$beta
    names(fit$beta) = rep(fit$param, length(fit$beta)/fit$p)

    if (fit$p > 1) {
      fit$vbeta = as.matrix(fit$parest[, paste0("vbeta.", 1:fit$p)])
      if (robust) {
        fit$vbeta_naive = as.matrix(fit$parest[, paste0("vbeta_naive.",
                                                        1:fit$p)])
      }
    } else {
      fit$vbeta = as.matrix(fit$parest[, "vbeta"])
      if (robust) {
        fit$vbeta_naive = as.matrix(fit$parest[, "vbeta_naive"])
      }
    }

    dimnames(fit$vbeta) = list(names(fit$beta), fit$param)
    if (robust) {
      dimnames(fit$vbeta_naive) = list(names(fit$beta), fit$param)
    }
  }

  fit$terms = t1
  if (fit$p > 0) fit$xlevels = xlevels

  fit$settings <- list(
    data = data, rep = rep, stratum = stratum, time = time,
    time2 = time2, event = event, covariates = covariates,
    weight = weight, offset = offset, id = id, ties = ties,
    iniy = init, robust = robust, est_basehaz = est_basehaz,
    est_resid = est_resid, firth = firth, plci = plci,
    alpha = alpha, maxiter = maxiter, eps = eps
  )

  class(fit) = "phregr"
  fit
}



#' @title Survival Curve for Proportional Hazards Regression Models
#' @description Obtains the predicted survivor function for a proportional
#' hazards regression model.
#'
#' @param object The output from the \code{phregr} call.
#' @param newdata A data frame with the same variable names as those that
#'   appear in the \code{phregr} call. For right-censored data, one curve is
#'   produced per row to represent a cohort whose covariates correspond to
#'   the values in \code{newdata}. For counting-process data, one curve is
#'   produced per \code{id} in \code{newdata} to present the survival curve
#'   along the path of time-dependent covariates at the observed event
#'   times in the data used to fit \code{phregr}.
#' @param sefit Whether to compute the standard error of the survival
#'   estimates.
#' @param conftype The type of the confidence interval. One of \code{"none"},
#'   \code{"plain"}, \code{"log"}, \code{"log-log"} (the default), or
#'   \code{"arcsin"}. The \code{arcsin} option bases the intervals on
#'   \code{asin(sqrt(surv))}.
#' @param conflev The level of the two-sided confidence interval for
#'   the survival probabilities. Defaults to 0.95.
#'
#' @details
#' If \code{newdata} is not provided and there is no covariate, survival
#' curves based on the \code{basehaz} data frame will be produced.
#'
#' @return A data frame with the following variables:
#'
#' * \code{id}: The id of the subject for counting-process data with
#'   time-dependent covariates.
#'
#' * \code{time}: The observed times in the data used to fit
#'   \code{phregr}.
#'
#' * \code{nrisk}: The number of patients at risk at the time point in the
#'   data used to fit \code{phregr}.
#'
#' * \code{nevent}: The number of patients having event at the time point
#'   in the data used to fit \code{phregr}.
#'
#' * \code{cumhaz}: The cumulative hazard at the time point.
#'
#' * \code{surv}: The estimated survival probability at the time point.
#'
#' * \code{sesurv}: The standard error of the estimated survival probability.
#'
#' * \code{lower}: The lower confidence limit for survival probability.
#'
#' * \code{upper}: The upper confidence limit for survival probability.
#'
#' * \code{conflev}: The level of the two-sided confidence interval.
#'
#' * \code{conftype}: The type of the confidence interval.
#'
#' * \code{covariates}: The values of covariates based on \code{newdata}.
#'
#' * \code{stratum}: The stratum of the subject.
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @references
#' Terry M. Therneau and Patricia M. Grambsch.
#' Modeling Survival Data: Extending the Cox Model.
#' Springer-Verlag, 2000.
#'
#' @examples
#'
#' library(dplyr)
#'
#' # Example 1 with right-censored data
#' fit1 <- phregr(data = rawdata %>% filter(iterationNumber == 1) %>%
#'                  mutate(treat = 1*(treatmentGroup == 1)),
#'                stratum = "stratum",
#'                time = "timeUnderObservation", event = "event",
#'                covariates = "treat")
#'
#' surv1 <- survfit_phregr(fit1,
#'                         newdata = data.frame(
#'                           stratum = as.integer(c(1,1,2,2)),
#'                           treat = c(1,0,1,0)))
#' head(surv1)
#'
#' # Example 2 with counting process data and robust variance estimate
#' fit2 <- phregr(data = heart %>% mutate(rx = as.numeric(transplant) - 1),
#'                time = "start", time2 = "stop", event = "event",
#'                covariates = c("rx", "age"), id = "id", robust = TRUE)
#'
#' surv2 <- survfit_phregr(fit2,
#'                         newdata = data.frame(
#'                           id = c(4,4,11,11),
#'                           age = c(-7.737,-7.737,-0.019,-0.019),
#'                           start = c(0,36,0,26),
#'                           stop = c(36,39,26,153),
#'                           rx = c(0,1,0,1)))
#' head(surv2)
#'
#' @export
survfit_phregr <- function(object, newdata, sefit = TRUE,
                           conftype = "log-log", conflev = 0.95) {

  if (!inherits(object, "phregr"))
    stop("object must be of class 'phregr'");

  p = object$p
  if (p == 0) {
    beta = 0
    vbeta = 0
  } else {
    beta = object$beta
    vbeta = object$vbeta
  }

  basehaz = object$basehaz

  covariates = object$settings$covariates
  stratum = object$settings$stratum
  offset = object$settings$offset
  id = object$settings$id

  if (id != "") {
    tstart = object$settings$time
    tstop = object$settings$time2
  } else {
    tstart = ""
    tstop = ""
  }

  nvar = length(covariates)
  if (missing(covariates) || is.null(covariates) || (nvar == 1 && (
    covariates[1] == "" || tolower(covariates[1]) == "none"))) {
    p3 = 0
  } else {
    fml1 = formula(paste("~", paste(covariates, collapse = "+")))
    p3 = length(rownames(attr(terms(fml1), "factors")))
  }

  if (p >= 1 && p3 >= 1 && !(missing(newdata) || is.null(newdata))) {
    df = newdata
    mf1 <- model.frame(fml1, data = df, na.action = na.pass,
                       xlev = object$xlevels)
    mm <- model.matrix(fml1, mf1)
    colnames(mm) = make.names(colnames(mm))
    varnames = colnames(mm)[-1]
    for (i in 1:length(varnames)) {
      if (!(varnames[i] %in% names(df))) {
        df[,varnames[i]] = mm[,varnames[i]]
      }
    }
  } else {
    beta = NA
    vbeta = NA
    varnames = ""
  }

  if (!is.matrix(vbeta)) vbeta = as.matrix(vbeta)

  if (missing(basehaz) || is.null(basehaz)) {
    stop("basehaz must be provided")
  }

  if (missing(newdata) || is.null(newdata)) {
    if (p > 0) {
      stop("newdata must be provided for Cox models with covariates")
    } else {
      p_stratum = length(stratum);
      if (p_stratum == 1 && (stratum[1] == "" || stratum[1] == "none")) {
        df = data.frame(dummy_x_ = 0)
      } else {
        df = unique(basehaz[, stratum, drop = FALSE])
      }
    }
  }

  survfit_phregcpp(p = p, beta = beta, vbeta = vbeta, basehaz = basehaz,
                   newdata = df, covariates = varnames,
                   stratum = stratum, offset = offset, id = id,
                   tstart = tstart, tstop = tstop, sefit = sefit,
                   conftype = conftype, conflev = conflev)
}



#' @title Residuals for Parametric Regression Models for Failure Time Data
#' @description Obtains the response, martingale, deviance, dfbeta, and
#' likelihood displacement residuals for a parametric regression model
#' for failure time data.
#'
#' @param object The output from the \code{phregr} call.
#' @param type The type of residuals desired, with options including
#'   \code{"response"}, \code{"martingale"}, \code{"deviance"},
#'   \code{"dfbeta"}, \code{"dfbetas"}, \code{"working"}, \code{"ldcase"},
#'   \code{"ldresp"}, \code{"ldshape"}, and \code{"matrix"}.
#' @param collapse Whether to collapse the residuals by \code{id}.
#' @param weighted Whether to compute weighted residuals.
#'
#' @details
#' The algorithms follow the \code{residuals.survreg} function in the
#' \code{survival} package, except for martingale residuals, which
#' are defined only for event or right-censored data for exponential,
#' weibull, lognormal, and loglogistic distributions.
#'
#' @return
#' Either a vector or a matrix of residuals, depending on the specified type:
#'
#' * \code{response} residuals are on the scale of the original data.
#'
#' * \code{martingale} residuals are event indicators minus the cumulative
#'   hazards for event or right-censored data.
#'
#' * \code{working} residuals are on the scale of the linear predictor.
#'
#' * \code{deviance} residuals are on the log-likelihood scale.
#'
#' * \code{dfbeta} residuals are returned as a matrix, where the
#'   \eqn{i}-th row represents the approximate change in the model
#'   coefficients resulting from the inclusion of subject \eqn{i}.
#'
#' * \code{dfbetas} residuals are similar to \code{dfbeta} residuals, but
#'   each column is scaled by the standard deviation of the
#'   corresponding coefficient.
#'
#' * \code{matrix} residuals are a matrix of derivatives of the
#'   log-likelihood function. Let \eqn{L} be the log-likelihood, \eqn{p} be
#'   the linear predictor (\eqn{X\beta}), and \eqn{s} be \eqn{log(\sigma)}.
#'   Then the resulting matrix contains six columns: \eqn{L},
#'   \eqn{\partial L/\partial p}, \eqn{\partial^2 L/\partial p^2},
#'   \eqn{\partial L/\partial s}, \eqn{\partial^2 L/\partial s^2}, and
#'   \eqn{\partial L^2/\partial p\partial s}.
#'
#' * \code{ldcase} residulas are likelihood displacement for case weight
#'   perturbation.
#'
#' * \code{ldresp} residuals are likelihood displacement for response value
#'   perturbation.
#'
#' * \code{ldshape} residuals are likelihood displacement related to the
#'   shape parameter.
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @references
#' Escobar, L. A. and Meeker, W. Q.
#' Assessing influence in regression analysis with censored data.
#' Biometrics 1992; 48:507-528.
#'
#' @examples
#'
#' library(dplyr)
#'
#' fit1 <- liferegr(
#'   data = tobin %>% mutate(time = ifelse(durable>0, durable, NA)),
#'   time = "time", time2 = "durable",
#'   covariates = c("age", "quant"), dist = "normal")
#'
#' resid <- residuals_liferegr(fit1, type = "response")
#' head(resid)
#'
#' @export
residuals_liferegr <- function(
    object, type=c("response", "martingale", "deviance", "dfbeta", "dfbetas",
                   "working", "ldcase", "ldresp", "ldshape", "matrix"),
    collapse=FALSE, weighted=(type %in% c("dfbeta", "dfbetas"))) {

  if (!inherits(object, "liferegr"))
    stop("object must be of class 'liferegr'");

  p = object$p
  beta = object$beta

  data = object$settings$data
  stratum = object$settings$stratum
  time = object$settings$time
  time2 = object$settings$time2
  event = object$settings$event
  covariates = object$settings$covariates
  weight = object$settings$weight
  offset = object$settings$offset
  id = object$settings$id
  dist = object$settings$dist
  param = object$param

  rownames(data) = NULL

  elements = c(stratum, covariates, weight, offset, id)
  elements = unique(elements[elements != "" & elements != "none"])
  if (!(length(elements) == 0)) {
    fml = formula(paste("~", paste(elements, collapse = "+")))
    mf = model.frame(fml, data = data, na.action = na.omit)
  } else {
    mf = model.frame(formula("~1"), data = data)
  }

  rownum = as.integer(rownames(mf))
  df = data[rownum,]

  nvar = length(covariates)
  if (missing(covariates) || is.null(covariates) || (nvar == 1 && (
    covariates[1] == "" || tolower(covariates[1]) == "none"))) {
    p3 = 0
  } else {
    fml1 = formula(paste("~", paste(covariates, collapse = "+")))
    p3 = length(rownames(attr(terms(fml1), "factors")))
  }

  if (p >= 1 && p3 >= 1) {
    mf1 <- model.frame(fml1, data = df, na.action = na.pass)
    mm <- model.matrix(fml1, mf1)
    colnames(mm) = make.names(colnames(mm))
    varnames = colnames(mm)[-1]
    for (i in 1:length(varnames)) {
      if (!(varnames[i] %in% names(df))) {
        df[,varnames[i]] = mm[,varnames[i]]
      }
    }
  } else {
    varnames = ""
  }

  type <- match.arg(type)

  if (type=='dfbeta' || type=='dfbetas') {
    if (missing(weighted))
      weighted <- TRUE  # different default for this case
  }

  n = nrow(data)

  vv <- drop(object$vbeta_naive)
  if (is.null(vv)) vv <- drop(object$vbeta)

  if (weight != "") {
    weights <- df[[weight]]
  } else {
    weights <- rep(1,n)
  }

  if (id != "") {
    idn <- df[[id]]
  } else {
    idn <- seq(1,n)
  }

  rr = residuals_liferegcpp(beta = beta,
                            vbeta = vv,
                            data = df,
                            stratum = stratum,
                            time = time,
                            time2 = time2,
                            event = event,
                            covariates = varnames,
                            weight = weight,
                            offset = offset,
                            id = id,
                            dist = dist,
                            type = type,
                            collapse = collapse,
                            weighted = weighted)

  if (type=="response" || type=="martingale" || type=="deviance" ||
      type=="working" || type=="ldcase" || type=="ldresp" ||
      type=="ldshape") {
    rr <- as.numeric(rr)
  } else if (type=="dfbeta" || type=="dfbetas") {
    colnames(rr) <- param
  } else {
    colnames(rr) <- c("g", "dg", "ddg", "ds", "dds", "dsg");
  }

  rr
}


#' @title Residuals for Proportional Hazards Regression Models
#' @description Obtains the martingale, deviance, score, or Schoenfeld
#' residuals for a proportional hazards regression model.
#'
#' @param object The output from the \code{phregr} call.
#' @param type The type of residuals desired, with options including
#'   \code{"martingale"}, \code{"deviance"}, \code{"score"},
#'   \code{"schoenfeld"}, \code{"dfbeta"}, \code{"dfbetas"}, and
#'   \code{"scaledsch"}.
#' @param collapse Whether to collapse the residuals by \code{id}.
#'   This is not applicable for Schoenfeld type residuals.
#' @param weighted Whether to compute weighted residuals.
#'
#' @details
#' For score and Schoenfeld type residuals, the proportional hazards model
#' must include at least one covariate. The algorithms for \code{deviance},
#' \code{dfbeta}, \code{dfbetas}, and \code{scaledsch} residuals follow
#' the \code{residuals.coxph} function in the \code{survival} package.
#'
#' @return For martingale and deviance residuals, the result is a vector
#' with one element corresponding to each subject (without \code{collapse}).
#' For score residuals, the result is a matrix where each row represents
#' a subject and each column corresponds to a variable. The row order
#' aligns with the input data used in the original fit. For Schoenfeld
#' residuals, the result is a matrix with one row for each event and
#' one column per variable. These rows are sorted by time within strata,
#' with the attributes \code{stratum} and \code{time} included.
#'
#' Score residuals represent each individual's contribution to the score
#' vector. Two commonly used transformations of this are \code{dfbeta},
#' which represents the approximate change in the coefficient vector
#' if the observation is excluded, and \code{dfbetas}, which gives the
#' approximate change in the coefficients scaled by the standard error
#' of the coefficients.
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @references
#' Terry M. Therneau, Patricia M. Grambsch, and Thomas M. Fleming.
#' Martingale based residuals for survival models.
#' Biometrika 1990; 77:147-160.
#'
#' Patricia M. Grambsch and Terry M. Therneau.
#' Proportional hazards tests and diagnostics based on weighted residuals.
#' Biometrika 1994; 81:515-26.
#'
#' @examples
#'
#' library(dplyr)
#'
#' # Example 1 with right-censored data
#' fit1 <- phregr(data = rawdata %>% filter(iterationNumber == 1) %>%
#'                  mutate(treat = 1*(treatmentGroup == 1)),
#'                stratum = "stratum",
#'                time = "timeUnderObservation", event = "event",
#'                covariates = "treat")
#'
#' ressco <- residuals_phregr(fit1, type = "score")
#' head(ressco)
#'
#' # Example 2 with counting process data
#' fit2 <- phregr(data = heart %>% mutate(rx = as.numeric(transplant) - 1),
#'                time = "start", time2 = "stop", event = "event",
#'                covariates = c("rx", "age"), id = "id", robust = TRUE)
#'
#' resssch <- residuals_phregr(fit2, type = "scaledsch")
#' head(resssch)
#'
#' @export
residuals_phregr <- function(
    object, type=c("martingale", "deviance", "score", "schoenfeld",
                   "dfbeta", "dfbetas", "scaledsch"),
    collapse=FALSE, weighted=(type %in% c("dfbeta", "dfbetas"))) {

  if (!inherits(object, "phregr"))
    stop("object must be of class 'phregr'");

  p = object$p
  beta = object$beta
  residuals = object$residuals

  data = object$settings$data
  stratum = object$settings$stratum
  time = object$settings$time
  time2 = object$settings$time2
  event = object$settings$event
  covariates = object$settings$covariates
  weight = object$settings$weight
  offset = object$settings$offset
  id = object$settings$id
  ties = object$settings$ties
  param = object$param

  rownames(data) = NULL

  elements = c(stratum, time, event, covariates, weight, offset, id)
  elements = unique(elements[elements != "" & elements != "none"])
  fml = formula(paste("~", paste(elements, collapse = "+")))
  mf = model.frame(fml, data = data, na.action = na.omit)

  rownum = as.integer(rownames(mf))
  df = data[rownum,]

  nvar = length(covariates)
  if (missing(covariates) || is.null(covariates) || (nvar == 1 && (
    covariates[1] == "" || tolower(covariates[1]) == "none"))) {
    p3 = 0
  } else {
    fml1 = formula(paste("~", paste(covariates, collapse = "+")))
    p3 = length(rownames(attr(terms(fml1), "factors")))
  }

  if (p >= 1 && p3 >= 1) {
    mf1 <- model.frame(fml1, data = df, na.action = na.pass)
    mm <- model.matrix(fml1, mf1)
    colnames(mm) = make.names(colnames(mm))
    varnames = colnames(mm)[-1]
    for (i in 1:length(varnames)) {
      if (!(varnames[i] %in% names(df))) {
        df[,varnames[i]] = mm[,varnames[i]]
      }
    }
  } else {
    varnames = ""
  }


  type <- match.arg(type)

  if (type=='dfbeta' || type=='dfbetas') {
    if (missing(weighted))
      weighted <- TRUE  # different default for this case
  }

  n <- length(residuals)

  vv <- object$vbeta_naive
  if (is.null(vv)) vv <- object$vbeta

  if (weight != "") {
    weights <- df[[weight]]
  } else {
    weights <- rep(1,n)
  }

  if (id != "") {
    idn <- df[[id]]
  } else {
    idn <- seq(1,n)
  }


  if (p == 0) { # null Cox model case
    beta = 0
    vv = matrix(0,1,1)
  }

  temp = residuals_phregcpp(p = p,
                            beta = beta,
                            vbeta = vv,
                            resmart = residuals,
                            data = df,
                            stratum = stratum,
                            time = time,
                            time2 = time2,
                            event = event,
                            covariates = varnames,
                            weight = weight,
                            offset = offset,
                            id = id,
                            ties = ties,
                            type = type,
                            collapse = collapse,
                            weighted = weighted)


  if (type == "martingale" || type == "deviance") {
    rr <- temp$resid
  } else {
    if (p == 1) {
      rr = c(temp$resid)
    } else {
      rr = temp$resid
    }

    if (type == "schoenfeld" || type == "scaledsch") {
      attr(rr, "time") = temp$time
      if (length(temp) == 3) {
        attr(rr, "strata") = temp$strata
      }
    }
  }

  if (is.matrix(rr)) colnames(rr) <- param

  rr
}


#' @title Logistic Regression Models for Binary Data
#' @description Obtains the parameter estimates from logistic regression
#' models with binary data.
#'
#' @param data The input data frame that contains the following variables:
#'
#'   * \code{rep}: The replication for by-group processing.
#'
#'   * \code{event}: The event indicator, 1=event, 0=no event.
#'
#'   * \code{covariates}: The values of baseline covariates.
#'
#'   * \code{freq}: The frequency for each observation.
#'
#'   * \code{weight}: The weight for each observation.
#'
#'   * \code{offset}: The offset for each observation.
#'
#'   * \code{id}: The optional subject ID to group the score residuals
#'     in computing the robust sandwich variance.
#'
#' @param rep The name(s) of the replication variable(s) in the input data.
#' @param event The name of the event variable in the input data.
#' @param covariates The vector of names of baseline covariates
#'   in the input data.
#' @param freq The name of the frequency variable in the input data.
#'   The frequencies must be the same for all observations within each
#'   cluster as indicated by the id. Thus freq is the cluster frequency.
#' @param weight The name of the weight variable in the input data.
#' @param offset The name of the offset variable in the input data.
#' @param id The name of the id variable in the input data.
#' @param link The link function linking the response probabilities to the
#'   linear predictors. Options include "logit" (default), "probit", and
#'   "cloglog" (complementary log-log).
#' @param init A vector of initial values for the model parameters.
#'   By default, initial values are derived from an intercept-only model.
#' @param robust Whether a robust sandwich variance estimate should be
#'   computed. In the presence of the id variable, the score residuals
#'   will be aggregated for each id when computing the robust sandwich
#'   variance estimate.
#' @param firth Whether the firth's bias reducing penalized likelihood
#'   should be used. The default is \code{FALSE}.
#' @param flic Whether to apply intercept correction to obtain more
#'   accurate predicted probabilities. The default is \code{FALSE}.
#' @param plci Whether to obtain profile likelihood confidence interval.
#' @param alpha The two-sided significance level.
#' @param maxiter The maximum number of iterations.
#' @param eps The tolerance to declare convergence.
#'
#' @details
#' Fitting a logistic regression model using Firth's bias reduction method
#' is equivalent to penalization of the log-likelihood by the Jeffreys prior.
#' Firth's penalized log-likelihood is given by
#' \deqn{l(\beta) + \frac{1}{2} \log(\mbox{det}(I(\beta)))}
#' and the components of the gradient \eqn{g(\beta)} are computed as
#' \deqn{g(\beta_j) + \frac{1}{2} \mbox{trace}\left(I(\beta)^{-1}
#' \frac{\partial I(\beta)}{\partial \beta_j}\right)}
#' The Hessian matrix is not modified by this penalty.
#'
#' Firth's method reduces bias in maximum likelihood estimates of
#' coefficients, but it introduces a bias toward one-half in the
#' predicted probabilities.
#'
#' A straightforward modification to Firth’s logistic regression to
#' achieve unbiased average predicted probabilities involves a post hoc
#' adjustment of the intercept. This approach, known as Firth’s logistic
#' regression with intercept correction (FLIC), preserves the
#' bias-corrected effect estimates. By excluding the intercept from
#' penalization, it ensures that we don't sacrifice the accuracy of
#' effect estimates to improve the predictions.
#'
#' @return A list with the following components:
#'
#' * \code{sumstat}: The data frame of summary statistics of model fit
#'   with the following variables:
#'
#'     - \code{n}: The number of subjects.
#'
#'     - \code{nevents}: The number of events.
#'
#'     - \code{loglik0}: The (penalized) log-likelihood under null.
#'
#'     - \code{loglik1}: The maximum (penalized) log-likelihood.
#'
#'     - \code{niter}: The number of Newton-Raphson iterations.
#'
#'     - \code{p}: The number of parameters, including the intercept,
#'       and regression coefficients associated with the covariates.
#'
#'     - \code{link}: The link function.
#'
#'     - \code{robust}: Whether a robust sandwich variance estimate should
#'       be computed.
#'
#'     - \code{firth}: Whether the firth's penalized likelihood is used.
#'
#'     - \code{flic}: Whether to apply intercept correction.
#'
#'     - \code{fail}: Whether the model fails to converge.
#'
#'     - \code{loglik0_unpenalized}: The unpenalized log-likelihood under null.
#'
#'     - \code{loglik1_unpenalized}: The maximum unpenalized log-likelihood.
#'
#'     - \code{rep}: The replication.
#'
#' * \code{parest}: The data frame of parameter estimates with the
#'   following variables:
#'
#'     - \code{param}: The name of the covariate for the parameter estimate.
#'
#'     - \code{beta}: The parameter estimate.
#'
#'     - \code{sebeta}: The standard error of parameter estimate.
#'
#'     - \code{z}: The Wald test statistic for the parameter.
#'
#'     - \code{expbeta}: The exponentiated parameter estimate.
#'
#'     - \code{vbeta}: The covariance matrix for parameter estimates.
#'
#'     - \code{lower}: The lower limit of confidence interval.
#'
#'     - \code{upper}: The upper limit of confidence interval.
#'
#'     - \code{p}: The p-value from the chi-square test.
#'
#'     - \code{method}: The method to compute the confidence interval and
#'       p-value.
#'
#'     - \code{sebeta_naive}: The naive standard error of parameter estimate.
#'
#'     - \code{vbeta_naive}: The naive covariance matrix of parameter
#'       estimates.
#'
#'     - \code{rep}: The replication.
#'
#' * \code{fitted}: The data frame with the following variables:
#'
#'     - \code{linear_predictors}: The linear fit on the link function scale.
#'
#'     - \code{fitted_values}: The fitted probabilities of having an event,
#'       obtained by transforming the linear predictors by the inverse of
#'       the link function.
#'
#'     - \code{rep}: The replication.
#'
#' * \code{p}: The number of parameters.
#'
#' * \code{link}: The link function.
#'
#' * \code{param}: The parameter names.
#'
#' * \code{beta}: The parameter estimate.
#'
#' * \code{vbeta}: The covariance matrix for parameter estimates.
#'
#' * \code{vbeta_naive}: The naive covariance matrix for parameter estimates.
#'
#' * \code{linear_predictors}: The linear fit on the link function scale.
#'
#' * \code{fitted_values}: The fitted probabilities of having an event.
#'
#' * \code{terms}: The terms object.
#'
#' * \code{xlevels}: A record of the levels of the factors used in fitting.
#'
#' * \code{settings}: A list containing the input parameter values.
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @references
#' David Firth.
#' Bias Reduction of Maximum Likelihood Estimates.
#' Biometrika 1993; 80:27–38.
#'
#' Georg Heinze and Michael Schemper.
#' A solution to the problem of separation in logistic regression.
#' Statistics in Medicine 2002;21:2409–2419.
#'
#' Rainer Puhr, Georg Heinze, Mariana Nold, Lara Lusa, and
#' Angelika Geroldinger.
#' Firth's logistic regression with rare events: accurate effect
#' estimates and predictions?
#' Statistics in Medicine 2017; 36:2302-2317.
#'
#' @examples
#'
#' (fit1 <- logisregr(
#'   ingots, event = "NotReady", covariates = "Heat*Soak", freq = "Freq"))
#'
#' @export
logisregr <- function(data, rep = "", event = "event", covariates = "",
                      freq = "", weight = "", offset = "", id = "",
                      link = "logit", init = NA_real_,
                      robust = FALSE, firth = FALSE,
                      flic = FALSE, plci = FALSE, alpha = 0.05,
                      maxiter = 50, eps = 1.0e-9) {

  rownames(data) = NULL

  elements = c(rep, event, covariates, freq, weight, offset)
  elements = unique(elements[elements != "" & elements != "none"])
  fml = formula(paste("~", paste(elements, collapse = "+")))
  mf = model.frame(fml, data = data, na.action = na.omit)

  rownum = as.integer(rownames(mf))
  df = data[rownum,]

  nvar = length(covariates)
  if (missing(covariates) || is.null(covariates) || (nvar == 1 && (
    covariates[1] == "" || tolower(covariates[1]) == "none"))) {
    p = 0
    t1 = terms(formula("~1"))
  } else {
    fml1 = formula(paste("~", paste(covariates, collapse = "+")))
    p = length(rownames(attr(terms(fml1), "factors")))
    t1 = terms(fml1)
  }

  if (p >= 1) {
    mf1 <- model.frame(fml1, data = df, na.action = na.pass)
    mm <- model.matrix(fml1, mf1)
    xlevels = mf1$xlev
    param = colnames(mm)
    colnames(mm) = make.names(colnames(mm))
    varnames = colnames(mm)[-1]
    for (i in 1:length(varnames)) {
      if (!(varnames[i] %in% names(df))) {
        df[,varnames[i]] = mm[,varnames[i]]
      }
    }
  } else {
    xlevels = NULL
    param = "(Intercept)"
    varnames = ""
  }

  fit <- logisregcpp(data = df, rep = rep, event = event,
                     covariates = varnames, freq = freq, weight = weight,
                     offset = offset, id = id, link = link,
                     init = init, robust = robust,
                     firth = firth, flic = flic, plci = plci,
                     alpha = alpha, maxiter = maxiter, eps = eps)

  fit$p <- fit$sumstat$p[1]
  fit$link <- fit$sumstat$link[1]

  if (fit$p > 0) {
    fit$param = param
    fit$beta = fit$parest$beta
    names(fit$beta) = rep(fit$param, length(fit$beta)/fit$p)

    if (fit$p > 1) {
      fit$vbeta = as.matrix(fit$parest[, paste0("vbeta.", 1:fit$p)])
      if (robust) {
        fit$vbeta_naive = as.matrix(fit$parest[, paste0("vbeta_naive.",
                                                        1:fit$p)])
      }
    } else {
      fit$vbeta = as.matrix(fit$parest[, "vbeta"])
      if (robust) {
        fit$vbeta_naive = as.matrix(fit$parest[, "vbeta_naive"])
      }
    }

    dimnames(fit$vbeta) = list(names(fit$beta), fit$param)
    if (robust) {
      dimnames(fit$vbeta_naive) = list(names(fit$beta), fit$param)
    }
  }

  fit$linear_predictors = fit$fitted$linear_predictors
  fit$fitted_values = fit$fitted$fitted_values

  fit$terms = t1
  if (fit$p > 0) fit$xlevels = xlevels

  fit$settings <- list(
    data = data, rep = rep, event = event, covariates = covariates,
    freq = freq, weight = weight, offset = offset, id = id,
    link = link, init = init, robust = robust, firth = firth, flic = flic,
    plci = plci, alpha = alpha, maxiter = maxiter, eps = eps
  )

  class(fit) = "logisregr"
  fit
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


#' @title Hazard Function for Progressive Disease (PD) Given Correlation
#' Between PFS and OS
#'
#' @description
#' Computes the hazard function of a piecewise exponential
#' distribution for progressive disease (PD), such that the
#' resulting hazard function for progression-free survival (PFS)
#' closely matches a given piecewise hazard for PFS
#' and the correlation coefficient between PFS and OS matches
#' the specified value.
#'
#' @inheritParams param_piecewiseSurvivalTime
#' @param hazard_pfs A scalar or numeric vector specifying the
#'   hazard(s) for PFS based on a piecewise exponential distribution.
#' @param hazard_os A scalar or numeric vector specifying the
#'   hazard(s) for overall survival (OS) based on a piecewise
#'   exponential distribution.
#' @param rho_pfs_os A numeric value specifying the correlation
#'   between PFS and OS times.
#'
#' @details
#' This function determines the hazard vector \eqn{\lambda_{\text{pd}}}
#' for the piecewise exponential distribution of PD, so that the
#' implied survival function for PFS time,
#' \eqn{T_{\text{pfs}} = \min(T_{\text{pd}}, T_{\text{os}})}, closely
#' matches the specified piecewise exponential distribution for PFS
#' with hazard vector \eqn{\lambda_{\text{pfs}}}, and the correlation
#' between PFS and OS times matches the specified value \code{rho_pfs_os}.
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
#' * \code{rho_pd_os}: The correlation between PD and OS times.
#'
#' * \code{rho_pfs_os}: The correlation between PFS and OS times (as input).
#'
#' @author
#' Kaifeng Lu (\email{kaifenglu@gmail.com})
#'
#' @examples
#' u <- c(0, 1, 3, 4)
#' lambda1 <- c(0.0151, 0.0403, 0.0501, 0.0558)
#' lambda2 <- 0.0145
#' rho_pfs_os <- 0.5
#' hazard_pd2(u, lambda1, lambda2, rho_pfs_os)
#'
#' @export
hazard_pd2 <- function(piecewiseSurvivalTime = 0,
                       hazard_pfs = 0.0578,
                       hazard_os = 0.02,
                       rho_pfs_os = 0.5) {

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

  if (rho_pfs_os <= -1 || rho_pfs_os >= 1) {
    stop("corr_pfs_os must lie between -1 and 1")
  }

  hazard_pd2cpp(piecewiseSurvivalTime, hazard_pfs, hazard_os, rho_pfs_os)
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


#' @title Random Multivariate Normal Generation
#'
#' @description Generates random samples from a multivariate normal
#' distribution with a specified mean vector and covariance matrix.
#'
#' @param n The number of samples to generate.
#' @param mean A numeric vector representing the mean of the distribution.
#' @param sigma A numeric matrix representing the covariance matrix.
#'
#' @details
#' This function generates samples from a multivariate normal distribution
#' using the Cholesky decomposition method. It first computes the Cholesky
#' factorization of the covariance matrix, then generates standard normal
#' random variables, and finally transforms them to the desired multivariate
#' normal distribution.
#'
#' @return A numeric matrix where each row represents a sample from the
#' multivariate normal distribution.
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @examples
#' # Generate 5 samples from a bivariate normal distribution with mean (0,0)
#' # and covariance matrix [[1, 0.5], [0.5, 1]]
#'
#' set.seed(314159)
#' rmvnorm(5, c(0, 0), matrix(c(1, 0.5, 0.5, 1), nrow=2))
#'
#' @export
rmvnorm <- function(n, mean = rep(0, nrow(sigma)),
                    sigma = diag(length(mean))) {

  if (n <= 0 || n != as.integer(n)) {
    stop("n must be a positive integer")
  }
  if (!is.numeric(mean) || !is.numeric(sigma)) {
    stop("mean and sigma must be numeric vectors and matrices")
  }
  if (!isSymmetric(sigma, tol = sqrt(.Machine$double.eps),
                   check.attributes = FALSE)) {
    stop("sigma must be a symmetric matrix")
  }
  if (length(mean) != nrow(sigma))
    stop("mean and sigma have non-conforming size")

  rmvnormcpp(n, mean, sigma)
}


#' @title Assess Proportional Hazards Assumption Based on Supremum Test
#' @description Obtains the standardized score processes and the simulated
#' distribution under the null hypothesis as well as the p-values for
#' the supremum tests.
#'
#' @param object The output from the \code{phregr} call.
#' @param resample The number of simulation samples for the supremem test.
#' @param seed The random seed for the simulations.
#'
#' @details
#' The supremum test corresponds to the ASSESS statement with \code{ph}
#' option of SAS PROC PHREG.
#'
#' @return A list with the following components:
#'
#' * \code{time} the unique event times.
#'
#' * \code{score_t} the observed standardized score process.
#'
#' * \code{score_t_list} a list of simulated standardized score processes
#'   under the null hypothesis.
#'
#' * \code{max_abs_value} the supremum of the absolute value of the observed
#'   standardized score process for each covariate and the supremum of
#'   the sum of absolute values of the observed standardized score processes
#'   across all covariates.
#'
#' * \code{p_value} the p-values for the supremum tests for each covariate
#'   and the global test.
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @references
#' D. Y. Lin, L. J. Wei, and Z. Ying.
#' Checking the Cox model with cumulative sums of martingale-based
#' residuals.
#' Biometrika 1993; 80:557-572.
#'
#' @examples
#'
#' library(dplyr)
#'
#' fit <- phregr(data = liver, time = "Time", event = "Status",
#'               covariates = c("log(Bilirubin)", "log(Protime)",
#'                              "log(Albumin)", "Age", "Edema"),
#'               ties = "breslow")
#'
#' aph <- assess_phregr(fit, resample = 1000, seed = 314159)
#'
#' aph
#'
#' @export
assess_phregr <- function(object, resample = 1000, seed = 12345) {

  if (!inherits(object, "phregr"))
    stop("object must be of class 'phregr'");

  p = object$p
  beta = object$beta
  vbeta = object$vbeta
  data = object$settings$data
  stratum = object$settings$stratum
  time = object$settings$time
  time2 = object$settings$time2
  event = object$settings$event
  covariates = object$settings$covariates
  weight = object$settings$weight
  offset = object$settings$offset
  ties = object$settings$ties

  rownames(data) = NULL

  elements = c(stratum, time, event, covariates, weight, offset)
  elements = unique(elements[elements != "" & elements != "none"])
  fml = formula(paste("~", paste(elements, collapse = "+")))
  mf = model.frame(fml, data = data, na.action = na.omit)

  rownum = as.integer(rownames(mf))
  df = data[rownum,]

  nvar = length(covariates)
  if (missing(covariates) || is.null(covariates) || (nvar == 1 && (
    covariates[1] == "" || tolower(covariates[1]) == "none"))) {
    p3 = 0
  } else {
    fml1 = formula(paste("~", paste(covariates, collapse = "+")))
    p3 = length(rownames(attr(terms(fml1), "factors")))
  }

  if (p >= 1 && p3 >= 1) {
    mf1 <- model.frame(fml1, data = df, na.action = na.pass)
    mm <- model.matrix(fml1, mf1)
    colnames(mm) = make.names(colnames(mm))
    varnames = colnames(mm)[-1]
    for (i in 1:length(varnames)) {
      if (!(varnames[i] %in% names(df))) {
        df[,varnames[i]] = mm[,varnames[i]]
      }
    }
  } else {
    varnames = ""
  }

  aph <- assess_phregcpp(p = p,
                         beta = beta,
                         vbeta = vbeta,
                         data = df,
                         stratum = stratum,
                         time = time,
                         time2 = time2,
                         event = event,
                         covariates = varnames,
                         weight = weight,
                         offset = offset,
                         ties = ties,
                         resample = resample,
                         seed = seed)

  aph$covariates <- varnames
  aph$resample <- resample
  aph$seed <- seed

  class(aph) <- "assess_phregr"
  aph
}


#' @title Assess Proportional Hazards Assumption Based on Scaled
#' Schoenfeld Residuals
#' @description Obtains the scaled Schoenfeld residuals and tests the
#' proportional hazards assumption using a score test for the interaction
#' between each covariate and a transformed time variable.
#'
#' @param object The output from the \code{phregr} call.
#' @param transform A character string indicating how survival times
#'   should be transformed before the test is performed. Supported values
#'   include "identity", "log", "rank", and "km" (default).
#'   You may also supply the name of a user-defined function that
#'   takes one argument.
#'
#' @details
#' This corresponds to the \code{cox.zph} function from the \code{survival}
#' package with \code{terms = FALSE} and \code{global = TRUE}.
#'
#' @return A list with the following components:
#'
#' * \code{table} A matrix with one row for each parameter and a final
#'   row for the global test. The columns contain the score test
#'   for adding the time-dependent term, the degrees of freedom,
#'   and the two-sided p-value.
#'
#' * \code{x} The transformed time values.
#'
#' * \code{time} The original (untransformed) event times, with tied event
#'   times repeated.
#'
#' * \code{strata} The stratum index for each event.
#'
#' * \code{y} The matrix of scaled Schoenfeld residuals, with one column
#'   for each parameter and one row for each event. Column names correspond
#'   to the parameter names.
#'
#' * \code{var} An approximate covariance matrix of the scaled Schoenfeld
#'   residuals, used to construct an approximate standard error band for
#'   plots.
#'
#' * \code{transform} the transformation applied to the time values.
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @references
#'
#' Patricia M. Grambsch and Terry M. Therneau.
#' Proportional hazards tests and diagnostics based on weighted residuals.
#' Biometrika 1994; 81:515-26.
#'
#' @examples
#'
#' library(dplyr)
#'
#' fit <- phregr(data = liver, time = "Time", event = "Status",
#'               covariates = c("log(Bilirubin)", "log(Protime)",
#'                              "log(Albumin)", "Age", "Edema"),
#'               ties = "breslow")
#'
#' zph <- zph_phregr(fit, transform = "log")
#'
#' zph$table
#'
#' @export
zph_phregr <- function(object, transform = "km") {

  if (!inherits(object, "phregr"))
    stop("object must be of class 'phregr'");

  p = object$p
  beta = object$beta
  vbeta = object$vbeta
  resmart = object$residuals
  data = object$settings$data
  stratum = object$settings$stratum
  time = object$settings$time
  time2 = object$settings$time2
  event = object$settings$event
  covariates = object$settings$covariates
  weight = object$settings$weight
  offset = object$settings$offset
  ties = object$settings$ties

  rownames(data) = NULL

  elements = c(stratum, time, event, covariates, weight, offset)
  elements = unique(elements[elements != "" & elements != "none"])
  fml = formula(paste("~", paste(elements, collapse = "+")))
  mf = model.frame(fml, data = data, na.action = na.omit)

  rownum = as.integer(rownames(mf))
  df = data[rownum,]

  nvar = length(covariates)
  if (missing(covariates) || is.null(covariates) || (nvar == 1 && (
    covariates[1] == "" || tolower(covariates[1]) == "none"))) {
    p3 = 0
  } else {
    fml1 = formula(paste("~", paste(covariates, collapse = "+")))
    p3 = length(rownames(attr(terms(fml1), "factors")))
  }

  if (p >= 1 && p3 >= 1) {
    mf1 <- model.frame(fml1, data = df, na.action = na.pass)
    mm <- model.matrix(fml1, mf1)
    colnames(mm) = make.names(colnames(mm))
    varnames = colnames(mm)[-1]
    for (i in 1:length(varnames)) {
      if (!(varnames[i] %in% names(df))) {
        df[,varnames[i]] = mm[,varnames[i]]
      }
    }
  } else {
    varnames = ""
  }

  zph <- zph_phregcpp(p = p,
                      beta = beta,
                      vbeta = vbeta,
                      resmart = resmart,
                      data = df,
                      stratum = stratum,
                      time = time,
                      time2 = time2,
                      event = event,
                      covariates = varnames,
                      weight = weight,
                      offset = offset,
                      ties = ties,
                      transform = transform)

  class(zph) <- "cox.zph"
  zph
}

