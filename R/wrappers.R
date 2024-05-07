#' @title Error spending
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
#' @return A vector of errors spent up to the interim look.
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @examples
#' errorSpent(t = 0.5, error = 0.025, sf = "sfOF")
#' errorSpent(t = c(0.5, 0.75, 1), error = 0.025, sf = "sfHSD", sfpar = -4)
#'
#' @export
errorSpent <- function(t, error, sf = "sfOF", sfpar = NA) {
  sapply(t, errorSpentcpp, error, sf, sfpar)
}


#' @title Distribution function of truncated piecewise exponential
#' distribution
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

  ptpwexpcpp(q, piecewiseSurvivalTime, lambda, lowerBound, lower.tail, log.p)
}


#' @title Quantile function of truncated piecewise exponential distribution
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

  qtpwexpcpp(p, piecewiseSurvivalTime, lambda, lowerBound, lower.tail, log.p)
}


#' @title Random number generation function of truncated piecewise
#' exponential distribution
#' @description Obtains random samples from a truncated piecewise
#' exponential distribution.
#'
#' @param n The number of observations.
#' @inheritParams param_piecewiseSurvivalTime
#' @inheritParams param_lambda
#' @param lowerBound The left truncation time point for the survival time.
#'   Defaults to 0 for no truncation.
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

  rtpwexpcpp(n, piecewiseSurvivalTime, lambda, lowerBound)
}


#' @title Stagewise exit probabilities
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
  exitprobcpp(b, a, theta, I)
}


#' @title Adjusted p-values for Bonferroni-based graphical approaches
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

  x = fadjpboncpp(w, G, p)
  if (nrow(x) == 1) {
    x = as.vector(x)
  }
  x
}


#' @title Adjusted p-values for Dunnett-based graphical approaches
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


#' @title Adjusted p-values for Simes-based graphical approaches
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

  x = fadjpsimcpp(wgtmat, p, family)
  if (nrow(x) == 1) {
    x = as.vector(x)
  }
  x
}


#' @title Adjusted p-values for Holm, Hochberg, and Hommel procedures
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

  x = ftrunccpp(p, test, gamma)
  if (nrow(x) == 1) {
    x = as.vector(x)
  }
  x
}


#' @title Repeated p-values for group sequential design
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
#'   increasing and less than or equal to 1. Defaults to NULL,
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

  repp1 = repeatedPValuecpp(kMax, typeAlphaSpending, parameterAlphaSpending,
                            maxInformation, p1, information1, spendingTime1)

  if (is.vector(p)) { # convert the result to a vector
    repp = c(repp1)
  } else {
    repp = repp1
  }

  repp
}


#' @title Group sequential trials using Bonferroni-based graphical
#' approaches
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
# [[Rcpp::export]]
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

  fseqboncpp(w, G, alpha, kMax, typeAlphaSpending,
             parameterAlphaSpending,
             incidenceMatrix, maxInformation,
             p, information, spendingTime)
}


#' @title Adjusted p-values for stepwise testing procedures for two sequences
#' @description Obtains the adjusted p-values for the stepwise gatekeeping
#' procedures for multiplicity problems involving two sequences of
#' hypotheses.
#'
#' @param p The raw p-values for elementary hypotheses.
#' @param gamma The truncation parameters for each family. The truncation
#'   parameter for the last family is automatically set to 1.
#' @param test The component multiple testing procedure. It is either "Holm"
#'   or "Hochberg", and it defaults to "Hochberg".
#' @param retest Whether to allow retesting. It defaults to TRUE.
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

  x = fstp2seqcpp(p, gamma, test, retest)
  if (nrow(x) == 1) {
    x = as.vector(x)
  }
  x
}

#' @title Adjusted p-values for standard mixture gatekeeping procedures
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
#'   for the last family with active hypotheses. It defaults to TRUE.
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
#' fstdmix(p, family, serial, parallel, gamma, test = "hommel", exhaust = 0)
#'
#' @export
fstdmix <- function(p, family = NULL, serial, parallel,
                    gamma, test = "hommel", exhaust = 1) {
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

  x = fstdmixcpp(p, family, serial, parallel, gamma, test, exhaust)
  if (nrow(x) == 1) {
    x = as.vector(x)
  }
  x
}


#' @title Adjusted p-values for modified mixture gatekeeping procedures
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
#'   for the last family with active hypotheses. It defaults to TRUE.
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
#' fmodmix(p, family, serial, parallel, gamma, test = "hommel", exhaust = 1)
#'
#' @export
fmodmix <- function(p, family = NULL, serial, parallel,
                    gamma, test = "hommel", exhaust = 1) {
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

  x = fmodmixcpp(p, family, serial, parallel, gamma, test, exhaust)
  if (nrow(x) == 1) {
    x = as.vector(x)
  }
  x
}



#' @title Efficacy boundaries for group sequential design
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
  getBoundcpp(k, informationRates, alpha,
              typeAlphaSpending, parameterAlphaSpending,
              userAlphaSpending, spendingTime,
              efficacyStopping)
}
