#' @title Error spending functions
#' @description Obtains the error spent at the given information fractions
#' for the specified error spending function.
#'
#' @param t A vector of information fractions.
#' @param error Total error to spend.
#' @param sf Spending function. One of the following: "sfOF" for
#'  O'Brien-Fleming type spending function, "sfP" for Pocock type spending
#'  function, "sfKD" for Kim & DeMets spending function, and "sfHSD" for
#'  Hwang, Shi & DeCani spending function. Defaults to "sfOF".
#' @param sfpar Parameter for the spending function. Corresponds to rho for
#'  "sfKD" and gamma for "sfHSD".
#'
#' @return A vector of errors spent up to the interim look.
#'
#' @examples
#' errorSpent(t = 0.5, error = 0.025, sf = "sfOF")
#' errorSpent(t = c(0.5, 0.75, 1), error = 0.025, sf = "sfHSD", sfpar = -4)
#'
#' @export
errorSpent <- function(t, error, sf = "sfOF", sfpar = NA) {
  sapply(t, errorSpentcpp, error, sf, sfpar)
}



#' @title Quantile function of truncated piecewise exponential distribution
#' @description Obtains the quantile of a piecewise expoenential distribution
#' given that it exceeds a specified lower bound.
#'
#' @param probability The scalar probability corresponding to the quantile.
#' @inheritParams param_piecewiseSurvivalTime
#' @inheritParams param_lambda
#' @param lowerBound The left truncation time point for the survival time.
#' Defaults to 0 for no truncation.
#'
#' @return The quantile x such that
#' P(X > x | X > lowerBound) = 1 - probability.
#'
#' @examples
#' qtpwexp(probability = c(0.3, 0.5), piecewiseSurvivalTime = c(0, 6, 9, 15),
#'         lambda = c(0.025, 0.04, 0.015, 0.007))
#'
#' @export
qtpwexp <- function(probability, piecewiseSurvivalTime = 0, 
                    lambda = 0.0578, lowerBound = 0) {
  sapply(probability, qtpwexpcpp, piecewiseSurvivalTime, lambda, lowerBound)
}



#' @title Stagewise exit probabilities
#' @description Obtains the stagewise exit probabilities for both efficacy and
#' futility stopping.
#'
#' @param b Upper boundaries on the z-test statistic scale.
#' @param a Lower boundaries on the z-test statistic scale. Defaults to
#' \code{c(rep(-6.0, kMax-1), b[kMax])} if left unspecified, where
#' \code{kMax = length(b)}.
#' @param theta Stagewise parameter of interest, e.g., \code{-U/V} for
#'   weighted log-rank test, where \code{U} is the mean and \code{V} is
#'   the variance of the weighted log-rank test score statistic at each stage.
#'   For proportional hazards and conventional log-rank test, use the
#'   scalar input, \code{theta = -log(HR)}. Defaults to 0 corresponding to 
#'   the null hypothesis.
#' @param I Stagewise cumulative information, e.g., \code{V}, the variance
#'   of the weighted log-rank test score statistic at each stage. For
#'   conventional log-rank test, information can be approximated by
#'   \code{phi*(1-phi)*D}, where \code{phi} is the probability of being
#'   allocated to the active arm, and \code{D} is the total number of events 
#'   at each stage. Defaults to \code{seq(1, kMax)} if left unspecified. 
#'
#' @return A list of stagewise exit probabilities: one vector for efficacy
#' stopping probabilities, and the other vector for futility stopping
#' probabilities.
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
#' procedures. Statistics in Medicine. 2009;28:586-604. 
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
#' parameter tests. Biometrical Journal. 2011;53:894-913.
#' 
#' Kaifeng Lu. Graphical approaches using a Bonferroni mixture of weighted 
#' Simes tests. Statistics in Medicine. 2016;35:4041-4055.
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
fadjpsim <- function(wgtmat, p, family) {
  m = ncol(wgtmat)
  
  if (!is.matrix(p)) {
    p = matrix(p, ncol=m)
  }
  
  x = fadjpsimcpp(wgtmat, p, family)
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
#'                p = matrix(c(0.0062, 0.017, 0.009, 0.13,
#'                             0.0002, 0.0035, 0.002, 0.06), 
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
#' Research. 2013;5:311-320. \doi{10.1080/19466315.2013.807748}
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
