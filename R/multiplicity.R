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
#' @param nthreads The number of threads to use in simulations (0 means
#'   the default RcppParallel behavior).
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
                           spendingTime = NULL,
                           nthreads = 0) {

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

  # Respect user-requested number of threads (best effort)
  if (nthreads > 0) {
    n_physical_cores <- parallel::detectCores(logical = FALSE)
    RcppParallel::setThreadOptions(min(nthreads, n_physical_cores))
  }

  repp1 = repeatedPValuecpp(
    kMax = kMax, typeAlphaSpending = typeAlphaSpending,
    parameterAlphaSpending = parameterAlphaSpending,
    maxInformation = maxInformation, p = p1,
    information = information1, spendingTime = spendingTime1)

  if (nrow(repp1) == 1) { # convert the result to a vector
    repp = as.vector(repp1)
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
#' @param typeAlphaSpending The vector of alpha spending functions for
#'   the hypotheses. Each element is one of the following:
#'   "sfOF" for O'Brien-Fleming type spending function,
#'   "sfP" for Pocock type spending function,
#'   "sfKD" for Kim & DeMets spending function,
#'   "sfHSD" for Hwang, Shi & DeCani spending function.
#'   Defaults to "sfOF" if not provided.
#' @param parameterAlphaSpending The vector of parameter values for the
#'   alpha spending functions for the hypotheses. Each element corresponds
#'   to the value of rho for "sfKD" or gamma for "sfHSD".
#'   Defaults to missing if not provided.
#' @param maxInformation The vector of target maximum information for each
#'   hypothesis. Defaults to a vector of 1s if not provided.
#' @param incidenceMatrix The kMax x m incidence matrix indicating whether the
#'   specific hypothesis will be tested at the given look.
#'   If not provided, defaults to testing each hypothesis at all study looks.
#' @param k1 The number of study looks at the interim analysis.
#' @param p The matrix of raw p-values for each hypothesis by study look.
#' @param information The matrix of observed information for each hypothesis
#'   by study look.
#' @param spendingTime The spending time for alpha spending by study look.
#'   If not provided, it is the same as \code{informationRates} calculated
#'   from \code{information} and \code{maxInformation}.
#' @param nthreads The number of threads to use in simulations (0 means
#'   the default RcppParallel behavior).
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
#'   k1 = 2,
#'   p = matrix(c(0.0062, 0.017, 0.009, 0.13,
#'                0.0002, 0.0035, 0.002, 0.06),
#'              nrow=2, ncol=4, byrow=TRUE),
#'   information = matrix(c(rep(1/3, 4), rep(2/3, 4)),
#'                        nrow=2, ncol=4, byrow=TRUE))
#'
#'
#' @export
fseqbon <- function(w, G, alpha = 0.025, kMax,
                    typeAlphaSpending = NULL,
                    parameterAlphaSpending = NULL,
                    maxInformation = NULL,
                    incidenceMatrix = NULL,
                    k1, p, information,
                    spendingTime = NULL,
                    nthreads = 0) {
  m = length(w)

  if (is.null(typeAlphaSpending)) {
    typeAlphaSpending = rep("sfOF", m)
  }

  if (is.null(parameterAlphaSpending)) {
    parameterAlphaSpending = rep(NA, m)
  }

  if (is.null(spendingTime)) {
    spendingTime = matrix(0, 1, 1)
  }

  if (is.null(maxInformation)) {
    maxInformation = rep(1, m)
  }

  if (is.null(incidenceMatrix)) {
    incidenceMatrix = matrix(1, nrow=kMax, ncol=m)
  }

  # Respect user-requested number of threads (best effort)
  if (nthreads > 0) {
    n_physical_cores <- parallel::detectCores(logical = FALSE)
    RcppParallel::setThreadOptions(min(nthreads, n_physical_cores))
  }

  reject1 <- fseqboncpp(
    w = w, G = G, alpha = alpha, kMax = kMax,
    typeAlphaSpending = typeAlphaSpending,
    parameterAlphaSpending = parameterAlphaSpending,
    maxInformation = maxInformation,
    incidenceMatrix = incidenceMatrix,
    k1 = k1, p = p, information = information,
    spendingTime = spendingTime)

  if (nrow(reject1) == 1) { # convert the result to a vector
    reject = as.vector(reject1)
  } else {
    reject = reject1
  }

  reject
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
  n = ncol(p)

  if (n %% 2 != 0) {
    stop("p must have an even number of columns")
  }

  if (!all(gamma >= 0 & gamma <= 1)) {
    stop("gamma must lie between 0 and 1");
  }

  if (n == 2) {
    gamma = 1
  } else if (length(gamma) == n/2 - 1) {
    gamma = c(gamma, 1)
  } else if (length(gamma) != n/2) {
    stop("The number of families must be half of the number of hypotheses")
  } else {
    gamma[n/2] = 1
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
fstdmix <- function(p, family = NULL, serial, parallel = NULL,
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

  if (is.null(parallel)) {
    parallel = matrix(0, m, m)
  } else if (!is.matrix(parallel)) {
    parallel = matrix(parallel, nrow=m, ncol=m)
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
fmodmix <- function(p, family = NULL, serial, parallel = NULL,
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

  if (is.null(parallel)) {
    parallel = matrix(0, m, m)
  } else if (!is.matrix(parallel)) {
    parallel = matrix(parallel, nrow=m, ncol=m)
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
