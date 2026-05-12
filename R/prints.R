#' @title Print Group Sequential Design
#' @description Prints the stopping boundaries and information inflation
#' factor for group sequential design.
#'
#' @param x The design object to print.
#' @param ... Ensures that all arguments starting from "..." are named.
#'
#' @return A tabular printout of the design elements.
#'
#' @keywords internal
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @export
print.design <- function(x, ...) {
  a <- x$overallResults
  s <- x$byStageResults
  k <- a$kMax

  if (k>1) {
    str1 <- paste0("Group-sequential design with ", k, " stages")
  } else {
    str1 <- "Fixed design"
  }

  str2 <- paste0("theta: ", round(a$theta, 3), ", ",
                 "maximum information: ", round(a$information, 2))

  str3 <- paste0("Overall power: ",
                 round(a$overallReject, 4), ", ",
                 "overall alpha (1-sided): ",
                 round(a$alpha, 4))

  if (x$settings$typeBetaSpending != 'none' ||
      (k > 1 && any(x$byStageResults$futilityBounds[1:(k-1)] > -8))) {
    str3 <- paste0(str3, ", ",
                   "attained alpha: ", round(a$attainedAlpha, 4))
  }

  str4 <- paste0("Drift parameter: ", round(a$drift, 3), ", ",
                 "inflation factor: ", round(a$inflationFactor, 3))

  if (k>1) {
    str5 <- paste0("Expected information under H1: ",
                   round(a$expectedInformationH1, 2), ", ",
                   "expected information under H0: ",
                   round(a$expectedInformationH0, 2))

    asf <- tolower(x$settings$typeAlphaSpending)
    asfpar <- round(x$settings$parameterAlphaSpending, 3)
    asfuser <- round(x$settings$userAlphaSpending, 4)

    bsf <- tolower(x$settings$typeBetaSpending)
    bsfpar <- round(x$settings$parameterBetaSpending, 3)
    bsfuser <- round(x$settings$userBetaSpending, 4)

    if (asf == "of") {
      str6 <- paste0("Alpha spending: O'Brien-Fleming")
    } else if (asf == "p") {
      str6 <- paste0("Alpha spending: Pocock")
    } else if (asf == "wt") {
      str6 <- paste0("Alpha spending: Wang-Tsiatis(Delta = ", asfpar, ")")
    } else if (asf == "sfof") {
      str6 <- paste0("Alpha spending: Lan-DeMets O'Brien-Fleming")
    } else if (asf == "sfp") {
      str6 <- paste0("Alpha spending: Lan-DeMets Pocock")
    } else if (asf == "sfkd") {
      str6 <- paste0("Alpha spending: KD(rho = ", asfpar, ")")
    } else if (asf == "sfhsd") {
      str6 <- paste0("Alpha spending: HSD(gamma = ", asfpar, ")")
    } else if (asf == "user") {
      str6 <- paste0("Alpha spending: User defined(",
                     paste(asfuser, collapse = ","), ")")
    } else {
      str6 <- "Alpha spending: None"
    }

    if (bsf == "of") {
      str7 <- paste0("beta spending: O'Brien-Fleming")
    } else if (bsf == "p") {
      str7 <- paste0("beta spending: Pocock")
    } else if (bsf == "wt") {
      str7 <- paste0("beta spending: Wang-Tsiatis(Delta = ", bsfpar, ")")
    } else if (bsf == "sfof") {
      str7 <- paste0("beta spending: Lan-DeMets O'Brien-Fleming")
    } else if (bsf == "sfp") {
      str7 <- paste0("beta spending: Lan-DeMets Pocock")
    } else if (bsf == "sfkd") {
      str7 <- paste0("beta spending: KD(rho = ", bsfpar, ")")
    } else if (bsf == "sfhsd") {
      str7 <- paste0("beta spending: HSD(gamma = ", bsfpar, ")")
    } else if (bsf == "user") {
      str7 <- paste0("beta spending: User defined(",
                     paste(bsfuser, collapse = ","), ")")
    } else {
      str7 <- "beta spending: None"
    }

    if (!any(is.na(x$settings$spendingTime)) &&
        !all(x$settings$spendingTime == s$informationRates)) {
      str8 <- paste0("Spending time: ",
                     paste(round(x$settings$spendingTime, 3), collapse = ","))
      df1 <- data.frame(x = rep("", 8))
      colnames(df1) <- NULL
      rownames(df1) <- c(str1, str2, str3, str4, str5,
                         paste(str6, str7, sep = ", "), str8, "")
    } else {
      df1 <- data.frame(x = rep("", 7))
      colnames(df1) <- NULL
      rownames(df1) <- c(str1, str2, str3, str4, str5,
                         paste(str6, str7, sep = ", "), "")
    }
  } else {
    df1 <- data.frame(x = rep("", 5))
    colnames(df1) <- NULL
    rownames(df1) <- c(str1, str2, str3, str4, "")
  }


  if (k>1) {
    b <- s[, c("informationRates", "efficacyBounds", "futilityBounds",
               "cumulativeRejection", "cumulativeFutility",
               "cumulativeAlphaSpent", "efficacyTheta", "futilityTheta",
               "efficacyP", "futilityP", "information",
               "cumulativeRejectionH0", "cumulativeFutilityH0")]

    # format number of digits after decimal for each column
    j2 <- 11
    j3 <- c(1,2,3,7,8)
    j4 <- c(4,5,6,9,10,12,13)

    b[j2] <- lapply(b[j2], formatC, format = "f", digits = 2)
    b[j3] <- lapply(b[j3], formatC, format = "f", digits = 3)
    b[j4] <- lapply(b[j4], formatC, format = "f", digits = 4)

    if (x$settings$typeBetaSpending != 'none' ||
        (k > 1 && any(x$byStageResults$futilityBounds[1:(k-1)] > -8))) {
      df <- t(b)
      rownames(df) <- c("Information rate",
                        "Efficacy boundary (Z)",
                        "Futility boundary (Z)",
                        "Cumulative rejection",
                        "Cumulative futility",
                        "Cumulative alpha spent",
                        "Efficacy boundary (theta)",
                        "Futility boundary (theta)",
                        "Efficacy boundary (p)",
                        "Futility boundary (p)",
                        "Information",
                        "Cumulative rejection under H0",
                        "Cumulative futility under H0")

    } else {
      df <- t(b[, c(1,2,4,6,7,9,11)])
      rownames(df) <- c("Information rate",
                        "Efficacy boundary (Z)",
                        "Cumulative rejection",
                        "Cumulative alpha spent",
                        "Efficacy boundary (theta)",
                        "Efficacy boundary (p)",
                        "Information")
    }

    colnames(df) <- paste("Stage", seq_len(ncol(df)), sep=" ")
  } else {
    b <- s[, c("efficacyBounds", "efficacyTheta", "efficacyP")]

    # format number of digits after decimal for each column
    j3 <- c(1,2)
    j4 <- 3

    b[j3] <- lapply(b[j3], formatC, format = "f", digits = 3)
    b[j4] <- lapply(b[j4], formatC, format = "f", digits = 4)

    df <- t(b)

    rownames(df) <- c("Efficacy boundary (Z)",
                      "Efficacy boundary (theta)",
                      "Efficacy boundary (p)")
    colnames(df) <- NA
  }

  print(df1, ..., na.print = "" , quote = FALSE )
  print(df, ..., na.print = "" , quote = FALSE )
  invisible(x)
}


#' @title Print Group Sequential Equivalence Design
#' @description Prints the stopping boundaries for group sequential
#' equivalence design.
#'
#' @param x The design object to print.
#' @param ... Ensures that all arguments starting from "..." are named.
#'
#' @return A tabular printout of the design elements.
#'
#' @keywords internal
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @export
print.designEquiv <- function(x, ...) {
  a <- x$overallResults
  s <- x$byStageResults
  k <- a$kMax

  if (k>1) {
    str1 <- paste0("Group-sequential design with ", k, " stages")
  } else {
    str1 <- "Fixed design"
  }

  str1 <- paste(str1, "for equivalence test")

  str2 <- paste0("Lower equivalence limit: ", round(a$thetaLower, 3), ", ",
                 "upper equivalence limit: ", round(a$thetaUpper, 3), ", ",
                 "parameter value: ", round(a$theta, 3))

  str3 <- paste0("Maximum information: ", round(a$information, 2))

  str4 <- paste0("Overall power: ",
                 round(a$overallReject, 4), ", ",
                 "overall alpha: ",
                 round(a$alpha, 4), ", ",
                 "attained under H10: ",
                 round(a$attainedAlphaH10, 4), ", ",
                 "under H20: ",
                 round(a$attainedAlphaH20, 4))

  if (k>1) {
    str5 <- paste0("Expected information under H1: ",
                   round(a$expectedInformationH1, 2), ", ",
                   "under H10: ",
                   round(a$expectedInformationH10, 2), ", ",
                   "under H20: ",
                   round(a$expectedInformationH20, 2))

    asf <- tolower(x$settings$typeAlphaSpending)
    asfpar <- round(x$settings$parameterAlphaSpending, 3)
    asfuser <- round(x$settings$userAlphaSpending, 4)

    if (asf == "of") {
      str6 <- paste0("Alpha spending: O'Brien-Fleming")
    } else if (asf == "p") {
      str6 <- paste0("Alpha spending: Pocock")
    } else if (asf == "wt") {
      str6 <- paste0("Alpha spending: Wang-Tsiatis(Delta = ", asfpar, ")")
    } else if (asf == "sfof") {
      str6 <- paste0("Alpha spending: Lan-DeMets O'Brien-Fleming")
    } else if (asf == "sfp") {
      str6 <- paste0("Alpha spending: Lan-DeMets Pocock")
    } else if (asf == "sfkd") {
      str6 <- paste0("Alpha spending: KD(rho = ", asfpar, ")")
    } else if (asf == "sfhsd") {
      str6 <- paste0("Alpha spending: HSD(gamma = ", asfpar, ")")
    } else if (asf == "user") {
      str6 <- paste0("Alpha spending: User defined(",
                     paste(asfuser, collapse = ","), ")")
    } else {
      str6 <- "Alpha spending: None"
    }

    if (!any(is.na(x$settings$spendingTime)) &&
        !all(x$settings$spendingTime == s$informationRates)) {
      str7 <- paste0("Spending time: ",
                     paste(round(x$settings$spendingTime, 3), collapse = ","))
      df1 <- data.frame(x = rep("", 8))
      colnames(df1) <- NULL
      rownames(df1) <- c(str1, str2, str3, str4, str5, str6, str7, "")
    } else {
      df1 <- data.frame(x = rep("", 7))
      colnames(df1) <- NULL
      rownames(df1) <- c(str1, str2, str3, str4, str5, str6, "")
    }
  } else {
    df1 <- data.frame(x = rep("", 5))
    colnames(df1) <- NULL
    rownames(df1) <- c(str1, str2, str3, str4, "")
  }


  if (k>1) {
    b <- s[, c("informationRates", "efficacyBounds",
               "cumulativeRejection", "cumulativeAlphaSpent",
               "cumulativeAttainedAlphaH10", "cumulativeAttainedAlphaH20",
               "efficacyThetaLower", "efficacyThetaUpper",
               "efficacyP", "information")]

    # format number of digits after decimal for each column
    j2 <- 10
    j3 <- c(1,2,7,8)
    j4 <- c(3,4,5,6,9)

    b[j2] <- lapply(b[j2], formatC, format = "f", digits = 2)
    b[j3] <- lapply(b[j3], formatC, format = "f", digits = 3)
    b[j4] <- lapply(b[j4], formatC, format = "f", digits = 4)

    df <- t(b)
    rownames(df) <- c("Information rate",
                      "Boundary for each 1-sided test (Z)",
                      "Cumulative rejection",
                      "Cumulative alpha for each 1-sided test",
                      "Cumulative alpha attained under H10",
                      "Cumulative alpha attained under H20",
                      "Boundary for lower limit (theta)",
                      "Boundary for upper limit (theta)",
                      "Boundary for each 1-sided test (p)",
                      "Information")

    colnames(df) <- paste("Stage", seq_len(ncol(df)), sep=" ")
  } else {
    b <- s[, c("efficacyBounds", "efficacyThetaLower",
               "efficacyThetaUpper", "efficacyP")]

    # format number of digits after decimal for each column
    j3 <- c(1,2,3)
    j4 <- 4

    b[j3] <- lapply(b[j3], formatC, format = "f", digits = 3)
    b[j4] <- lapply(b[j4], formatC, format = "f", digits = 4)

    df <- t(b)

    rownames(df) <- c("Boundary for each 1-sided test (Z)",
                      "Boundary for lower limit (theta)",
                      "Boundary for upper limit (theta)",
                      "Boundary for each 1-sided test (p)")
    colnames(df) <- NA
  }

  print(df1, ..., na.print = "" , quote = FALSE )
  print(df, ..., na.print = "" , quote = FALSE )
  invisible(x)
}


#' @title Print Adaptive Group Sequential Design
#' @description Prints the primary and second trial information for
#' an adaptive group sequential design.
#'
#' @param x The adaptDesign object to print.
#' @param ... Ensures that all arguments starting from "..." are named.
#'
#' @return A tabular printout of the design elements.
#'
#' @keywords internal
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @export
print.adaptDesign <- function(x, ...) {
  des1 <- x$primaryTrial

  str1 <- "Primary trial:"
  str2 <- paste0("Group-sequential design with ", des1$kMax, " stages")
  str3 <- paste0("Max information: ", round(des1$maxInformation, 2))
  str4 <- paste0("Interim adaptation look: ",  des1$L, ", ",
                 "z-statistic value: ", round(des1$zL, 3))
  str5 <- paste0("theta: ", round(des1$theta, 3))
  str6 <- paste0("Conditional type I error: ",
                 round(des1$conditionalAlpha, 4))
  str7 <- paste0("Conditional power: ", round(des1$conditionalPower, 3),
                 ", ", "predictive power: ",
                 round(des1$predictivePower, 4))
  str8 <- paste0("Muller & Schafer method for secondary trial: ",
                 des1$MullerSchafer)

  df1a <- data.frame(x = rep("", 9))
  colnames(df1a) <- NULL
  rownames(df1a) <- c(str1, str2, str3, str4, str5, str6, str7, str8, "")

  b <- data.frame(informationRates = des1$informationRates,
                  efficacyBounds = des1$efficacyBounds,
                  futilityBounds = des1$futilityBounds,
                  information = des1$information)

  b[1:3] <- lapply(b[1:3], formatC, format = "f", digits = 3)
  b[4] <- lapply(b[4], formatC, format = "f", digits = 2)

  if (des1$kMax > 1 && any(des1$futilityBounds[1:(des1$kMax-1)] > -8)) {
    df1b <- t(b)
    rownames(df1b) <- c("Information rate",
                        "Efficacy boundary (Z)",
                        "Futility boundary (Z)",
                        "Information")
  } else {
    df1b <- t(b[, c(1,2,4)])
    rownames(df1b) <- c("Information rate",
                        "Efficacy boundary (Z)",
                        "Information")
  }
  colnames(df1b) <- paste("Stage", seq_len(ncol(df1b)), sep=" ")


  des2 <- x$secondaryTrial

  str1 <- "Secondary trial:"

  if (des2$kMax > 1) {
    str2 <- paste0("Group-sequential design with ", des2$kMax, " stages")
  } else {
    str2 <- "Fixed design"
  }

  str3 <- paste0("Maximum information: ", round(des2$maxInformation, 2))
  str4 <- paste0("Overall power: ",
                 round(des2$overallReject, 4), ", ",
                 "overall significance level (1-sided): ",
                 round(des2$alpha, 4))

  df2a <- data.frame(x = rep("", 5))
  colnames(df2a) <- NULL
  rownames(df2a) <- c(str1, str2, str3, str4, "")

  b <- data.frame(informationRates = des2$informationRates,
                  efficacyBounds = des2$efficacyBounds,
                  futilityBounds = des2$futilityBounds,
                  cumulativeRejection = des2$cumulativeRejection,
                  cumulativeFutility = des2$cumulativeFutility,
                  cumulativeAlphaSpent = des2$cumulativeAlphaSpent,
                  information = des2$information)

  # format number of digits after decimal for each column
  j2 <- 7
  j3 <- c(1,2,3)
  j4 <- c(4,5,6)

  b[j2] <- lapply(b[j2], formatC, format = "f", digits = 2)
  b[j3] <- lapply(b[j3], formatC, format = "f", digits = 3)
  b[j4] <- lapply(b[j4], formatC, format = "f", digits = 4)

  if (des2$kMax > 1 && any(des2$futilityBounds[1:(des2$kMax-1)] > -8)) {
    df2b <- t(b)
    rownames(df2b) <- c("Information rate",
                        "Efficacy boundary (Z)",
                        "Futility boundary (Z)",
                        "Cumulative rejection",
                        "Cumulative futility",
                        "Cumulative alpha spent",
                        "Information")

  } else {
    df2b <- t(b[, c(1,2,4,6,7)])
    rownames(df2b) <- c("Information rate",
                        "Efficacy boundary (Z)",
                        "Cumulative rejection",
                        "Cumulative alpha spent",
                        "Information")
  }
  colnames(df2b) <- paste("Stage", seq_len(ncol(df2b)), sep=" ")

  des3 <- x$integratedTrial
  str1 <- "Integrated trial:"
  str2 <- paste0("Group-sequential design with ", des3$kMax, " stages")
  str3 <- paste0("Maximum information: ", round(des3$maxInformation, 2))
  str4 <- paste0("Interim adaptation look: ",  des3$L, ", ",
                 "z-statistic value: ", round(des3$zL, 3))
  df3a <- data.frame(x = rep("", 5))
  colnames(df3a) <- NULL
  rownames(df3a) <- c(str1, str2, str3, str4, "")

  b <- data.frame(informationRates = des3$informationRates,
                  efficacyBounds = des3$efficacyBounds,
                  futilityBounds = des3$futilityBounds,
                  information = des3$information)

  b[1:3] <- lapply(b[1:3], formatC, format = "f", digits = 3)
  b[4] <- lapply(b[4], formatC, format = "f", digits = 2)

  if (des3$kMax > 1 && any(des3$futilityBounds[1:(des3$kMax-1)] > -8)) {
    df3b <- t(b)
    rownames(df3b) <- c("Information rate",
                        "Efficacy boundary (Z)",
                        "Futility boundary (Z)",
                        "Information")
  } else {
    df3b <- t(b[, c(1,2,4)])
    rownames(df3b) <- c("Information rate",
                        "Efficacy boundary (Z)",
                        "Information")
  }

  colnames(df3b) <- paste("Stage", seq_len(ncol(df3b)), sep=" ")

  print(df1a, ..., na.print = "" , quote = FALSE )
  print(df1b, ..., na.print = "" , quote = FALSE )
  print(df2a, ..., na.print = "" , quote = FALSE )
  print(df2b, ..., na.print = "" , quote = FALSE )
  print(df3a, ..., na.print = "" , quote = FALSE )
  print(df3b, ..., na.print = "" , quote = FALSE )
  invisible(x)
}


#' @title Print Power and Sample Size Results for Log-Rank Tests
#' @description Prints the summary statistics from power calculation.
#'
#' @param x The lrpower object to print.
#' @param ... Ensures that all arguments starting from "..." are named.
#'
#' @return A tabular printout of the summary statistics from power
#' calculation.
#'
#' @keywords internal
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @export
print.lrpower <- function(x, ...) {
  a <- x$overallResults
  s <- x$byStageResults
  k <- a$kMax

  if (k>1) {
    str1 <- paste0("Group-sequential design with ", k, " stages")
  } else {
    str1 <- "Fixed design"
  }

  if (a$rho1 != 0 || a$rho2 != 0) {
    str1 <- paste0(str1, " for weighted log-rank test, FH(",
                   a$rho1, ", ", a$rho2, ")")
  } else {
    str1 <- paste0(str1, " for log-rank test")
  }

  str2 <- paste0("Overall power: ",
                 round(a$overallReject, 4), ", ",
                 "overall significance level (1-sided): ",
                 round(a$alpha, 4))

  if (k>1) {
    str3 <- paste0("Maximum # events: ",
                   round(a$numberOfEvents, 1), ", ",
                   "expected # events: ",
                   round(a$expectedNumberOfEvents, 1))

    str4 <- paste0("Maximum # dropouts: ",
                   round(a$numberOfDropouts, 1), ", ",
                   "expected # dropouts: ",
                   round(a$expectedNumberOfDropouts, 1))

    str5 <- paste0("Maximum # subjects: ",
                   round(a$numberOfSubjects, 1), ", ",
                   "expected # subjects: ",
                   round(a$expectedNumberOfSubjects, 1))

    str6 <- paste0("Maximum information: ",
                   round(a$information, 2), ", ",
                   "expected information: ",
                   round(a$expectedInformation, 2))

    str7 <- paste0("Total study duration: ",
                   round(a$studyDuration, 1), ", ",
                   "expected study duration: ",
                   round(a$expectedStudyDuration, 1))

  } else {
    str3 <- paste0("Number of events: ",
                   round(a$numberOfEvents, 1))

    str4 <- paste0("Number of dropouts: ",
                   round(a$numberOfDropouts, 1))

    str5 <- paste0("Number of subjects: ",
                   round(a$numberOfSubjects, 1))

    str6 <- paste0("Information: ",
                   round(a$information, 2))

    str7 <- paste0("Study duration: ",
                   round(a$studyDuration, 1))
  }

  str8 <- paste0("Accrual duration: ",
                 round(a$accrualDuration, 1), ", ",
                 "follow-up duration: ",
                 round(a$followupTime, 1), ", ",
                 "fixed follow-up: ", a$fixedFollowup)

  str9 <- paste0("Allocation ratio: ",
                 round(x$settings$allocationRatioPlanned, 3))

  if (k > 1) {
    asf <- tolower(x$settings$typeAlphaSpending)
    asfpar <- round(x$settings$parameterAlphaSpending, 3)
    asfuser <- round(x$settings$userAlphaSpending, 4)

    bsf <- tolower(x$settings$typeBetaSpending)
    bsfpar <- round(x$settings$parameterBetaSpending, 3)

    if (asf == "of") {
      str10 <- paste0("Alpha spending: O'Brien-Fleming")
    } else if (asf == "p") {
      str10 <- paste0("Alpha spending: Pocock")
    } else if (asf == "wt") {
      str10 <- paste0("Alpha spending: Wang-Tsiatis(Delta = ", asfpar, ")")
    } else if (asf == "sfof") {
      str10 <- paste0("Alpha spending: Lan-DeMets O'Brien-Fleming")
    } else if (asf == "sfp") {
      str10 <- paste0("Alpha spending: Lan-DeMets Pocock")
    } else if (asf == "sfkd") {
      str10 <- paste0("Alpha spending: KD(rho = ", asfpar, ")")
    } else if (asf == "sfhsd") {
      str10 <- paste0("Alpha spending: HSD(gamma = ", asfpar, ")")
    } else if (asf == "user") {
      str10 <- paste0("Alpha spending: User defined(",
                      paste(asfuser, collapse = ","), ")")
    } else {
      str10 <- "Alpha spending: None"
    }

    if (bsf == "of") {
      str11 <- paste0("beta spending: O'Brien-Fleming")
    } else if (bsf == "p") {
      str11 <- paste0("beta spending: Pocock")
    } else if (bsf == "wt") {
      str11 <- paste0("beta spending: Wang-Tsiatis(Delta = ", bsfpar, ")")
    } else if (bsf == "sfof") {
      str11 <- paste0("beta spending: Lan-DeMets O'Brien-Fleming")
    } else if (bsf == "sfp") {
      str11 <- paste0("beta spending: Lan-DeMets Pocock")
    } else if (bsf == "sfkd") {
      str11 <- paste0("beta spending: KD(rho = ", bsfpar, ")")
    } else if (bsf == "sfhsd") {
      str11 <- paste0("beta spending: HSD(gamma = ", bsfpar, ")")
    } else if (bsf == "user") {
      bsfuser <- round(x$settings$userBetaSpending, 4)
      str11 <- paste0("beta spending: User defined(",
                      paste(bsfuser, collapse = ","), ")")
    } else {
      str11 <- "beta spending: None"
    }

    if (!any(is.na(x$settings$spendingTime)) &&
        !all(x$settings$spendingTime == s$informationRates)) {
      str12 <- paste0("Spending time: ",
                      paste(round(x$settings$spendingTime, 3), collapse = ","))
      df1 <- data.frame(x = rep("", 12))
      colnames(df1) <- NULL
      rownames(df1) <- c(str1, str2, str3, str4, str5, str6, str7,
                         str8, str9, paste(str10, str11, sep = ", "),
                         str12, "")
    } else {
      df1 <- data.frame(x = rep("", 11))
      colnames(df1) <- NULL
      rownames(df1) <- c(str1, str2, str3, str4, str5, str6, str7,
                         str8, str9, paste(str10, str11, sep = ", "), "")
    }
  } else {
    df1 <- data.frame(x = rep("", 10))
    colnames(df1) <- NULL
    rownames(df1) <- c(str1, str2, str3, str4, str5, str6, str7,
                       str8, str9, "")
  }

  if (k>1) {

    b <- s[, c("informationRates", "efficacyBounds", "futilityBounds",
               "cumulativeRejection", "cumulativeFutility",
               "cumulativeAlphaSpent",
               "numberOfEvents", "numberOfDropouts", "numberOfSubjects",
               "analysisTime", "efficacyHR", "futilityHR",
               "efficacyP", "futilityP", "information", "HR")]

    # format number of digits after decimal for each column
    j1 <- c(7,8,9,10)
    j2 <- 15
    j3 <- c(1,2,3,11,12,16)
    j4 <- c(4,5,6,13,14)

    b[j1] <- lapply(b[j1], formatC, format = "f", digits = 1)
    b[j2] <- lapply(b[j2], formatC, format = "f", digits = 2)
    b[j3] <- lapply(b[j3], formatC, format = "f", digits = 3)
    b[j4] <- lapply(b[j4], formatC, format = "f", digits = 4)

    if (x$settings$typeBetaSpending != 'none' ||
        (k > 1 && any(x$byStageResults$futilityBounds[1:(k-1)] > -8))) {
      df <- t(b)
      rownames(df) <- c("Information rate",
                        "Efficacy boundary (Z)",
                        "Futility boundary (Z)",
                        "Cumulative rejection",
                        "Cumulative futility",
                        "Cumulative alpha spent",
                        "Number of events",
                        "Number of dropouts",
                        "Number of subjects",
                        "Analysis time",
                        "Efficacy boundary (HR)",
                        "Futility boundary (HR)",
                        "Efficacy boundary (p)",
                        "Futility boundary (p)",
                        "Information",
                        "HR")

    } else {
      df <- t(b[, c(1,2,4,6,7,8,9,10,11,13,15,16)])
      rownames(df) <- c("Information rate",
                        "Efficacy boundary (Z)",
                        "Cumulative rejection",
                        "Cumulative alpha spent",
                        "Number of events",
                        "Number of dropouts",
                        "Number of subjects",
                        "Analysis time",
                        "Efficacy boundary (HR)",
                        "Efficacy boundary (p)",
                        "Information",
                        "HR")
    }



    colnames(df) <- paste("Stage", seq_len(ncol(df)), sep=" ")
  } else {

    b <- s[, c("efficacyBounds", "efficacyHR", "efficacyP",
               "HR")]

    # format number of digits after decimal for each column
    j3 <- c(1,2,4)
    j4 <- 3

    b[j3] <- lapply(b[j3], formatC, format = "f", digits = 3)
    b[j4] <- lapply(b[j4], formatC, format = "f", digits = 4)

    df <- t(b)

    rownames(df) <- c("Efficacy boundary (Z)",
                      "Efficacy boundary (HR)",
                      "Efficacy boundary (p)",
                      "HR")


    colnames(df) <- NA
  }

  print(df1, ..., na.print = "" , quote = FALSE )
  print(df, ..., na.print = "" , quote = FALSE )
  invisible(x)
}


#' @title Print Simulation Results for Log-Rank Tests
#' @description Prints the summary statistics from simulation.
#'
#' @param x The lrsim object to print.
#' @param ... Ensures that all arguments starting from "..." are named.
#'
#' @return A tabular printout of the summary statistics from simulation runs.
#'
#' @keywords internal
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @export
print.lrsim <- function(x, ...) {
  a <- x$overview
  k <- a$kMax

  if (k>1) {
    str1 <- paste0("Group-sequential design with ", k, " stages")
  } else {
    str1 <- "Fixed design"
  }

  if (a$rho1 != 0 || a$rho2 != 0) {
    str1 <- paste0(str1, " for weighted log-rank test, FH(",
                   a$rho1, ", ", a$rho2, ")")
  } else {
    str1 <- paste0(str1, " for log-rank test")
  }

  str2 <- paste0("Empirical power: ", round(a$overallReject, 4))

  str3 <- paste0("Expected # events: ",
                 round(a$expectedNumberOfEvents, 1))

  str4 <- paste0("Expected # dropouts: ",
                 round(a$expectedNumberOfDropouts, 1))

  str5 <- paste0("Expected # subjects: ",
                 round(a$expectedNumberOfSubjects, 1))

  str6 <- paste0("Expected study duration: ",
                 round(a$expectedStudyDuration, 1))

  str7 <- paste0("n: ", a$n, ", ",
                 "fixed follow-up: ", a$fixedFollowup)

  str8 <- paste0("Number of simulations: ", a$numberOfIterations)

  df1 <- data.frame(x = rep("", 9))
  colnames(df1) <- NULL
  rownames(df1) <- c(str1, str2, str3, str4, str5, str6, str7, str8, "")


  if (k>1) {
    b <- data.frame(a$cumulativeRejection,
                    a$cumulativeFutility,
                    a$numberOfEvents,
                    a$numberOfDropouts,
                    a$numberOfSubjects,
                    a$analysisTime)

    # format number of digits after decimal for each column
    j1 <- c(3,4,5,6)
    j4 <- c(1,2)

    b[j1] <- lapply(b[j1], formatC, format = "f", digits = 1)
    b[j4] <- lapply(b[j4], formatC, format = "f", digits = 4)

    df <- t(b)
    rownames(df) <- c("Cumulative rejection",
                      "Cumulative futility",
                      "Number of events",
                      "Number of dropouts",
                      "Number of subjects",
                      "Analysis time")
    colnames(df) <- paste("Stage", seq_len(ncol(df)), sep=" ")
  }

  print(df1, ..., na.print = "" , quote = FALSE )

  if (k>1) {
    print(df, ..., na.print = "" , quote = FALSE )
  }

  invisible(x)
}


#' @title Print Power and Sample Size Results for Negative Binomial Rate
#' Ratio
#' @description Prints the summary statistics from power calculation of
#' negative binomial rate ratio.
#'
#' @param x The nbpower object to print.
#' @param ... Ensures that all arguments starting from "..." are named.
#'
#' @return A tabular printout of the summary statistics from power
#' calculation.
#'
#' @keywords internal
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @export
print.nbpower <- function(x, ...) {
  a <- x$overallResults
  s <- x$byStageResults
  k <- a$kMax

  if (k>1) {
    str1 <- paste0("Group-sequential design with ", k, " stages")
  } else {
    str1 <- "Fixed design"
  }

  str1 <- paste(str1, "for negative binomial rate ratio")


  str2 <- paste0("Rate ratio under H0: ",
                 round(a$rateRatioH0, 3), ", ",
                 "rate ratio under H1: ",
                 round(a$rateRatio, 3))


  if (length(x$settings$stratumFraction) > 1) {
    str3a <- paste0("Stratum fraction: ",
                    paste(round(x$settings$stratumFraction, 3),
                          collapse = " "))

  }

  str3 <- paste0("Event rate for treatment: ",
                 paste(round(x$settings$lambda1, 4),
                       collapse = " "), ", ",
                 "event rate for control: ",
                 paste(round(x$settings$lambda2, 4),
                       collapse = " "))

  str4 <- paste0("Dispersion for treatment: ",
                 paste(round(x$settings$kappa1, 3),
                       collapse = " "), ", ",
                 "dispersion for control: ",
                 paste(round(x$settings$kappa2, 3),
                       collase = " "))

  str5 <- paste0("Overall power: ",
                 round(a$overallReject, 4), ", ",
                 "overall significance level (1-sided): ",
                 round(a$alpha, 4))

  if (k>1) {
    str6 <- paste0("Maximum # events: ",
                   round(a$numberOfEvents, 1), ", ",
                   "expected # events: ",
                   round(a$expectedNumberOfEvents, 1))

    str7 <- paste0("Maximum # dropouts: ",
                   round(a$numberOfDropouts, 1), ", ",
                   "expected # dropouts: ",
                   round(a$expectedNumberOfDropouts, 1))

    str8 <- paste0("Maximum # subjects: ",
                   round(a$numberOfSubjects, 1), ", ",
                   "expected # subjects: ",
                   round(a$expectedNumberOfSubjects, 1))

    str9 <- paste0("Maximum exposure: ",
                   round(a$exposure, 1), ", ",
                   "expected exposure: ",
                   round(a$expectedExposure, 1))

    str10 <- paste0("Maximum information: ",
                    round(a$information, 2), ", ",
                    "expected information: ",
                    round(a$expectedInformation, 2))

    str11 <- paste0("Total study duration: ",
                    round(a$studyDuration, 1), ", ",
                    "expected study duration: ",
                    round(a$expectedStudyDuration, 1))
  } else {
    str6 <- paste0("Number of events: ",
                   round(a$numberOfEvents, 1))

    str7 <- paste0("Number of dropouts: ",
                   round(a$numberOfDropouts, 1))

    str8 <- paste0("Number of subjects: ",
                   round(a$numberOfSubjects, 1))

    str9 <- paste0("Exposure: ",
                   round(a$exposure, 1))

    str10 <- paste0("Information: ",
                    round(a$information, 2))

    str11 <- paste0("Study duration: ",
                    round(a$studyDuration, 1))
  }

  str12 <- paste0("Accrual duration: ",
                  round(a$accrualDuration, 1), ", ",
                  "follow-up duration: ",
                  round(a$followupTime, 1), ", ",
                  "fixed follow-up: ", a$fixedFollowup)

  str13 <- paste0("Allocation ratio: ",
                  round(x$settings$allocationRatioPlanned, 3), ", ",
                  "variance of standardized test statistic: ",
                  ifelse(x$settings$nullVariance, "under H0", "under H1"))


  if (k > 1) {
    asf <- tolower(x$settings$typeAlphaSpending)
    asfpar <- round(x$settings$parameterAlphaSpending, 3)
    asfuser <- round(x$settings$userAlphaSpending, 4)

    bsf <- tolower(x$settings$typeBetaSpending)
    bsfpar <- round(x$settings$parameterBetaSpending, 3)

    if (asf == "of") {
      str14 <- paste0("Alpha spending: O'Brien-Fleming")
    } else if (asf == "p") {
      str14 <- paste0("Alpha spending: Pocock")
    } else if (asf == "wt") {
      str14 <- paste0("Alpha spending: Wang-Tsiatis(Delta = ", asfpar, ")")
    } else if (asf == "sfof") {
      str14 <- paste0("Alpha spending: Lan-DeMets O'Brien-Fleming")
    } else if (asf == "sfp") {
      str14 <- paste0("Alpha spending: Lan-DeMets Pocock")
    } else if (asf == "sfkd") {
      str14 <- paste0("Alpha spending: KD(rho = ", asfpar, ")")
    } else if (asf == "sfhsd") {
      str14 <- paste0("Alpha spending: HSD(gamma = ", asfpar, ")")
    } else if (asf == "user") {
      str14 <- paste0("Alpha spending: User defined(",
                      paste(asfuser, collapse = ","), ")")
    } else {
      str14 <- "Alpha spending: None"
    }

    if (bsf == "of") {
      str15 <- paste0("beta spending: O'Brien-Fleming")
    } else if (bsf == "p") {
      str15 <- paste0("beta spending: Pocock")
    } else if (bsf == "wt") {
      str15 <- paste0("beta spending: Wang-Tsiatis(Delta = ", bsfpar, ")")
    } else if (bsf == "sfof") {
      str15 <- paste0("beta spending: Lan-DeMets O'Brien-Fleming")
    } else if (bsf == "sfp") {
      str15 <- paste0("beta spending: Lan-DeMets Pocock")
    } else if (bsf == "sfkd") {
      str15 <- paste0("beta spending: KD(rho = ", bsfpar, ")")
    } else if (bsf == "sfhsd") {
      str15 <- paste0("beta spending: HSD(gamma = ", bsfpar, ")")
    } else if (bsf == "user") {
      str15 <- paste0("beta spending: User defined")
    } else {
      str15 <- "beta spending: None"
    }

    if (!any(is.na(x$settings$spendingTime)) &&
        !all(x$settings$spendingTime == s$informationRates)) {
      str16 <- paste0("Spending time: ",
                      paste(round(x$settings$spendingTime, 3), collapse = ","))

      if (length(x$settings$stratumFraction) > 1) {
        df1 <- data.frame(x = rep("", 17))
        colnames(df1) <- NULL
        rownames(df1) <- c(str1, str2, str3a, str3, str4, str5, str6, str7,
                           str8, str9, str10, str11, str12, str13,
                           paste(str14, str15, sep = ", "), str16, "")

      } else {
        df1 <- data.frame(x = rep("", 16))
        colnames(df1) <- NULL
        rownames(df1) <- c(str1, str2, str3, str4, str5, str6, str7,
                           str8, str9, str10, str11, str12, str13,
                           paste(str14, str15, sep = ", "), str16, "")
      }
    } else {
      if (length(x$settings$stratumFraction) > 1) {
        df1 <- data.frame(x = rep("", 16))
        colnames(df1) <- NULL
        rownames(df1) <- c(str1, str2, str3a, str3, str4, str5, str6, str7,
                           str8, str9, str10, str11, str12, str13,
                           paste(str14, str15, sep = ", "), "")
      } else {
        df1 <- data.frame(x = rep("", 15))
        colnames(df1) <- NULL
        rownames(df1) <- c(str1, str2, str3, str4, str5, str6, str7,
                           str8, str9, str10, str11, str12, str13,
                           paste(str14, str15, sep = ", "), "")
      }
    }
  } else {
    if (length(x$settings$stratumFraction) > 1) {
      df1 <- data.frame(x = rep("", 15))
      colnames(df1) <- NULL
      rownames(df1) <- c(str1, str2, str3a, str3, str4, str5, str6, str7,
                         str8, str9, str10, str11, str12, str13, "")
    } else {
      df1 <- data.frame(x = rep("", 14))
      colnames(df1) <- NULL
      rownames(df1) <- c(str1, str2, str3, str4, str5, str6, str7,
                         str8, str9, str10, str11, str12, str13, "")
    }
  }

  if (k>1) {
    b <- s[, c("informationRates", "efficacyBounds", "futilityBounds",
               "cumulativeRejection", "cumulativeFutility",
               "cumulativeAlphaSpent", "numberOfEvents",
               "numberOfDropouts", "numberOfSubjects", "exposure",
               "analysisTime", "efficacyRateRatio", "futilityRateRatio",
               "efficacyP", "futilityP", "information")]

    # format number of digits after decimal for each column
    j1 <- c(7,8,9,10,11)
    j2 <- 16
    j3 <- c(1,2,3,12,13)
    j4 <- c(4,5,6,14,15)

    b[j1] <- lapply(b[j1], formatC, format = "f", digits = 1)
    b[j2] <- lapply(b[j2], formatC, format = "f", digits = 2)
    b[j3] <- lapply(b[j3], formatC, format = "f", digits = 3)
    b[j4] <- lapply(b[j4], formatC, format = "f", digits = 4)

    if (x$settings$typeBetaSpending != 'none' ||
        (k > 1 && any(x$byStageResults$futilityBounds[1:(k-1)] > -8))) {
      df <- t(b)
      rownames(df) <- c("Information rate",
                        "Efficacy boundary (Z)",
                        "Futility boundary (Z)",
                        "Cumulative rejection",
                        "Cumulative futility",
                        "Cumulative alpha spent",
                        "Number of events",
                        "Number of dropouts",
                        "Number of subjects",
                        "Exposure",
                        "Analysis time",
                        "Efficacy boundary (rate ratio)",
                        "Futility boundary (rate ratio)",
                        "Efficacy boundary (p)",
                        "Futility boundary (p)",
                        "Information")

    } else {
      df <- t(b[, c(1,2,4,6,7,8,9,10,11,12,14,16)])
      rownames(df) <- c("Information rate",
                        "Efficacy boundary (Z)",
                        "Cumulative rejection",
                        "Cumulative alpha spent",
                        "Number of events",
                        "Number of dropouts",
                        "Number of subjects",
                        "Exposure",
                        "Analysis time",
                        "Efficacy boundary (rate ratio)",
                        "Efficacy boundary (p)",
                        "Information")
    }

    colnames(df) <- paste("Stage", seq_len(ncol(df)), sep=" ")
  } else {
    b <- s[, c("efficacyBounds", "efficacyRateRatio", "efficacyP")]

    # format number of digits after decimal for each column
    j3 <- c(1,2)
    j4 <- 3

    b[j3] <- lapply(b[j3], formatC, format = "f", digits = 3)
    b[j4] <- lapply(b[j4], formatC, format = "f", digits = 4)

    df <- t(b)

    rownames(df) <- c("Efficacy boundary (Z)",
                      "Efficacy boundary (rate ratio)",
                      "Efficacy boundary (p)")
    colnames(df) <- NA
  }

  print(df1, ..., na.print = "" , quote = FALSE )
  print(df, ..., na.print = "" , quote = FALSE )
  invisible(x)
}


#' @title Print Power and Sample Size Results for Equivalence in Negative
#' Binomial Rate Ratio
#' @description Prints the summary statistics from power calculation of
#' equivalence in negative binomial rate ratio.
#'
#' @param x The nbpowerequiv object to print.
#' @param ... Ensures that all arguments starting from "..." are named.
#'
#' @return A tabular printout of the summary statistics from power
#' calculation.
#'
#' @keywords internal
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @export
print.nbpowerequiv <- function(x, ...) {
  a <- x$overallResults
  s <- x$byStageResults
  k <- a$kMax

  if (k>1) {
    str1 <- paste0("Group-sequential design with ", k, " stages")
  } else {
    str1 <- "Fixed design"
  }

  str1 <- paste(str1, "for equivalence in negative binomial rate ratio")

  str2 <- paste0("Lower limit for rate ratio: ",
                 round(a$rateRatioLower, 3), ", ",
                 "upper limit for rate ratio: ",
                 round(a$rateRatioUpper, 3), ", ",
                 "rate ratio: ",
                 round(a$rateRatio, 3))


  if (length(x$settings$stratumFraction) > 1) {
    str3a <- paste0("Stratum fraction: ",
                    paste(round(x$settings$stratumFraction, 3),
                          collapse = " "))

  }

  str3 <- paste0("Event rate for treatment: ",
                 paste(round(x$settings$lambda1, 4),
                       collapse = " "), ", ",
                 "event rate for control: ",
                 paste(round(x$settings$lambda2, 4),
                       collapse = " "))

  str4 <- paste0("Dispersion for treatment: ",
                 paste(round(x$settings$kappa1, 3),
                       collapse = " "), ", ",
                 "dispersion for control: ",
                 paste(round(x$settings$kappa2, 3),
                       collapse = " "))

  str5 <- paste0("Overall power: ",
                 round(a$overallReject, 4), ", ",
                 "overall alpha: ",
                 round(a$alpha, 4))

  if (k>1) {
    str6 <- paste0("Maximum # events: ",
                   round(a$numberOfEvents, 1), ", ",
                   "expected # events: ",
                   round(a$expectedNumberOfEvents, 1))

    str7 <- paste0("Maximum # subjects: ",
                   round(a$numberOfSubjects, 1), ", ",
                   "expected # subjects: ",
                   round(a$expectedNumberOfSubjects, 1))

    str8 <- paste0("Maximum exposure: ",
                   round(a$exposure, 1), ", ",
                   "expected exposure: ",
                   round(a$expectedExposure, 1))

    str9 <- paste0("Maximum information: ",
                   round(a$information, 2), ", ",
                   "expected information: ",
                   round(a$expectedInformation, 2))

    str10 <- paste0("Total study duration: ",
                    round(a$studyDuration, 1), ", ",
                    "expected study duration: ",
                    round(a$expectedStudyDuration, 1))
  } else {
    str6 <- paste0("Number of events: ",
                   round(a$numberOfEvents, 1))

    str7 <- paste0("Number of subjects: ",
                   round(a$numberOfSubjects, 1))

    str8 <- paste0("Exposure: ",
                   round(a$exposure, 1))

    str9 <- paste0("Information: ",
                   round(a$information, 2))

    str10 <- paste0("Study duration: ",
                    round(a$studyDuration, 1))
  }

  str11 <- paste0("Accrual duration: ",
                  round(a$accrualDuration, 1), ", ",
                  "follow-up duration: ",
                  round(a$followupTime, 1), ", ",
                  "fixed follow-up: ", a$fixedFollowup)

  str12 <- paste0("Allocation ratio: ",
                  round(x$settings$allocationRatioPlanned, 3))

  if (k > 1) {
    asf <- tolower(x$settings$typeAlphaSpending)
    asfpar <- round(x$settings$parameterAlphaSpending, 3)
    asfuser <- round(x$settings$userAlphaSpending, 4)

    if (asf == "of") {
      str13 <- paste0("Alpha spending: O'Brien-Fleming")
    } else if (asf == "p") {
      str13 <- paste0("Alpha spending: Pocock")
    } else if (asf == "wt") {
      str13 <- paste0("Alpha spending: Wang-Tsiatis(Delta = ", asfpar, ")")
    } else if (asf == "sfof") {
      str13 <- paste0("Alpha spending: Lan-DeMets O'Brien-Fleming")
    } else if (asf == "sfp") {
      str13 <- paste0("Alpha spending: Lan-DeMets Pocock")
    } else if (asf == "sfkd") {
      str13 <- paste0("Alpha spending: KD(rho = ", asfpar, ")")
    } else if (asf == "sfhsd") {
      str13 <- paste0("Alpha spending: HSD(gamma = ", asfpar, ")")
    } else if (asf == "user") {
      str13 <- paste0("Alpha spending: User defined(",
                      paste(asfuser, collapse = ","), ")")
    } else {
      str13 <- "Alpha spending: None"
    }

    if (!any(is.na(x$settings$spendingTime)) &&
        !all(x$settings$spendingTime == s$informationRates)) {
      str14 <- paste0("Spending time: ",
                      paste(round(x$settings$spendingTime, 3), collapse = ","))
      if (length(x$settings$stratumFraction) > 1) {
        df1 <- data.frame(x = rep("", 16))
        colnames(df1) <- NULL
        rownames(df1) <- c(str1, str2, str3a, str3, str4, str5, str6, str7,
                           str8, str9, str10, str11, str12,
                           str13, str14, "")
      } else {
        df1 <- data.frame(x = rep("", 15))
        colnames(df1) <- NULL
        rownames(df1) <- c(str1, str2, str3, str4, str5, str6, str7,
                           str8, str9, str10, str11, str12,
                           str13, str14, "")
      }
    } else {
      if (length(x$settings$stratumFraction) > 1) {
        df1 <- data.frame(x = rep("", 15))
        colnames(df1) <- NULL
        rownames(df1) <- c(str1, str2, str3a, str3, str4, str5, str6, str7,
                           str8, str9, str10, str11, str12,
                           str13, "")
      } else {
        df1 <- data.frame(x = rep("", 14))
        colnames(df1) <- NULL
        rownames(df1) <- c(str1, str2, str3, str4, str5, str6, str7,
                           str8, str9, str10, str11, str12,
                           str13, "")
      }
    }
  } else {
    if (length(x$settings$stratumFraction) > 1) {
      df1 <- data.frame(x = rep("", 14))
      colnames(df1) <- NULL
      rownames(df1) <- c(str1, str2, str3a, str3, str4, str5, str6, str7,
                         str8, str9, str10, str11, str12, "")
    } else {
      df1 <- data.frame(x = rep("", 13))
      colnames(df1) <- NULL
      rownames(df1) <- c(str1, str2, str3, str4, str5, str6, str7,
                         str8, str9, str10, str11, str12, "")
    }
  }

  if (k>1) {
    b <- s[, c("informationRates", "efficacyBounds",
               "cumulativeRejection", "cumulativeAlphaSpent",
               "cumulativeAttainedAlphaH10",
               "cumulativeAttainedAlphaH20", "numberOfEvents",
               "numberOfDropouts", "numberOfSubjects",
               "exposure", "analysisTime",
               "efficacyRateRatioLower", "efficacyRateRatioUpper",
               "efficacyP", "information")]

    # format number of digits after decimal for each column
    j1 <- c(7,8,9,10,11)
    j2 <- 15
    j3 <- c(1,2,12,13)
    j4 <- c(3,4,5,6,14)

    b[j1] <- lapply(b[j1], formatC, format = "f", digits = 1)
    b[j2] <- lapply(b[j2], formatC, format = "f", digits = 2)
    b[j3] <- lapply(b[j3], formatC, format = "f", digits = 3)
    b[j4] <- lapply(b[j4], formatC, format = "f", digits = 4)

    df <- t(b)
    rownames(df) <- c("Information rate",
                      "Boundary for each 1-sided test (Z)",
                      "Cumulative rejection",
                      "Cumulative alpha for each 1-sided test",
                      "Cumulative alpha attained under H10",
                      "Cumulative alpha attained under H20",
                      "Number of events",
                      "Number of dropouts",
                      "Number of subjects",
                      "Exposure",
                      "Analysis time",
                      "Boundary for lower limit (rate ratio)",
                      "Boundary for upper limit (rate ratio)",
                      "Boundary for each 1-sided test (p)",
                      "Information")

    colnames(df) <- paste("Stage", seq_len(ncol(df)), sep=" ")
  } else {
    b <- s[, c("efficacyBounds", "efficacyRateRatioLower",
               "efficacyRateRatioUpper",  "efficacyP")]

    # format number of digits after decimal for each column
    j3 <- c(1,2,3)
    j4 <- 4

    b[j3] <- lapply(b[j3], formatC, format = "f", digits = 3)
    b[j4] <- lapply(b[j4], formatC, format = "f", digits = 4)

    df <- t(b)

    rownames(df) <- c("Boundary for each 1-sided test (Z)",
                      "Boundary for lower limit (rate ratio)",
                      "Boundary for upper limit (rate ratio)",
                      "Boundary for each 1-sided test (p)")
    colnames(df) <- NA
  }

  print(df1, ..., na.print = "" , quote = FALSE )
  print(df, ..., na.print = "" , quote = FALSE )
  invisible(x)
}


#' @title Print Power and Sample Size Results for One-Sample Negative
#' Binomial Rate
#' @description Prints the summary statistics from power calculation of
#' one-sample negative binomial rate.
#'
#' @param x The nbpower1s object to print.
#' @param ... Ensures that all arguments starting from "..." are named.
#'
#' @return A tabular printout of the summary statistics from power
#' calculation.
#'
#' @keywords internal
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @export
print.nbpower1s <- function(x, ...) {
  a <- x$overallResults
  s <- x$byStageResults
  k <- a$kMax

  if (k>1) {
    str1 <- paste0("Group-sequential design with ", k, " stages")
  } else {
    str1 <- "Fixed design"
  }

  str1 <- paste(str1, "for one-sample negative binomial rate")

  str2 <- paste0("Rate under H0: ",
                 round(a$lambdaH0, 4), ", ",
                 "rate under H1: ",
                 round(a$lambda, 4))

  if (length(x$settings$stratumFraction) > 1) {
    str3 <- paste0("Stratum fraction: ",
                   paste(round(x$settings$stratumFraction, 3),
                         collapse = " "), ", ",
                   "event rate: ",
                   paste(round(x$settings$lambda, 4),
                         collapse = " "), ", ",
                   "dispersion: ",
                   paste(round(x$settings$kappa, 3),
                         collapse = " "))
  } else {
    str3 <- paste0("Dispersion: ",
                   round(x$settings$kappa, 3))
  }

  str4 <- paste0("Overall power: ",
                 round(a$overallReject, 4), ", ",
                 "overall significance level (1-sided): ",
                 round(a$alpha, 4))

  if (k>1) {
    str5 <- paste0("Maximum # events: ",
                   round(a$numberOfEvents, 1), ", ",
                   "expected # events: ",
                   round(a$expectedNumberOfEvents, 1))

    str6 <- paste0("Maximum # dropouts: ",
                   round(a$numberOfDropouts, 1), ", ",
                   "expected # dropouts: ",
                   round(a$expectedNumberOfDropouts, 1))

    str7 <- paste0("Maximum # subjects: ",
                   round(a$numberOfSubjects, 1), ", ",
                   "expected # subjects: ",
                   round(a$expectedNumberOfSubjects, 1))

    str8 <- paste0("Maximum exposure: ",
                   round(a$exposure, 1), ", ",
                   "expected exposure: ",
                   round(a$expectedExposure, 1))

    str9 <- paste0("Maximum information: ",
                   round(a$information, 2), ", ",
                   "expected information: ",
                   round(a$expectedInformation, 2))

    str10 <- paste0("Total study duration: ",
                    round(a$studyDuration, 1), ", ",
                    "expected study duration: ",
                    round(a$expectedStudyDuration, 1))
  } else {
    str5 <- paste0("Number of events: ",
                   round(a$numberOfEvents, 1))

    str6 <- paste0("Number of dropouts: ",
                   round(a$numberOfDropouts, 1))

    str7 <- paste0("Number of subjects: ",
                   round(a$numberOfSubjects, 1))

    str8 <- paste0("Exposure: ",
                   round(a$exposure, 1))

    str9 <- paste0("Information: ",
                   round(a$information, 2))

    str10 <- paste0("Study duration: ",
                    round(a$studyDuration, 1))
  }

  str11 <- paste0("Accrual duration: ",
                  round(a$accrualDuration, 1), ", ",
                  "follow-up duration: ",
                  round(a$followupTime, 1), ", ",
                  "fixed follow-up: ", a$fixedFollowup)


  if (k > 1) {
    asf <- tolower(x$settings$typeAlphaSpending)
    asfpar <- round(x$settings$parameterAlphaSpending, 3)
    asfuser <- round(x$settings$userAlphaSpending, 4)

    bsf <- tolower(x$settings$typeBetaSpending)
    bsfpar <- round(x$settings$parameterBetaSpending, 3)

    if (asf == "of") {
      str12 <- paste0("Alpha spending: O'Brien-Fleming")
    } else if (asf == "p") {
      str12 <- paste0("Alpha spending: Pocock")
    } else if (asf == "wt") {
      str12 <- paste0("Alpha spending: Wang-Tsiatis(Delta = ", asfpar, ")")
    } else if (asf == "sfof") {
      str12 <- paste0("Alpha spending: Lan-DeMets O'Brien-Fleming")
    } else if (asf == "sfp") {
      str12 <- paste0("Alpha spending: Lan-DeMets Pocock")
    } else if (asf == "sfkd") {
      str12 <- paste0("Alpha spending: KD(rho = ", asfpar, ")")
    } else if (asf == "sfhsd") {
      str12 <- paste0("Alpha spending: HSD(gamma = ", asfpar, ")")
    } else if (asf == "user") {
      str12 <- paste0("Alpha spending: User defined(",
                      paste(asfuser, collapse = ","), ")")
    } else {
      str12 <- "Alpha spending: None"
    }

    if (bsf == "of") {
      str13 <- paste0("beta spending: O'Brien-Fleming")
    } else if (bsf == "p") {
      str13 <- paste0("beta spending: Pocock")
    } else if (bsf == "wt") {
      str13 <- paste0("beta spending: Wang-Tsiatis(Delta = ", bsfpar, ")")
    } else if (bsf == "sfof") {
      str13 <- paste0("beta spending: Lan-DeMets O'Brien-Fleming")
    } else if (bsf == "sfp") {
      str13 <- paste0("beta spending: Lan-DeMets Pocock")
    } else if (bsf == "sfkd") {
      str13 <- paste0("beta spending: KD(rho = ", bsfpar, ")")
    } else if (bsf == "sfhsd") {
      str13 <- paste0("beta spending: HSD(gamma = ", bsfpar, ")")
    } else if (bsf == "user") {
      str13 <- paste0("beta spending: User defined")
    } else {
      str13 <- "beta spending: None"
    }

    if (!any(is.na(x$settings$spendingTime)) &&
        !all(x$settings$spendingTime == s$informationRates)) {
      str14 <- paste0("Spending time: ",
                      paste(round(x$settings$spendingTime, 3), collapse = ","))
      df1 <- data.frame(x = rep("", 14))
      colnames(df1) <- NULL
      rownames(df1) <- c(str1, str2, str3, str4, str5, str6, str7,
                         str8, str9, str10, str11,
                         paste(str12, str13, sep = ", "), str14, "")
    } else {
      df1 <- data.frame(x = rep("", 13))
      colnames(df1) <- NULL
      rownames(df1) <- c(str1, str2, str3, str4, str5, str6, str7,
                         str8, str9, str10, str11,
                         paste(str12, str13, sep = ", "), "")
    }
  } else {
    df1 <- data.frame(x = rep("", 12))
    colnames(df1) <- NULL
    rownames(df1) <- c(str1, str2, str3, str4, str5, str6, str7,
                       str8, str9, str10, str11, "")
  }

  if (k>1) {
    b <- s[, c("informationRates", "efficacyBounds", "futilityBounds",
               "cumulativeRejection", "cumulativeFutility",
               "cumulativeAlphaSpent",
               "numberOfEvents", "numberOfDropouts", "numberOfSubjects",
               "exposure", "analysisTime", "efficacyRate", "futilityRate",
               "efficacyP", "futilityP", "information")]

    # format number of digits after decimal for each column
    j1 <- c(7,8,9,10,11)
    j2 <- 16
    j3 <- c(1,2,3)
    j4 <- c(4,5,6,12,13,14,15)

    b[j1] <- lapply(b[j1], formatC, format = "f", digits = 1)
    b[j2] <- lapply(b[j2], formatC, format = "f", digits = 2)
    b[j3] <- lapply(b[j3], formatC, format = "f", digits = 3)
    b[j4] <- lapply(b[j4], formatC, format = "f", digits = 4)

    if (x$settings$typeBetaSpending != 'none' ||
        (k > 1 && any(x$byStageResults$futilityBounds[1:(k-1)] > -8))) {
      df <- t(b)
      rownames(df) <- c("Information rate",
                        "Efficacy boundary (Z)",
                        "Futility boundary (Z)",
                        "Cumulative rejection",
                        "Cumulative futility",
                        "Cumulative alpha spent",
                        "Number of events",
                        "Number of dropouts",
                        "Number of subjects",
                        "Exposure",
                        "Analysis time",
                        "Efficacy boundary (rate)",
                        "Futility boundary (rate)",
                        "Efficacy boundary (p)",
                        "Futility boundary (p)",
                        "Information")

    } else {
      df <- t(b[, c(1,2,4,6,7,8,9,10,11,12,14,16)])
      rownames(df) <- c("Information rate",
                        "Efficacy boundary (Z)",
                        "Cumulative rejection",
                        "Cumulative alpha spent",
                        "Number of events",
                        "Number of dropouts",
                        "Number of subjects",
                        "Exposure",
                        "Analysis time",
                        "Efficacy boundary (rate)",
                        "Efficacy boundary (p)",
                        "Information")
    }

    colnames(df) <- paste("Stage", seq_len(ncol(df)), sep=" ")
  } else {
    b <- s[, c("efficacyBounds", "efficacyRate", "efficacyP")]

    # format number of digits after decimal for each column
    j3 <- c(1,2)
    j4 <- 3

    b[j3] <- lapply(b[j3], formatC, format = "f", digits = 3)
    b[j4] <- lapply(b[j4], formatC, format = "f", digits = 4)

    df <- t(b)

    rownames(df) <- c("Efficacy boundary (Z)",
                      "Efficacy boundary (rate)",
                      "Efficacy boundary (p)")
    colnames(df) <- NA
  }

  print(df1, ..., na.print = "" , quote = FALSE )
  print(df, ..., na.print = "" , quote = FALSE )
  invisible(x)
}


#' @title Print Power and Sample Size Results for Milestone Survival
#' Difference
#' @description Prints the summary statistics from power calculation.
#'
#' @param x The kmpower object to print.
#' @param ... Ensures that all arguments starting from "..." are named.
#'
#' @return A tabular printout of the summary statistics from power
#' calculation.
#'
#' @keywords internal
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @export
print.kmpower <- function(x, ...) {
  a <- x$overallResults
  s <- x$byStageResults
  k <- a$kMax

  if (k>1) {
    str1 <- paste0("Group-sequential design with ", k, " stages")
  } else {
    str1 <- "Fixed design"
  }

  str1 <- paste0(str1, " for difference in milestone survival")

  str2 <- paste0("Milestone: ", round(a$milestone, 3), ", ",
                 "survival difference under H0: ",
                 round(a$survDiffH0, 3))

  str3 <- paste0("Milestone survival on treatment: ",
                 round(a$surv1, 3), ", ",
                 "on control: ", round(a$surv2, 3))

  str4 <- paste0("Overall power: ",
                 round(a$overallReject, 4), ", ",
                 "overall significance level (1-sided): ",
                 round(a$alpha, 4))

  if (k>1) {
    str5 <- paste0("Maximum # events: ",
                   round(a$numberOfEvents, 1), ", ",
                   "expected # events: ",
                   round(a$expectedNumberOfEvents, 1))

    str6 <- paste0("Maximum # subjects: ",
                   round(a$numberOfSubjects, 1), ", ",
                   "expected # subjects: ",
                   round(a$expectedNumberOfSubjects, 1))

    str7 <- paste0("Maximum information: ",
                   round(a$information, 2), ", ",
                   "expected information: ",
                   round(a$expectedInformation, 2))

    str8 <- paste0("Total study duration: ",
                   round(a$studyDuration, 1), ", ",
                   "expected study duration: ",
                   round(a$expectedStudyDuration, 1))

  } else {
    str5 <- paste0("Number of events: ",
                   round(a$numberOfEvents, 1))

    str6 <- paste0("Number of subjects: ",
                   round(a$numberOfSubjects, 1))

    str7 <- paste0("Information: ",
                   round(a$information, 2))

    str8 <- paste0("Study duration: ",
                   round(a$studyDuration, 1))
  }

  str9 <- paste0("Accrual duration: ",
                 round(a$accrualDuration, 1), ", ",
                 "follow-up duration: ",
                 round(a$followupTime, 1), ", ",
                 "fixed follow-up: ", a$fixedFollowup)

  str10 <- paste0("Allocation ratio: ",
                  round(x$settings$allocationRatioPlanned, 3))

  if (k > 1) {
    asf <- tolower(x$settings$typeAlphaSpending)
    asfpar <- round(x$settings$parameterAlphaSpending, 3)
    asfuser <- round(x$settings$userAlphaSpending, 4)

    bsf <- tolower(x$settings$typeBetaSpending)
    bsfpar <- round(x$settings$parameterBetaSpending, 3)

    if (asf == "of") {
      str11 <- paste0("Alpha spending: O'Brien-Fleming")
    } else if (asf == "p") {
      str11 <- paste0("Alpha spending: Pocock")
    } else if (asf == "wt") {
      str11 <- paste0("Alpha spending: Wang-Tsiatis(Delta = ", asfpar, ")")
    } else if (asf == "sfof") {
      str11 <- paste0("Alpha spending: Lan-DeMets O'Brien-Fleming")
    } else if (asf == "sfp") {
      str11 <- paste0("Alpha spending: Lan-DeMets Pocock")
    } else if (asf == "sfkd") {
      str11 <- paste0("Alpha spending: KD(rho = ", asfpar, ")")
    } else if (asf == "sfhsd") {
      str11 <- paste0("Alpha spending: HSD(gamma = ", asfpar, ")")
    } else if (asf == "user") {
      str11 <- paste0("Alpha spending: User defined(",
                      paste(asfuser, collapse = ","), ")")
    } else {
      str11 <- "Alpha spending: None"
    }

    if (bsf == "of") {
      str12 <- paste0("beta spending: O'Brien-Fleming")
    } else if (bsf == "p") {
      str12 <- paste0("beta spending: Pocock")
    } else if (bsf == "wt") {
      str12 <- paste0("beta spending: Wang-Tsiatis(Delta = ", bsfpar, ")")
    } else if (bsf == "sfof") {
      str12 <- paste0("beta spending: Lan-DeMets O'Brien-Fleming")
    } else if (bsf == "sfp") {
      str12 <- paste0("beta spending: Lan-DeMets Pocock")
    } else if (bsf == "sfkd") {
      str12 <- paste0("beta spending: KD(rho = ", bsfpar, ")")
    } else if (bsf == "sfhsd") {
      str12 <- paste0("beta spending: HSD(gamma = ", bsfpar, ")")
    } else if (bsf == "user") {
      bsfuser <- round(x$settings$userBetaSpending, 4)
      str12 <- paste0("beta spending: User defined(",
                      paste(bsfuser, collapse = ","), ")")
    } else {
      str12 <- "beta spending: None"
    }

    if (!any(is.na(x$settings$spendingTime)) &&
        !all(x$settings$spendingTime == s$informationRates)) {
      str13 <- paste0("Spending time: ",
                      paste(round(x$settings$spendingTime, 3), collapse = ","))
      df1 <- data.frame(x = rep("", 13))
      colnames(df1) <- NULL
      rownames(df1) <- c(str1, str2, str3, str4, str5, str6, str7,
                         str8, str9, str10, paste(str11, str12, sep = ", "),
                         str13, "")
    } else {
      df1 <- data.frame(x = rep("", 12))
      colnames(df1) <- NULL
      rownames(df1) <- c(str1, str2, str3, str4, str5, str6, str7,
                         str8, str9, str10, paste(str11, str12, sep = ", "),
                         "")
    }
  } else {
    df1 <- data.frame(x = rep("", 11))
    colnames(df1) <- NULL
    rownames(df1) <- c(str1, str2, str3, str4, str5, str6, str7,
                       str8, str9, str10, "")
  }

  if (k>1) {
    b <- s[, c("informationRates", "efficacyBounds", "futilityBounds",
               "cumulativeRejection", "cumulativeFutility",
               "cumulativeAlphaSpent", "numberOfEvents",
               "numberOfDropouts", "numberOfSubjects", "numberOfMilestone",
               "analysisTime", "efficacySurvDiff", "futilitySurvDiff",
               "efficacyP", "futilityP", "information")]

    # format number of digits after decimal for each column
    j1 <- c(7,8,9,10,11)
    j2 <- 16
    j3 <- c(1,2,3,12,13)
    j4 <- c(4,5,6,14,15)

    b[j1] <- lapply(b[j1], formatC, format = "f", digits = 1)
    b[j2] <- lapply(b[j2], formatC, format = "f", digits = 2)
    b[j3] <- lapply(b[j3], formatC, format = "f", digits = 3)
    b[j4] <- lapply(b[j4], formatC, format = "f", digits = 4)

    if (x$settings$typeBetaSpending != 'none' ||
        (k > 1 && any(x$byStageResults$futilityBounds[1:(k-1)] > -8))) {
      df <- t(b)
      rownames(df) <- c("Information rate",
                        "Efficacy boundary (Z)",
                        "Futility boundary (Z)",
                        "Cumulative rejection",
                        "Cumulative futility",
                        "Cumulative alpha spent",
                        "Number of events",
                        "Number of dropouts",
                        "Number of subjects",
                        "Number of milestone subjects",
                        "Analysis time",
                        "Efficacy boundary (surv diff)",
                        "Futility boundary (surv diff)",
                        "Efficacy boundary (p)",
                        "Futility boundary (p)",
                        "Information")

    } else {
      df <- t(b[, c(1,2,4,6,7,8,9,10,11,12,14,16)])
      rownames(df) <- c("Information rate",
                        "Efficacy boundary (Z)",
                        "Cumulative rejection",
                        "Cumulative alpha spent",
                        "Number of events",
                        "Number of dropouts",
                        "Number of subjects",
                        "Number of milestone subjects",
                        "Analysis time",
                        "Efficacy boundary (surv diff)",
                        "Efficacy boundary (p)",
                        "Information")
    }

    colnames(df) <- paste("Stage", seq_len(ncol(df)), sep=" ")
  } else {

    b <- s[, c("efficacyBounds", "efficacySurvDiff", "efficacyP")]

    # format number of digits after decimal for each column
    j3 <- c(1,2)
    j4 <- 3

    b[j3] <- lapply(b[j3], formatC, format = "f", digits = 3)
    b[j4] <- lapply(b[j4], formatC, format = "f", digits = 4)

    df <- t(b)

    rownames(df) <- c("Efficacy boundary (Z)",
                      "Efficacy boundary (surv diff)",
                      "Efficacy boundary (p)")

    colnames(df) <- NA
  }

  print(df1, ..., na.print = "" , quote = FALSE )
  print(df, ..., na.print = "" , quote = FALSE )
  invisible(x)
}


#' @title Print Power and Sample Size Results for Restricted Mean Survival
#' Time Difference
#' @description Prints the summary statistics from power calculation.
#'
#' @param x The rmpower object to print.
#' @param ... Ensures that all arguments starting from "..." are named.
#'
#' @return A tabular printout of the summary statistics from power
#' calculation.
#'
#' @keywords internal
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @export
print.rmpower <- function(x, ...) {
  a <- x$overallResults
  s <- x$byStageResults
  k <- a$kMax

  if (k>1) {
    str1 <- paste0("Group-sequential design with ", k, " stages")
  } else {
    str1 <- "Fixed design"
  }

  str1 <- paste0(str1, " for difference in restricted mean survival time")

  str2 <- paste0("Milestone: ", round(a$milestone, 3), ", ",
                 "restricted mean survival time difference under H0: ",
                 round(a$rmstDiffH0, 3))

  str3 <- paste0("Restricted mean survival time on treatment: ",
                 round(a$rmst1, 3), ", ",
                 "on control: ", round(a$rmst2, 3))

  str4 <- paste0("Overall power: ",
                 round(a$overallReject, 4), ", ",
                 "overall significance level (1-sided): ",
                 round(a$alpha, 4))

  if (k>1) {
    str5 <- paste0("Maximum # events: ",
                   round(a$numberOfEvents, 1), ", ",
                   "expected # events: ",
                   round(a$expectedNumberOfEvents, 1))

    str6 <- paste0("Maximum # subjects: ",
                   round(a$numberOfSubjects, 1), ", ",
                   "expected # subjects: ",
                   round(a$expectedNumberOfSubjects, 1))

    str7 <- paste0("Maximum information: ",
                   round(a$information, 2), ", ",
                   "expected information: ",
                   round(a$expectedInformation, 2))

    str8 <- paste0("Total study duration: ",
                   round(a$studyDuration, 1), ", ",
                   "expected study duration: ",
                   round(a$expectedStudyDuration, 1))

  } else {
    str5 <- paste0("Number of events: ",
                   round(a$numberOfEvents, 1))

    str6 <- paste0("Number of subjects: ",
                   round(a$numberOfSubjects, 1))

    str7 <- paste0("Information: ",
                   round(a$information, 2))

    str8 <- paste0("Study duration: ",
                   round(a$studyDuration, 1))
  }

  str9 <- paste0("Accrual duration: ",
                 round(a$accrualDuration, 1), ", ",
                 "follow-up duration: ",
                 round(a$followupTime, 1), ", ",
                 "fixed follow-up: ", a$fixedFollowup)

  str10 <- paste0("Allocation ratio: ",
                  round(x$settings$allocationRatioPlanned, 3))

  if (k > 1) {
    asf <- tolower(x$settings$typeAlphaSpending)
    asfpar <- round(x$settings$parameterAlphaSpending, 3)
    asfuser <- round(x$settings$userAlphaSpending, 4)

    bsf <- tolower(x$settings$typeBetaSpending)
    bsfpar <- round(x$settings$parameterBetaSpending, 3)

    if (asf == "of") {
      str11 <- paste0("Alpha spending: O'Brien-Fleming")
    } else if (asf == "p") {
      str11 <- paste0("Alpha spending: Pocock")
    } else if (asf == "wt") {
      str11 <- paste0("Alpha spending: Wang-Tsiatis(Delta = ", asfpar, ")")
    } else if (asf == "sfof") {
      str11 <- paste0("Alpha spending: Lan-DeMets O'Brien-Fleming")
    } else if (asf == "sfp") {
      str11 <- paste0("Alpha spending: Lan-DeMets Pocock")
    } else if (asf == "sfkd") {
      str11 <- paste0("Alpha spending: KD(rho = ", asfpar, ")")
    } else if (asf == "sfhsd") {
      str11 <- paste0("Alpha spending: HSD(gamma = ", asfpar, ")")
    } else if (asf == "user") {
      str11 <- paste0("Alpha spending: User defined(",
                      paste(asfuser, collapse = ","), ")")
    } else {
      str11 <- "Alpha spending: None"
    }

    if (bsf == "of") {
      str12 <- paste0("beta spending: O'Brien-Fleming")
    } else if (bsf == "p") {
      str12 <- paste0("beta spending: Pocock")
    } else if (bsf == "wt") {
      str12 <- paste0("beta spending: Wang-Tsiatis(Delta = ", bsfpar, ")")
    } else if (bsf == "sfof") {
      str12 <- paste0("beta spending: Lan-DeMets O'Brien-Fleming")
    } else if (bsf == "sfp") {
      str12 <- paste0("beta spending: Lan-DeMets Pocock")
    } else if (bsf == "sfkd") {
      str12 <- paste0("beta spending: KD(rho = ", bsfpar, ")")
    } else if (bsf == "sfhsd") {
      str12 <- paste0("beta spending: HSD(gamma = ", bsfpar, ")")
    } else if (bsf == "user") {
      str12 <- paste0("beta spending: User defined")
    } else {
      str12 <- "beta spending: None"
    }

    if (!any(is.na(x$settings$spendingTime)) &&
        !all(x$settings$spendingTime == s$informationRates)) {
      str13 <- paste0("Spending time: ",
                      paste(round(x$settings$spendingTime, 3), collapse = ","))
      df1 <- data.frame(x = rep("", 13))
      colnames(df1) <- NULL
      rownames(df1) <- c(str1, str2, str3, str4, str5, str6, str7,
                         str8, str9, str10, paste(str11, str12, sep = ", "),
                         str13, "")
    } else {
      df1 <- data.frame(x = rep("", 12))
      colnames(df1) <- NULL
      rownames(df1) <- c(str1, str2, str3, str4, str5, str6, str7,
                         str8, str9, str10, paste(str11, str12, sep = ", "),
                         "")
    }
  } else {
    df1 <- data.frame(x = rep("", 11))
    colnames(df1) <- NULL
    rownames(df1) <- c(str1, str2, str3, str4, str5, str6, str7,
                       str8, str9, str10, "")
  }

  if (k>1) {
    b <- s[, c("informationRates", "efficacyBounds", "futilityBounds",
               "cumulativeRejection", "cumulativeFutility",
               "cumulativeAlphaSpent", "numberOfEvents",
               "numberOfDropouts", "numberOfSubjects", "numberOfMilestone",
               "analysisTime", "efficacyRmstDiff", "futilityRmstDiff",
               "efficacyP", "futilityP", "information")]

    # format number of digits after decimal for each column
    j1 <- c(7,8,9,10,11)
    j2 <- 16
    j3 <- c(1,2,3,12,13)
    j4 <- c(4,5,6,14,15)

    b[j1] <- lapply(b[j1], formatC, format = "f", digits = 1)
    b[j2] <- lapply(b[j2], formatC, format = "f", digits = 2)
    b[j3] <- lapply(b[j3], formatC, format = "f", digits = 3)
    b[j4] <- lapply(b[j4], formatC, format = "f", digits = 4)

    if (x$settings$typeBetaSpending != 'none' ||
        (k > 1 && any(x$byStageResults$futilityBounds[1:(k-1)] > -8))) {
      df <- t(b)
      rownames(df) <- c("Information rate",
                        "Efficacy boundary (Z)",
                        "Futility boundary (Z)",
                        "Cumulative rejection",
                        "Cumulative futility",
                        "Cumulative alpha spent",
                        "Number of events",
                        "Number of dropouts",
                        "Number of subjects",
                        "Number of milestone subjects",
                        "Analysis time",
                        "Efficacy boundary (rmst diff)",
                        "Futility boundary (rmst diff)",
                        "Efficacy boundary (p)",
                        "Futility boundary (p)",
                        "Information")

    } else {
      df <- t(b[, c(1,2,4,6,7,8,9,10,11,12,14,16)])
      rownames(df) <- c("Information rate",
                        "Efficacy boundary (Z)",
                        "Cumulative rejection",
                        "Cumulative alpha spent",
                        "Number of events",
                        "Number of dropouts",
                        "Number of subjects",
                        "Number of milestone subjects",
                        "Analysis time",
                        "Efficacy boundary (rmst diff)",
                        "Efficacy boundary (p)",
                        "Information")
    }

    colnames(df) <- paste("Stage", seq_len(ncol(df)), sep=" ")
  } else {

    b <- s[, c("efficacyBounds", "efficacyRmstDiff", "efficacyP")]

    # format number of digits after decimal for each column
    j3 <- c(1,2)
    j4 <- 3

    b[j3] <- lapply(b[j3], formatC, format = "f", digits = 3)
    b[j4] <- lapply(b[j4], formatC, format = "f", digits = 4)

    df <- t(b)

    rownames(df) <- c("Efficacy boundary (Z)",
                      "Efficacy boundary (rmst diff)",
                      "Efficacy boundary (p)")

    colnames(df) <- NA
  }

  print(df1, ..., na.print = "" , quote = FALSE )
  print(df, ..., na.print = "" , quote = FALSE )
  invisible(x)
}


#' @title Print Power and Sample Size Results for Equivalence in Milestone
#' Survival Probability Difference
#' @description Prints the summary statistics from power calculation of
#' equivalence in milestone survival probability difference.
#'
#' @param x The kmpowerequiv object to print.
#' @param ... Ensures that all arguments starting from "..." are named.
#'
#' @return A tabular printout of the summary statistics from power
#' calculation.
#'
#' @keywords internal
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @export
print.kmpowerequiv <- function(x, ...) {
  a <- x$overallResults
  s <- x$byStageResults
  k <- a$kMax

  if (k>1) {
    str1 <- paste0("Group-sequential design with ", k, " stages")
  } else {
    str1 <- "Fixed design"
  }

  str1 <- paste(str1, "for equivalence in milestone survival difference")

  str2 <- paste0("Milestone: ", round(a$milestone, 3), ", ",
                 "lower limit for survival difference: ",
                 round(a$survDiffLower, 3), ", ",
                 "upper limit: ",
                 round(a$survDiffUpper, 3))

  str3 <- paste0("Milestone survival on treatment: ",
                 round(a$surv1, 3), ", ",
                 "on control: ", round(a$surv2, 3), ", ",
                 "difference: ", round(a$survDiff, 3))

  str4 <- paste0("Overall power: ",
                 round(a$overallReject, 4), ", ",
                 "overall alpha: ",
                 round(a$alpha, 4))

  if (k>1) {
    str5 <- paste0("Maximum # events: ",
                   round(a$numberOfEvents, 1), ", ",
                   "expected # events: ",
                   round(a$expectedNumberOfEvents, 1))

    str6 <- paste0("Maximum # subjects: ",
                   round(a$numberOfSubjects, 1), ", ",
                   "expected # subjects: ",
                   round(a$expectedNumberOfSubjects, 1))

    str7 <- paste0("Maximum information: ",
                   round(a$information, 2), ", ",
                   "expected information: ",
                   round(a$expectedInformation, 2))

    str8 <- paste0("Total study duration: ",
                   round(a$studyDuration, 1), ", ",
                   "expected study duration: ",
                   round(a$expectedStudyDuration, 1))
  } else {
    str5 <- paste0("Number of events: ",
                   round(a$numberOfEvents, 1))

    str6 <- paste0("Number of subjects: ",
                   round(a$numberOfSubjects, 1))

    str7 <- paste0("Information: ",
                   round(a$information, 2))

    str8 <- paste0("Study duration: ",
                   round(a$studyDuration, 1))
  }

  str9 <- paste0("Accrual duration: ",
                 round(a$accrualDuration, 1), ", ",
                 "follow-up duration: ",
                 round(a$followupTime, 1), ", ",
                 "fixed follow-up: ", a$fixedFollowup)

  str10 <- paste0("Allocation ratio: ",
                  round(x$settings$allocationRatioPlanned, 3))

  if (k > 1) {
    asf <- tolower(x$settings$typeAlphaSpending)
    asfpar <- round(x$settings$parameterAlphaSpending, 3)
    asfuser <- round(x$settings$userAlphaSpending, 4)

    if (asf == "of") {
      str11 <- paste0("Alpha spending: O'Brien-Fleming")
    } else if (asf == "p") {
      str11 <- paste0("Alpha spending: Pocock")
    } else if (asf == "wt") {
      str11 <- paste0("Alpha spending: Wang-Tsiatis(Delta = ", asfpar, ")")
    } else if (asf == "sfof") {
      str11 <- paste0("Alpha spending: Lan-DeMets O'Brien-Fleming")
    } else if (asf == "sfp") {
      str11 <- paste0("Alpha spending: Lan-DeMets Pocock")
    } else if (asf == "sfkd") {
      str11 <- paste0("Alpha spending: KD(rho = ", asfpar, ")")
    } else if (asf == "sfhsd") {
      str11 <- paste0("Alpha spending: HSD(gamma = ", asfpar, ")")
    } else if (asf == "user") {
      str11 <- paste0("Alpha spending: User defined(",
                      paste(asfuser, collapse = ","), ")")
    } else {
      str11 <- "Alpha spending: None"
    }

    if (!any(is.na(x$settings$spendingTime)) &&
        !all(x$settings$spendingTime == s$informationRates)) {
      str12 <- paste0("Spending time: ",
                      paste(round(x$settings$spendingTime, 3), collapse = ","))
      df1 <- data.frame(x = rep("", 13))
      colnames(df1) <- NULL
      rownames(df1) <- c(str1, str2, str3, str4, str5, str6, str7,
                         str8, str9, str10, str11, str12, "")
    } else {
      df1 <- data.frame(x = rep("", 12))
      colnames(df1) <- NULL
      rownames(df1) <- c(str1, str2, str3, str4, str5, str6, str7,
                         str8, str9, str10, str11, "")
    }
  } else {
    df1 <- data.frame(x = rep("", 11))
    colnames(df1) <- NULL
    rownames(df1) <- c(str1, str2, str3, str4, str5, str6, str7,
                       str8, str9, str10, "")
  }

  if (k>1) {
    b <- s[, c("informationRates", "efficacyBounds",
               "cumulativeRejection", "cumulativeAlphaSpent",
               "cumulativeAttainedAlphaH10",
               "cumulativeAttainedAlphaH20", "numberOfEvents",
               "numberOfDropouts", "numberOfSubjects",
               "numberOfMilestone", "analysisTime",
               "efficacySurvDiffLower", "efficacySurvDiffUpper",
               "efficacyP", "information")]

    # format number of digits after decimal for each column
    j1 <- c(7,8,9,10,11)
    j2 <- 15
    j3 <- c(1,2,12,13)
    j4 <- c(3,4,5,6,14)

    b[j1] <- lapply(b[j1], formatC, format = "f", digits = 1)
    b[j2] <- lapply(b[j2], formatC, format = "f", digits = 2)
    b[j3] <- lapply(b[j3], formatC, format = "f", digits = 3)
    b[j4] <- lapply(b[j4], formatC, format = "f", digits = 4)

    df <- t(b)
    rownames(df) <- c("Information rate",
                      "Boundary for each 1-sided test (Z)",
                      "Cumulative rejection",
                      "Cumulative alpha for each 1-sided test",
                      "Cumulative alpha attained under H10",
                      "Cumulative alpha attained under H20",
                      "Number of events",
                      "Number of dropouts",
                      "Number of subjects",
                      "Number of milestone subjects",
                      "Analysis time",
                      "Boundary for lower limit (surv diff)",
                      "Boundary for upper limit (surv diff)",
                      "Boundary for each 1-sided test (p)",
                      "Information")

    colnames(df) <- paste("Stage", seq_len(ncol(df)), sep=" ")
  } else {
    b <- s[, c("efficacyBounds", "efficacySurvDiffLower",
               "efficacySurvDiffUpper",  "efficacyP")]

    # format number of digits after decimal for each column
    j3 <- c(1,2,3)
    j4 <- 4

    b[j3] <- lapply(b[j3], formatC, format = "f", digits = 3)
    b[j4] <- lapply(b[j4], formatC, format = "f", digits = 4)

    df <- t(b)

    rownames(df) <- c("Boundary for each 1-sided test (Z)",
                      "Boundary for lower limit (surv diff)",
                      "Boundary for upper limit (surv diff)",
                      "Boundary for each 1-sided test (p)")
    colnames(df) <- NA
  }

  print(df1, ..., na.print = "" , quote = FALSE )
  print(df, ..., na.print = "" , quote = FALSE )
  invisible(x)
}


#' @title Print Power and Sample Size Results for Equivalence in Restricted
#' Mean Survival Time Difference
#' @description Prints the summary statistics from power calculation of
#' equivalence in restricted mean survival time difference.
#'
#' @param x The rmpowerequiv object to print.
#' @param ... Ensures that all arguments starting from "..." are named.
#'
#' @return A tabular printout of the summary statistics from power
#' calculation.
#'
#' @keywords internal
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @export
print.rmpowerequiv <- function(x, ...) {
  a <- x$overallResults
  s <- x$byStageResults
  k <- a$kMax

  if (k>1) {
    str1 <- paste0("Group-sequential design with ", k, " stages")
  } else {
    str1 <- "Fixed design"
  }

  str1 <- paste(str1, paste("for equivalence in RMST difference"))

  str2 <- paste0("Milestone: ", round(a$milestone, 3), ", ",
                 "lower limit for RMST difference: ",
                 round(a$rmstDiffLower, 3), ", ",
                 "upper limit: ",
                 round(a$rmstDiffUpper, 3))

  str3 <- paste0("RMST on treatment: ",
                 round(a$rmst1, 3), ", ",
                 "on control: ", round(a$rmst2, 3), ", ",
                 "difference: ", round(a$rmstDiff, 3))

  str4 <- paste0("Overall power: ",
                 round(a$overallReject, 4), ", ",
                 "overall alpha: ",
                 round(a$alpha, 4))

  if (k>1) {
    str5 <- paste0("Maximum # events: ",
                   round(a$numberOfEvents, 1), ", ",
                   "expected # events: ",
                   round(a$expectedNumberOfEvents, 1))

    str6 <- paste0("Maximum # subjects: ",
                   round(a$numberOfSubjects, 1), ", ",
                   "expected # subjects: ",
                   round(a$expectedNumberOfSubjects, 1))

    str7 <- paste0("Maximum information: ",
                   round(a$information, 2), ", ",
                   "expected information: ",
                   round(a$expectedInformation, 2))

    str8 <- paste0("Total study duration: ",
                   round(a$studyDuration, 1), ", ",
                   "expected study duration: ",
                   round(a$expectedStudyDuration, 1))
  } else {
    str5 <- paste0("Number of events: ",
                   round(a$numberOfEvents, 1))

    str6 <- paste0("Number of subjects: ",
                   round(a$numberOfSubjects, 1))

    str7 <- paste0("Information: ",
                   round(a$information, 2))

    str8 <- paste0("Study duration: ",
                   round(a$studyDuration, 1))
  }

  str9 <- paste0("Accrual duration: ",
                 round(a$accrualDuration, 1), ", ",
                 "follow-up duration: ",
                 round(a$followupTime, 1), ", ",
                 "fixed follow-up: ", a$fixedFollowup)

  str10 <- paste0("Allocation ratio: ",
                  round(x$settings$allocationRatioPlanned, 3))

  if (k > 1) {
    asf <- tolower(x$settings$typeAlphaSpending)
    asfpar <- round(x$settings$parameterAlphaSpending, 3)
    asfuser <- round(x$settings$userAlphaSpending, 4)

    if (asf == "of") {
      str11 <- paste0("Alpha spending: O'Brien-Fleming")
    } else if (asf == "p") {
      str11 <- paste0("Alpha spending: Pocock")
    } else if (asf == "wt") {
      str11 <- paste0("Alpha spending: Wang-Tsiatis(Delta = ", asfpar, ")")
    } else if (asf == "sfof") {
      str11 <- paste0("Alpha spending: Lan-DeMets O'Brien-Fleming")
    } else if (asf == "sfp") {
      str11 <- paste0("Alpha spending: Lan-DeMets Pocock")
    } else if (asf == "sfkd") {
      str11 <- paste0("Alpha spending: KD(rho = ", asfpar, ")")
    } else if (asf == "sfhsd") {
      str11 <- paste0("Alpha spending: HSD(gamma = ", asfpar, ")")
    } else if (asf == "user") {
      str11 <- paste0("Alpha spending: User defined(",
                      paste(asfuser, collapse = ","), ")")
    } else {
      str11 <- "Alpha spending: None"
    }

    if (!any(is.na(x$settings$spendingTime)) &&
        !all(x$settings$spendingTime == s$informationRates)) {
      str12 <- paste0("Spending time: ",
                      paste(round(x$settings$spendingTime, 3), collapse = ","))
      df1 <- data.frame(x = rep("", 13))
      colnames(df1) <- NULL
      rownames(df1) <- c(str1, str2, str3, str4, str5, str6, str7,
                         str8, str9, str10, str11, str12, "")
    } else {
      df1 <- data.frame(x = rep("", 12))
      colnames(df1) <- NULL
      rownames(df1) <- c(str1, str2, str3, str4, str5, str6, str7,
                         str8, str9, str10, str11, "")
    }
  } else {
    df1 <- data.frame(x = rep("", 11))
    colnames(df1) <- NULL
    rownames(df1) <- c(str1, str2, str3, str4, str5, str6, str7,
                       str8, str9, str10, "")
  }

  if (k>1) {
    b <- s[, c("informationRates", "efficacyBounds",
               "cumulativeRejection", "cumulativeAlphaSpent",
               "cumulativeAttainedAlphaH10",
               "cumulativeAttainedAlphaH20", "numberOfEvents",
               "numberOfDropouts", "numberOfSubjects",
               "numberOfMilestone", "analysisTime",
               "efficacyRmstDiffLower", "efficacyRmstDiffUpper",
               "efficacyP", "information")]

    # format number of digits after decimal for each column
    j1 <- c(7,8,9,10,11)
    j2 <- 15
    j3 <- c(1,2,12,13)
    j4 <- c(3,4,5,6,14)

    b[j1] <- lapply(b[j1], formatC, format = "f", digits = 1)
    b[j2] <- lapply(b[j2], formatC, format = "f", digits = 2)
    b[j3] <- lapply(b[j3], formatC, format = "f", digits = 3)
    b[j4] <- lapply(b[j4], formatC, format = "f", digits = 4)

    df <- t(b)
    rownames(df) <- c("Information rate",
                      "Boundary for each 1-sided test (Z)",
                      "Cumulative rejection",
                      "Cumulative alpha for each 1-sided test",
                      "Cumulative alpha attained under H10",
                      "Cumulative alpha attained under H20",
                      "Number of events",
                      "Number of dropouts",
                      "Number of subjects",
                      "Number of milestone subjects",
                      "Analysis time",
                      "Boundary for lower limit (rmst diff)",
                      "Boundary for upper limit (rmst diff)",
                      "Boundary for each 1-sided test (p)",
                      "Information")

    colnames(df) <- paste("Stage", seq_len(ncol(df)), sep=" ")
  } else {
    b <- s[, c("efficacyBounds", "efficacyRmstDiffLower",
               "efficacyRmstDiffUpper",  "efficacyP")]

    # format number of digits after decimal for each column
    j3 <- c(1,2,3)
    j4 <- 4

    b[j3] <- lapply(b[j3], formatC, format = "f", digits = 3)
    b[j4] <- lapply(b[j4], formatC, format = "f", digits = 4)

    df <- t(b)

    rownames(df) <- c("Boundary for each 1-sided test (Z)",
                      "Boundary for lower limit (rmst diff)",
                      "Boundary for upper limit (rmst diff)",
                      "Boundary for each 1-sided test (p)")
    colnames(df) <- NA
  }

  print(df1, ..., na.print = "" , quote = FALSE )
  print(df, ..., na.print = "" , quote = FALSE )
  invisible(x)
}


#' @title Print Power and Sample Size Results for Equivalence in Hazard
#' Ratio
#' @description Prints the summary statistics from power calculation of
#' equivalence in hazard ratio.
#'
#' @param x The lrpowerequiv object to print.
#' @param ... Ensures that all arguments starting from "..." are named.
#'
#' @return A tabular printout of the summary statistics from power
#' calculation.
#'
#' @keywords internal
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @export
print.lrpowerequiv <- function(x, ...) {
  a <- x$overallResults
  s <- x$byStageResults
  k <- a$kMax

  if (k>1) {
    str1 <- paste0("Group-sequential design with ", k, " stages")
  } else {
    str1 <- "Fixed design"
  }

  str1 <- paste(str1, "for equivalence in hazard ratio")

  str2 <- paste0("Lower limit for hazard ratio: ",
                 round(a$hazardRatioLower, 3), ", ",
                 "upper limit for hazard ratio: ",
                 round(a$hazardRatioUpper, 3))

  str3 <- paste0("Overall power: ",
                 round(a$overallReject, 4), ", ",
                 "overall alpha: ",
                 round(a$alpha, 4))

  if (k>1) {
    str4 <- paste0("Maximum # events: ",
                   round(a$numberOfEvents, 1), ", ",
                   "expected # events: ",
                   round(a$expectedNumberOfEvents, 1))

    str5 <- paste0("Maximum # dropouts: ",
                   round(a$numberOfDropouts, 1), ", ",
                   "expected # dropouts: ",
                   round(a$expectedNumberOfDropouts, 1))

    str6 <- paste0("Maximum # subjects: ",
                   round(a$numberOfSubjects, 1), ", ",
                   "expected # subjects: ",
                   round(a$expectedNumberOfSubjects, 1))

    str7 <- paste0("Maximum information: ",
                   round(a$information, 2), ", ",
                   "expected information: ",
                   round(a$expectedInformation, 2))

    str8 <- paste0("Total study duration: ",
                   round(a$studyDuration, 1), ", ",
                   "expected study duration: ",
                   round(a$expectedStudyDuration, 1))
  } else {
    str4 <- paste0("Number of events: ",
                   round(a$numberOfEvents, 1))

    str5 <- paste0("Number of dropouts: ",
                   round(a$numberOfDropouts, 1))

    str6 <- paste0("Number of subjects: ",
                   round(a$numberOfSubjects, 1))

    str7 <- paste0("Information: ",
                   round(a$information, 2))

    str8 <- paste0("Study duration: ",
                   round(a$studyDuration, 1))
  }

  str9 <- paste0("Accrual duration: ",
                 round(a$accrualDuration, 1), ", ",
                 "follow-up duration: ",
                 round(a$followupTime, 1), ", ",
                 "fixed follow-up: ", a$fixedFollowup)

  str10 <- paste0("Allocation ratio: ",
                  round(x$settings$allocationRatioPlanned, 3))

  if (k > 1) {
    asf <- tolower(x$settings$typeAlphaSpending)
    asfpar <- round(x$settings$parameterAlphaSpending, 3)
    asfuser <- round(x$settings$userAlphaSpending, 4)

    if (asf == "of") {
      str11 <- paste0("Alpha spending: O'Brien-Fleming")
    } else if (asf == "p") {
      str11 <- paste0("Alpha spending: Pocock")
    } else if (asf == "wt") {
      str11 <- paste0("Alpha spending: Wang-Tsiatis(Delta = ", asfpar, ")")
    } else if (asf == "sfof") {
      str11 <- paste0("Alpha spending: Lan-DeMets O'Brien-Fleming")
    } else if (asf == "sfp") {
      str11 <- paste0("Alpha spending: Lan-DeMets Pocock")
    } else if (asf == "sfkd") {
      str11 <- paste0("Alpha spending: KD(rho = ", asfpar, ")")
    } else if (asf == "sfhsd") {
      str11 <- paste0("Alpha spending: HSD(gamma = ", asfpar, ")")
    } else if (asf == "user") {
      str11 <- paste0("Alpha spending: User defined(",
                      paste(asfuser, collapse = ","), ")")
    } else {
      str11 <- "Alpha spending: None"
    }

    if (!any(is.na(x$settings$spendingTime)) &&
        !all(x$settings$spendingTime == s$informationRates)) {
      str12 <- paste0("Spending time: ",
                      paste(round(x$settings$spendingTime, 3), collapse = ","))
      df1 <- data.frame(x = rep("", 13))
      colnames(df1) <- NULL
      rownames(df1) <- c(str1, str2, str3, str4, str5, str6, str7,
                         str8, str9, str10, str11, str12, "")
    } else {
      df1 <- data.frame(x = rep("", 12))
      colnames(df1) <- NULL
      rownames(df1) <- c(str1, str2, str3, str4, str5, str6, str7,
                         str8, str9, str10, str11, "")
    }
  } else {
    df1 <- data.frame(x = rep("", 11))
    colnames(df1) <- NULL
    rownames(df1) <- c(str1, str2, str3, str4, str5, str6, str7,
                       str8, str9, str10, "")
  }

  if (k>1) {
    b <- s[, c("informationRates", "efficacyBounds",
               "cumulativeRejection", "cumulativeAlphaSpent",
               "cumulativeAttainedAlphaH10",
               "cumulativeAttainedAlphaH20", "numberOfEvents",
               "numberOfDropouts", "numberOfSubjects",
               "analysisTime", "efficacyHRLower", "efficacyHRUpper",
               "efficacyP", "information", "HR")]

    # format number of digits after decimal for each column
    j1 <- c(7,8,9,10)
    j2 <- c(14,15)
    j3 <- c(1,2,11,12)
    j4 <- c(3,4,5,6,13)

    b[j1] <- lapply(b[j1], formatC, format = "f", digits = 1)
    b[j2] <- lapply(b[j2], formatC, format = "f", digits = 2)
    b[j3] <- lapply(b[j3], formatC, format = "f", digits = 3)
    b[j4] <- lapply(b[j4], formatC, format = "f", digits = 4)

    df <- t(b)
    rownames(df) <- c("Information rate",
                      "Boundary for each 1-sided test (Z)",
                      "Cumulative rejection",
                      "Cumulative alpha for each 1-sided test",
                      "Cumulative alpha attained under H10",
                      "Cumulative alpha attained under H20",
                      "Number of events",
                      "Number of dropouts",
                      "Number of subjects",
                      "Analysis time",
                      "Boundary for lower limit (HR)",
                      "Boundary for upper limit (HR)",
                      "Boundary for each 1-sided test (p)",
                      "Information",
                      "HR")

    colnames(df) <- paste("Stage", seq_len(ncol(df)), sep=" ")
  } else {
    b <- s[, c("efficacyBounds", "efficacyHRLower",
               "efficacyHRUpper",  "efficacyP")]

    # format number of digits after decimal for each column
    j3 <- c(1,2,3)
    j4 <- 4

    b[j3] <- lapply(b[j3], formatC, format = "f", digits = 3)
    b[j4] <- lapply(b[j4], formatC, format = "f", digits = 4)

    df <- t(b)

    rownames(df) <- c("Boundary for each 1-sided test (Z)",
                      "Boundary for lower limit (HR)",
                      "Boundary for upper limit (HR)",
                      "Boundary for each 1-sided test (p)")
    colnames(df) <- NA
  }

  print(df1, ..., na.print = "" , quote = FALSE )
  print(df, ..., na.print = "" , quote = FALSE )
  invisible(x)
}


#' @title Print Power and Sample Size Results for One-Sample Milestone
#' Survival Probability
#' @description Prints the summary statistics from power calculation of
#' one-sample milestone survival probability.
#'
#' @param x The kmpower1s object to print.
#' @param ... Ensures that all arguments starting from "..." are named.
#'
#' @return A tabular printout of the summary statistics from power
#' calculation.
#'
#' @keywords internal
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @export
print.kmpower1s <- function(x, ...) {
  a <- x$overallResults
  s <- x$byStageResults
  k <- a$kMax

  if (k>1) {
    str1 <- paste0("Group-sequential design with ", k, " stages")
  } else {
    str1 <- "Fixed design"
  }

  str1 <- paste(str1, "for one-sample milestone survival probability")

  str2 <- paste0("Milestone: ", round(a$milestone, 3), ", ",
                 "survival probability under H0: ",
                 round(a$survH0, 3), ", ",
                 "under H1: ", round(a$surv, 3))

  str3 <- paste0("Overall power: ",
                 round(a$overallReject, 4), ", ",
                 "overall significance level (1-sided): ",
                 round(a$alpha, 4))

  if (k>1) {
    str4 <- paste0("Maximum # events: ",
                   round(a$numberOfEvents, 1), ", ",
                   "expected # events: ",
                   round(a$expectedNumberOfEvents, 1))

    str5 <- paste0("Maximum # subjects: ",
                   round(a$numberOfSubjects, 1), ", ",
                   "expected # subjects: ",
                   round(a$expectedNumberOfSubjects, 1))

    str6 <- paste0("Maximum # milestone subjects: ",
                   round(a$numberOfMilestone, 1), ", ",
                   "expected # milestone subjects: ",
                   round(a$expectedNumberOfMilestone, 1))

    str7 <- paste0("Maximum information: ",
                   round(a$information, 2), ", ",
                   "expected information: ",
                   round(a$expectedInformation, 2))

    str8 <- paste0("Total study duration: ",
                   round(a$studyDuration, 1), ", ",
                   "expected study duration: ",
                   round(a$expectedStudyDuration, 1))
  } else {
    str4 <- paste0("Number of events: ",
                   round(a$numberOfEvents, 1))

    str5 <- paste0("Number of subjects: ",
                   round(a$numberOfSubjects, 1))

    str6 <- paste0("Number of milestone subjects: ",
                   round(a$numberOfMilestone, 1))

    str7 <- paste0("Information: ",
                   round(a$information, 2))

    str8 <- paste0("Study duration: ",
                   round(a$studyDuration, 1))
  }

  str9 <- paste0("Accrual duration: ",
                 round(a$accrualDuration, 1), ", ",
                 "follow-up duration: ",
                 round(a$followupTime, 1), ", ",
                 "fixed follow-up: ", a$fixedFollowup)


  if (k > 1) {
    asf <- tolower(x$settings$typeAlphaSpending)
    asfpar <- round(x$settings$parameterAlphaSpending, 3)
    asfuser <- round(x$settings$userAlphaSpending, 4)

    bsf <- tolower(x$settings$typeBetaSpending)
    bsfpar <- round(x$settings$parameterBetaSpending, 3)

    if (asf == "of") {
      str10 <- paste0("Alpha spending: O'Brien-Fleming")
    } else if (asf == "p") {
      str10 <- paste0("Alpha spending: Pocock")
    } else if (asf == "wt") {
      str10 <- paste0("Alpha spending: Wang-Tsiatis(Delta = ", asfpar, ")")
    } else if (asf == "sfof") {
      str10 <- paste0("Alpha spending: Lan-DeMets O'Brien-Fleming")
    } else if (asf == "sfp") {
      str10 <- paste0("Alpha spending: Lan-DeMets Pocock")
    } else if (asf == "sfkd") {
      str10 <- paste0("Alpha spending: KD(rho = ", asfpar, ")")
    } else if (asf == "sfhsd") {
      str10 <- paste0("Alpha spending: HSD(gamma = ", asfpar, ")")
    } else if (asf == "user") {
      str10 <- paste0("Alpha spending: User defined(",
                      paste(asfuser, collapse = ","), ")")
    } else {
      str10 <- "Alpha spending: None"
    }

    if (bsf == "of") {
      str11 <- paste0("beta spending: O'Brien-Fleming")
    } else if (bsf == "p") {
      str11 <- paste0("beta spending: Pocock")
    } else if (bsf == "wt") {
      str11 <- paste0("beta spending: Wang-Tsiatis(Delta = ", bsfpar, ")")
    } else if (bsf == "sfof") {
      str11 <- paste0("beta spending: Lan-DeMets O'Brien-Fleming")
    } else if (bsf == "sfp") {
      str11 <- paste0("beta spending: Lan-DeMets Pocock")
    } else if (bsf == "sfkd") {
      str11 <- paste0("beta spending: KD(rho = ", bsfpar, ")")
    } else if (bsf == "sfhsd") {
      str11 <- paste0("beta spending: HSD(gamma = ", bsfpar, ")")
    } else if (bsf == "user") {
      bsfuser <- round(x$settings$userBetaSpending, 4)
      str11 <- paste0("beta spending: User defined(",
                      paste(bsfuser, collapse = ","), ")")
    } else {
      str11 <- "beta spending: None"
    }

    if (!any(is.na(x$settings$spendingTime)) &&
        !all(x$settings$spendingTime == s$informationRates)) {
      str12 <- paste0("Spending time: ",
                      paste(round(x$settings$spendingTime, 3), collapse = ","))
      df1 <- data.frame(x = rep("", 12))
      colnames(df1) <- NULL
      rownames(df1) <- c(str1, str2, str3, str4, str5, str6, str7,
                         str8, str9, paste(str10, str11, sep = ", "),
                         str12, "")
    } else {
      df1 <- data.frame(x = rep("", 11))
      colnames(df1) <- NULL
      rownames(df1) <- c(str1, str2, str3, str4, str5, str6, str7,
                         str8, str9, paste(str10, str11, sep = ", "), "")
    }
  } else {
    df1 <- data.frame(x = rep("", 10))
    colnames(df1) <- NULL
    rownames(df1) <- c(str1, str2, str3, str4, str5, str6, str7,
                       str8, str9, "")
  }

  if (k>1) {
    b <- s[, c("informationRates", "efficacyBounds", "futilityBounds",
               "cumulativeRejection", "cumulativeFutility",
               "cumulativeAlphaSpent", "numberOfEvents",
               "numberOfDropouts", "numberOfSubjects", "numberOfMilestone",
               "analysisTime", "efficacySurv", "futilitySurv",
               "efficacyP", "futilityP", "information")]

    # format number of digits after decimal for each column
    j1 <- c(7,8,9,10,11)
    j2 <- 16
    j3 <- c(1,2,3,12,13)
    j4 <- c(4,5,6,14,15)

    b[j1] <- lapply(b[j1], formatC, format = "f", digits = 1)
    b[j2] <- lapply(b[j2], formatC, format = "f", digits = 2)
    b[j3] <- lapply(b[j3], formatC, format = "f", digits = 3)
    b[j4] <- lapply(b[j4], formatC, format = "f", digits = 4)

    if (x$settings$typeBetaSpending != 'none' ||
        (k > 1 && any(x$byStageResults$futilityBounds[1:(k-1)] > -8))) {
      df <- t(b)
      rownames(df) <- c("Information rate",
                        "Efficacy boundary (Z)",
                        "Futility boundary (Z)",
                        "Cumulative rejection",
                        "Cumulative futility",
                        "Cumulative alpha spent",
                        "Number of events",
                        "Number of dropouts",
                        "Number of subjects",
                        "Number of milestone subjects",
                        "Analysis time",
                        "Efficacy boundary (surv)",
                        "Futility boundary (surv)",
                        "Efficacy boundary (p)",
                        "Futility boundary (p)",
                        "Information")

    } else {
      df <- t(b[, c(1,2,4,6,7,8,9,10,11,12,14,16)])
      rownames(df) <- c("Information rate",
                        "Efficacy boundary (Z)",
                        "Cumulative rejection",
                        "Cumulative alpha spent",
                        "Number of events",
                        "Number of dropouts",
                        "Number of subjects",
                        "Number of milestone subjects",
                        "Analysis time",
                        "Efficacy boundary (surv)",
                        "Efficacy boundary (p)",
                        "Information")
    }

    colnames(df) <- paste("Stage", seq_len(ncol(df)), sep=" ")
  } else {
    b <- s[, c("efficacyBounds", "efficacySurv", "efficacyP")]

    # format number of digits after decimal for each column
    j3 <- c(1,2)
    j4 <- 3

    b[j3] <- lapply(b[j3], formatC, format = "f", digits = 3)
    b[j4] <- lapply(b[j4], formatC, format = "f", digits = 4)

    df <- t(b)

    rownames(df) <- c("Efficacy boundary (Z)",
                      "Efficacy boundary (surv)",
                      "Efficacy boundary (p)")
    colnames(df) <- NA
  }

  print(df1, ..., na.print = "" , quote = FALSE )
  print(df, ..., na.print = "" , quote = FALSE )
  invisible(x)
}

#' @title Print Power and Sample Size Results for One-Sample Restricted
#' Mean Survival Time
#' @description Prints the summary statistics from power calculation of
#' one-sample restricted mean survival time.
#'
#' @param x The rmpower1s object to print.
#' @param ... Ensures that all arguments starting from "..." are named.
#'
#' @return A tabular printout of the summary statistics from power
#' calculation.
#'
#' @keywords internal
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @export
print.rmpower1s <- function(x, ...) {
  a <- x$overallResults
  s <- x$byStageResults
  k <- a$kMax

  if (k>1) {
    str1 <- paste0("Group-sequential design with ", k, " stages")
  } else {
    str1 <- "Fixed design"
  }

  str1 <- paste(str1, "for one-sample restricted mean survival time")


  str2 <- paste0("Milestone: ", round(a$milestone, 3), ", ",
                 "restricted mean survival time under H0: ",
                 round(a$rmstH0, 3), ", ",
                 "under H1: ", round(a$rmst, 3))

  str3 <- paste0("Overall power: ",
                 round(a$overallReject, 4), ", ",
                 "overall significance level (1-sided): ",
                 round(a$alpha, 4))

  if (k>1) {
    str4 <- paste0("Maximum # events: ",
                   round(a$numberOfEvents, 1), ", ",
                   "expected # events: ",
                   round(a$expectedNumberOfEvents, 1))

    str5 <- paste0("Maximum # subjects: ",
                   round(a$numberOfSubjects, 1), ", ",
                   "expected # subjects: ",
                   round(a$expectedNumberOfSubjects, 1))

    str6 <- paste0("Maximum # milestone subjects: ",
                   round(a$numberOfMilestone, 1), ", ",
                   "expected # milestone subjects: ",
                   round(a$expectedNumberOfMilestone, 1))

    str7 <- paste0("Maximum information: ",
                   round(a$information, 2), ", ",
                   "expected information: ",
                   round(a$expectedInformation, 2))

    str8 <- paste0("Total study duration: ",
                   round(a$studyDuration, 1), ", ",
                   "expected study duration: ",
                   round(a$expectedStudyDuration, 1))
  } else {
    str4 <- paste0("Number of events: ",
                   round(a$numberOfEvents, 1))

    str5 <- paste0("Number of subjects: ",
                   round(a$numberOfSubjects, 1))

    str6 <- paste0("Number of milestone subjects: ",
                   round(a$numberOfMilestone, 1))

    str7 <- paste0("Information: ",
                   round(a$information, 2))

    str8 <- paste0("Study duration: ",
                   round(a$studyDuration, 1))
  }

  str9 <- paste0("Accrual duration: ",
                 round(a$accrualDuration, 1), ", ",
                 "follow-up duration: ",
                 round(a$followupTime, 1), ", ",
                 "fixed follow-up: ", a$fixedFollowup)


  if (k > 1) {
    asf <- tolower(x$settings$typeAlphaSpending)
    asfpar <- round(x$settings$parameterAlphaSpending, 3)
    asfuser <- round(x$settings$userAlphaSpending, 4)

    bsf <- tolower(x$settings$typeBetaSpending)
    bsfpar <- round(x$settings$parameterBetaSpending, 3)

    if (asf == "of") {
      str10 <- paste0("Alpha spending: O'Brien-Fleming")
    } else if (asf == "p") {
      str10 <- paste0("Alpha spending: Pocock")
    } else if (asf == "wt") {
      str10 <- paste0("Alpha spending: Wang-Tsiatis(Delta = ", asfpar, ")")
    } else if (asf == "sfof") {
      str10 <- paste0("Alpha spending: Lan-DeMets O'Brien-Fleming")
    } else if (asf == "sfp") {
      str10 <- paste0("Alpha spending: Lan-DeMets Pocock")
    } else if (asf == "sfkd") {
      str10 <- paste0("Alpha spending: KD(rho = ", asfpar, ")")
    } else if (asf == "sfhsd") {
      str10 <- paste0("Alpha spending: HSD(gamma = ", asfpar, ")")
    } else if (asf == "user") {
      str10 <- paste0("Alpha spending: User defined(",
                      paste(asfuser, collapse = ","), ")")
    } else {
      str10 <- "Alpha spending: None"
    }

    if (bsf == "of") {
      str11 <- paste0("beta spending: O'Brien-Fleming")
    } else if (bsf == "p") {
      str11 <- paste0("beta spending: Pocock")
    } else if (bsf == "wt") {
      str11 <- paste0("beta spending: Wang-Tsiatis(Delta = ", bsfpar, ")")
    } else if (bsf == "sfof") {
      str11 <- paste0("beta spending: Lan-DeMets O'Brien-Fleming")
    } else if (bsf == "sfp") {
      str11 <- paste0("beta spending: Lan-DeMets Pocock")
    } else if (bsf == "sfkd") {
      str11 <- paste0("beta spending: KD(rho = ", bsfpar, ")")
    } else if (bsf == "sfhsd") {
      str11 <- paste0("beta spending: HSD(gamma = ", bsfpar, ")")
    } else if (bsf == "user") {
      str11 <- paste0("beta spending: User defined")
    } else {
      str11 <- "beta spending: None"
    }

    if (!any(is.na(x$settings$spendingTime)) &&
        !all(x$settings$spendingTime == s$informationRates)) {
      str12 <- paste0("Spending time: ",
                      paste(round(x$settings$spendingTime, 3), collapse = ","))
      df1 <- data.frame(x = rep("", 12))
      colnames(df1) <- NULL
      rownames(df1) <- c(str1, str2, str3, str4, str5, str6, str7,
                         str8, str9, paste(str10, str11, sep = ", "),
                         str12, "")
    } else {
      df1 <- data.frame(x = rep("", 11))
      colnames(df1) <- NULL
      rownames(df1) <- c(str1, str2, str3, str4, str5, str6, str7,
                         str8, str9, paste(str10, str11, sep = ", "), "")
    }
  } else {
    df1 <- data.frame(x = rep("", 10))
    colnames(df1) <- NULL
    rownames(df1) <- c(str1, str2, str3, str4, str5, str6, str7,
                       str8, str9, "")
  }

  if (k>1) {
    b <- s[, c("informationRates", "efficacyBounds", "futilityBounds",
               "cumulativeRejection", "cumulativeFutility",
               "cumulativeAlphaSpent", "numberOfEvents",
               "numberOfDropouts", "numberOfSubjects", "numberOfMilestone",
               "analysisTime", "efficacyRmst", "futilityRmst",
               "efficacyP", "futilityP", "information")]

    # format number of digits after decimal for each column
    j1 <- c(7,8,9,10,11)
    j2 <- 16
    j3 <- c(1,2,3,12,13)
    j4 <- c(4,5,6,14,15)

    b[j1] <- lapply(b[j1], formatC, format = "f", digits = 1)
    b[j2] <- lapply(b[j2], formatC, format = "f", digits = 2)
    b[j3] <- lapply(b[j3], formatC, format = "f", digits = 3)
    b[j4] <- lapply(b[j4], formatC, format = "f", digits = 4)

    if (x$settings$typeBetaSpending != 'none' ||
        (k > 1 && any(x$byStageResults$futilityBounds[1:(k-1)] > -8))) {
      df <- t(b)
      rownames(df) <- c("Information rate",
                        "Efficacy boundary (Z)",
                        "Futility boundary (Z)",
                        "Cumulative rejection",
                        "Cumulative futility",
                        "Cumulative alpha spent",
                        "Number of events",
                        "Number of dropouts",
                        "Number of subjects",
                        "Number of milestone subjects",
                        "Analysis time",
                        "Efficacy boundary (rmst)",
                        "Futility boundary (rmst)",
                        "Efficacy boundary (p)",
                        "Futility boundary (p)",
                        "Information")

    } else {
      df <- t(b[, c(1,2,4,6,7,8,9,10,11,12,14,16)])
      rownames(df) <- c("Information rate",
                        "Efficacy boundary (Z)",
                        "Cumulative rejection",
                        "Cumulative alpha spent",
                        "Number of events",
                        "Number of dropouts",
                        "Number of subjects",
                        "Number of milestone subjects",
                        "Analysis time",
                        "Efficacy boundary (rmst)",
                        "Efficacy boundary (p)",
                        "Information")
    }

    colnames(df) <- paste("Stage", seq_len(ncol(df)), sep=" ")
  } else {
    b <- s[, c("efficacyBounds", "efficacyRmst", "efficacyP")]

    # format number of digits after decimal for each column
    j3 <- c(1,2)
    j4 <- 3

    b[j3] <- lapply(b[j3], formatC, format = "f", digits = 3)
    b[j4] <- lapply(b[j4], formatC, format = "f", digits = 4)

    df <- t(b)

    rownames(df) <- c("Efficacy boundary (Z)",
                      "Efficacy boundary (rmst)",
                      "Efficacy boundary (p)")
    colnames(df) <- NA
  }

  print(df1, ..., na.print = "" , quote = FALSE )
  print(df, ..., na.print = "" , quote = FALSE )
  invisible(x)
}


#' @title Print liferegr Object
#' @description Prints the concise information of liferegr fit.
#'
#' @param x The liferegr object to print.
#' @param ... Ensures that all arguments starting from "..." are named.
#'
#' @return A printout from the fit of an accelerated failue time model.
#'
#' @keywords internal
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @export
print.liferegr <- function(x, ...) {
  lrchisq <- -2*(x$sumstat$loglik0 - x$sumstat$loglik1)
  degrees <- x$sumstat$nvar
  pvalue <- sapply(1:nrow(x$sumstat), function(i) {
    ifelse(degrees[i] > 0,
           pchisq(lrchisq[i], degrees[i], 0, lower.tail = FALSE),
           NA)
  })
  df1 <- cbind(x$sumstat[, c("n", "nevents", "loglik0", "loglik1")],
               lrchisq = lrchisq, df = degrees, pvalue = pvalue,
               x$sumstat[, c("niter", "dist")])

  p <- x$p
  if (p > 0) {
    if (!x$settings$robust) {
      if (x$settings$plci) {
        df <- data.frame(param = x$param,
                         coef = x$parest$beta,
                         expcoef = x$parest$expbeta,
                         se = x$parest$sebeta,
                         z = x$parest$z,
                         lower = x$parest$lower,
                         upper = x$parest$upper,
                         p = x$parest$p,
                         method = x$parest$method)

        colnames(df) <- c("param", "coef", "exp(coef)", "se(coef)", "z",
                          paste("lower", 1-x$settings$alpha),
                          paste("upper", 1-x$settings$alpha),
                          "p", "method")
      } else {
        df <- data.frame(param = x$param,
                         coef = x$parest$beta,
                         expcoef = x$parest$expbeta,
                         se = x$parest$sebeta,
                         z = x$parest$z,
                         p = x$parest$p)

        colnames(df) <- c("param", "coef", "exp(coef)", "se(coef)", "z", "p")
      }
    } else {
      if (x$settings$plci) {
        df <- data.frame(param = x$param,
                         coef = x$parest$beta,
                         expcoef = x$parest$expbeta,
                         nse = x$parest$sebeta_naive,
                         se = x$parest$sebeta,
                         z = x$parest$z,
                         lower = x$parest$lower,
                         upper = x$parest$upper,
                         p = x$parest$p,
                         method = x$parest$method)

        colnames(df) <- c("param", "coef", "exp(coef)", "se(coef)",
                          "robust se", "z",
                          paste("lower", 1-x$settings$alpha),
                          paste("upper", 1-x$settings$alpha),
                          "p", "method")
      } else {
        df <- data.frame(param = x$param,
                         coef = x$parest$beta,
                         expcoef = x$parest$expbeta,
                         nse = x$parest$sebeta_naive,
                         se = x$parest$sebeta,
                         z = x$parest$z,
                         p = x$parest$p)

        colnames(df) <- c("param", "coef", "exp(coef)", "se(coef)",
                          "robust se", "z", "p")
      }
    }
  }

  print(df1, ..., na.print = "" , quote = FALSE )
  cat("\n")
  print(df, ..., na.print = "" , quote = FALSE )
  invisible(x)
}


#' @title Print phregr Object
#' @description Prints the concise information of phregr fit.
#'
#' @param x The phregr object to print.
#' @param ... Ensures that all arguments starting from "..." are named.
#'
#' @return A printout from the fit of a Cox proportional hazards model.
#'
#' @keywords internal
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @export
print.phregr <- function(x, ...) {
  lrchisq <- -2*(x$sumstat$loglik0 - x$sumstat$loglik1)
  degrees <- x$sumstat$p
  pvalue <- sapply(1:nrow(x$sumstat), function(i) {
    ifelse(degrees[i] > 0,
           pchisq(lrchisq[i], degrees[i], 0, lower.tail = FALSE),
           NA)
  })
  df1 <- cbind(x$sumstat[, c("n", "nevents", "loglik0", "loglik1")],
               lrchisq = lrchisq, df = degrees, pvalue = pvalue,
               x$sumstat[, c("scoretest", "niter", "ties")])
  print(df1, ..., na.print = "" , quote = FALSE )
  cat("\n")

  p <- x$p
  if (p > 0) {
    if (!x$settings$robust) {
      if (x$settings$plci) {
        df <- data.frame(param = x$param,
                         coef = x$parest$beta,
                         expcoef = x$parest$expbeta,
                         se = x$parest$sebeta,
                         z = x$parest$z,
                         lower = x$parest$lower,
                         upper = x$parest$upper,
                         p = x$parest$p,
                         method = x$parest$method)

        colnames(df) <- c("param", "coef", "exp(coef)", "se(coef)", "z",
                          paste("lower", 1-x$settings$alpha),
                          paste("upper", 1-x$settings$alpha), "p", "method")

      } else {
        df <- data.frame(param = x$param,
                         coef = x$parest$beta,
                         expcoef = x$parest$expbeta,
                         se = x$parest$sebeta,
                         z = x$parest$z,
                         p = x$parest$p)

        colnames(df) <- c("param", "coef", "exp(coef)", "se(coef)", "z", "p")
      }
    } else {
      if (x$settings$plci) {
        df <- data.frame(param = x$param,
                         coef = x$parest$beta,
                         expcoef = x$parest$expbeta,
                         nse = x$parest$sebeta_naive,
                         se = x$parest$sebeta,
                         z = x$parest$z,
                         lower = x$parest$lower,
                         upper = x$parest$upper,
                         p = x$parest$p,
                         method = x$parest$method)

        colnames(df) <- c("param", "coef", "exp(coef)", "se(coef)",
                          "robust se", "z",
                          paste("lower", 1-x$settings$alpha),
                          paste("upper", 1-x$settings$alpha),
                          "p", "method")
      } else {
        df <- data.frame(param = x$param,
                         coef = x$parest$beta,
                         expcoef = x$parest$expbeta,
                         nse = x$parest$sebeta_naive,
                         se = x$parest$sebeta,
                         z = x$parest$z,
                         p = x$parest$p)

        colnames(df) <- c("param", "coef", "exp(coef)", "se(coef)",
                          "robust se", "z", "p")
      }
    }

    print(df, ..., na.print = "" , quote = FALSE )
  }

  invisible(x)
}



#' @title Print logisregr Object
#' @description Prints the concise information of logisregr fit.
#'
#' @param x The logisregr object to print.
#' @param ... Ensures that all arguments starting from "..." are named.
#'
#' @return A printout from the fit of a logistic regression model.
#'
#' @keywords internal
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @export
print.logisregr <- function(x, ...) {
  lrchisq <- -2*(x$sumstat$loglik0 - x$sumstat$loglik1)
  degrees <- x$sumstat$p - 1
  pvalue <- sapply(1:nrow(x$sumstat), function(i) {
    ifelse(degrees[i] > 0,
           pchisq(lrchisq[i], degrees[i], 0, lower.tail = FALSE),
           NA)
  })
  df1 <- cbind(x$sumstat[, c("n", "nevents", "loglik0", "loglik1")],
               lrchisq = lrchisq, df = degrees, pvalue = pvalue,
               x$sumstat[, c("niter", "link", "firth", "flic")])

  p <- x$p
  if (p > 0) {
    if (!x$settings$robust) {
      if (x$settings$plci) {
        df <- data.frame(param = x$param,
                         coef = x$parest$beta,
                         expcoef = x$parest$expbeta,
                         se = x$parest$sebeta,
                         z = x$parest$z,
                         lower = x$parest$lower,
                         upper = x$parest$upper,
                         p = x$parest$p,
                         method = x$parest$method)

        colnames(df) <- c("param", "coef", "exp(coef)", "se(coef)", "z",
                          paste("lower", 1-x$settings$alpha),
                          paste("upper", 1-x$settings$alpha),
                          "p", "method")

      } else {
        df <- data.frame(param = x$param,
                         coef = x$parest$beta,
                         expcoef = x$parest$expbeta,
                         se = x$parest$sebeta,
                         z = x$parest$z,
                         p = x$parest$p)

        colnames(df) <- c("param", "coef", "exp(coef)", "se(coef)", "z", "p")
      }
    } else {
      if (x$settings$plci) {
        df <- data.frame(param = x$param,
                         coef = x$parest$beta,
                         expcoef = x$parest$expbeta,
                         nse = x$parest$sebeta_naive,
                         se = x$parest$sebeta,
                         z = x$parest$z,
                         lower = x$parest$lower,
                         upper = x$parest$upper,
                         p = x$parest$p,
                         method = x$parest$method)

        colnames(df) <- c("param", "coef", "exp(coef)", "se(coef)",
                          "robust se", "z",
                          paste("lower", 1-x$settings$alpha),
                          paste("upper", 1-x$settings$alpha),
                          "p", "method")
      } else {
        df <- data.frame(param = x$param,
                         coef = x$parest$beta,
                         expcoef = x$parest$expbeta,
                         nse = x$parest$sebeta_naive,
                         se = x$parest$sebeta,
                         z = x$parest$z,
                         p = x$parest$p)

        colnames(df) <- c("param", "coef", "exp(coef)", "se(coef)",
                          "robust se", "z", "p")
      }
    }
  }

  print(df1, ..., na.print = "" , quote = FALSE )
  cat("\n")
  print(df, ..., na.print = "" , quote = FALSE )
  invisible(x)
}


#' @title Print method for assess_phregr objects
#' @description Prints the concise information of an assess_phregr fit.
#'
#' @param x An object of class \code{assess_phregr}.
#' @param ... Ensures that all arguments starting from "..." are named.
#'
#' @return A printout from the fit of an assessment of proportional hazards
#' assumption of a Cox model.
#'
#' @keywords internal
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @export
print.assess_phregr <- function(x, ...) {

  format_pvalue <- function(p) {
    # Handle the case of p < 0.0001
    ifelse(p < 0.0001,
           "<.0001",
           # Handle the case of p > 0.9999
           ifelse(p > 0.9999,
                  ">.9999",
                  # For all other cases, format to 4 decimal places
                  sprintf("%.4f", p)))
  }

  df <- data.frame(covariate = c(x$covariates, "GLOBAL"),
                   max_abs_value = x$max_abs_value,
                   resample = x$resample,
                   seed = x$seed,
                   p_value = format_pvalue(x$p_value))

  j0 <- 2
  df[j0] <- lapply(df[j0], formatC, format = "f", digits = 4)
  print(df, ..., na.print = "", quote = FALSE)

  invisible(x)
}


#' @title Print Phase 2/3 Seamless Design
#' @description Prints the stopping boundaries and power for
#' phase 2/3 seamless design.
#'
#' @param x The design object to print.
#' @param ... Ensures that all arguments starting from "..." are named.
#'
#' @return A tabular printout of the design elements.
#'
#' @keywords internal
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @export
print.seamless <- function(x, ...) {
  a <- x$overallResults
  s <- x$byStageResults
  t <- x$byArmResults
  k <- a$K

  # overall results
  if (k>1) {
    str1 <- "Phase 2/3 seamless group-sequential design"
  } else {
    str1 <- "Phase 2/3 seamless design"
  }

  str2 <- paste0("Overall power: ",
                 round(a$overallReject, 4), ", ",
                 "overall alpha (1-sided): ",
                 round(a$alpha, 4))

  if (x$settings$typeBetaSpending != 'none' ||
      (any(x$byStageResults$futilityBounds[1:k] > -8))) {
    str3 <- paste0(str2, ", ",
                   "attained alpha: ", round(a$attainedAlpha, 4))
  }

  str3 <- paste0("Number of active arms in phase 2: ", a$M)

  str4 <- paste0("Randomization ratio of each active vs. control: ", a$r)

  str5 <- paste0("Using correlation for critical value calculation: ",
                 a$corr_known)

  str6 <- paste0("Number of looks in phase 3: ", a$K)

  str7 <- paste0("Max information for pairwise comparion: ",
                 round(a$information, 2))

  str8 <- paste0("Expected information under H1: ",
                 round(a$expectedInformationH1, 2), ", ",
                 "expected information under H0: ",
                 round(a$expectedInformationH0, 2))

  str9 <- paste0("Max information for oveall study: ",
                 round(a$informationOverall, 2))

  str10 <- paste0("Expected overall info under H1: ",
                  round(a$expectedInformationOverallH1, 2), ", ",
                  "expected overall info under H0: ",
                  round(a$expectedInformationOverallH0, 2))

  asf <- tolower(x$settings$typeAlphaSpending)
  asfpar <- round(x$settings$parameterAlphaSpending, 3)
  asfuser <- round(x$settings$userAlphaSpending, 4)

  bsf <- tolower(x$settings$typeBetaSpending)
  bsfpar <- round(x$settings$parameterBetaSpending, 3)
  bsfuser <- round(x$settings$userBetaSpending, 4)

  if (asf == "of") {
    str11 <- paste0("Alpha spending: O'Brien-Fleming")
  } else if (asf == "p") {
    str11 <- paste0("Alpha spending: Pocock")
  } else if (asf == "wt") {
    str11 <- paste0("Alpha spending: Wang-Tsiatis(Delta = ", asfpar, ")")
  } else if (asf == "sfof") {
    str11 <- paste0("Alpha spending: Lan-DeMets O'Brien-Fleming")
  } else if (asf == "sfp") {
    str11 <- paste0("Alpha spending: Lan-DeMets Pocock")
  } else if (asf == "sfkd") {
    str11 <- paste0("Alpha spending: KD(rho = ", asfpar, ")")
  } else if (asf == "sfhsd") {
    str11 <- paste0("Alpha spending: HSD(gamma = ", asfpar, ")")
  } else if (asf == "user") {
    str11 <- paste0("Alpha spending: User defined(",
                   paste(asfuser, collapse = ","), ")")
  } else {
    str11 <- "Alpha spending: None"
  }

  if (bsf == "of") {
    str12 <- paste0("beta spending: O'Brien-Fleming")
  } else if (bsf == "p") {
    str12 <- paste0("beta spending: Pocock")
  } else if (bsf == "wt") {
    str12 <- paste0("beta spending: Wang-Tsiatis(Delta = ", bsfpar, ")")
  } else if (bsf == "sfof") {
    str12 <- paste0("beta spending: Lan-DeMets O'Brien-Fleming")
  } else if (bsf == "sfp") {
    str12 <- paste0("beta spending: Lan-DeMets Pocock")
  } else if (bsf == "sfkd") {
    str12 <- paste0("beta spending: KD(rho = ", bsfpar, ")")
  } else if (bsf == "sfhsd") {
    str12 <- paste0("beta spending: HSD(gamma = ", bsfpar, ")")
  } else if (bsf == "user") {
    str12 <- paste0("beta spending: User defined(",
                    paste(bsfuser, collapse = ","), ")")
  } else {
    str12 <- "beta spending: None"
  }

  if (!any(is.na(x$settings$spendingTime)) &&
      !all(x$settings$spendingTime == s$informationRates)) {
    str13 <- paste0("Spending time: ",
                   paste(round(x$settings$spendingTime, 3), collapse = ","))
    df1 <- data.frame(x = rep("", 14))
    colnames(df1) <- NULL
    rownames(df1) <- c(str1, str2, str3, str4, str5, str6, str7,
                       str8, str9, str10,
                       paste(str11, str12, sep = ", "), str13, "")
  } else {
    df1 <- data.frame(x = rep("", 12))
    colnames(df1) <- NULL
    rownames(df1) <- c(str1, str2, str3, str4, str5, str6, str7,
                       str8, str9, str10,
                       paste(str11, str12, sep = ", "), "")
  }

  # by stage results
  b <- s[, c("informationRates",
             "efficacyBounds",
             "futilityBounds",
             "cumulativeRejection",
             "cumulativeFutility",
             "cumulativeAlphaSpent",
             "efficacyTheta",
             "futilityTheta",
             "efficacyP",
             "futilityP",
             "information",
             "informationOverall",
             "cumulativeRejectionH0",
             "cumulativeFutilityH0")]

  # format number of digits after decimal for each column
  j2 <- c(11,12)
  j3 <- c(1,2,3,7,8)
  j4 <- c(4,5,6,9,10,13,14)

  b[j2] <- lapply(b[j2], formatC, format = "f", digits = 2)
  b[j3] <- lapply(b[j3], formatC, format = "f", digits = 3)
  b[j4] <- lapply(b[j4], formatC, format = "f", digits = 4)

  if (x$settings$typeBetaSpending != 'none' ||
      (any(x$byStageResults$futilityBounds[1:k] > -8))) {

    df2 <- t(b)
    rownames(df2) <- c("Information rate",
                       "Efficacy boundary (Z)",
                       "Futility boundary (Z)",
                       "Cumulative rejection",
                       "Cumulative futility",
                       "Cumulative alpha spent",
                       "Efficacy boundary (theta)",
                       "Futility boundary (theta)",
                       "Efficacy boundary (p)",
                       "Futility boundary (p)",
                       "Information for pairwise comp",
                       "Information for overall study",
                       "Cumulative rejection under H0",
                       "Cumulative futility under H0")

    colnames(df2) <- paste("Stage", seq_len(ncol(df2)), sep=" ")
  } else {
    df2 <- t(b[, c(1,2,4,6,7,9,11,12)])
    rownames(df2) <- c("Information rate",
                       "Efficacy boundary (Z)",
                       "Cumulative rejection",
                       "Cumulative alpha spent",
                       "Efficacy boundary (theta)",
                       "Efficacy boundary (p)",
                       "Information for pairwise comp",
                       "Information for overall study")

    colnames(df2) <- paste("Stage", seq_len(ncol(df2)), sep=" ")
  }

  # by arm results
  j3 <- 1
  j4 <- c(2,3,4)
  t[j3] <- lapply(t[j3], formatC, format = "f", digits = 3)
  t[j4] <- lapply(t[j4], formatC, format = "f", digits = 4)

  df3 <- t(t)
  rownames(df3) <- c("Treatment effect (theta)",
                     "Being the best in phase 2",
                     "Power",
                     "Conditional power")
  colnames(df3) <- paste("Arm", seq_len(ncol(df3)), sep=" ")

  print(df1, ..., na.print = "" , quote = FALSE )
  print(df2, ..., na.print = "" , quote = FALSE )
  cat("\n")
  print(df3, ..., na.print = "" , quote = FALSE )
  invisible(x)
}


#' @title Print Adaptive Phase 2/3 Seamless Design
#' @description Prints the stopping boundaries and power for adaptive
#' phase 2/3 seamless design.
#'
#' @param x The design object to print.
#' @param ... Ensures that all arguments starting from "..." are named.
#'
#' @return A tabular printout of the design elements.
#'
#' @keywords internal
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @export
print.adaptDesign_seamless <- function(x, ...) {
  des1 <- x$primaryTrial

  str1 <- "Primary trial:"

  if (des1$K>1) {
    str2 <- "Phase 2/3 seamless group-sequential design"
  } else {
    str2 <- "Phase 2/3 seamless design"
  }

  str3 <- paste0("Number of active arms in phase 2: ", des1$M)

  str4 <- paste0("Randomization ratio of each active vs. control: ", des1$r)

  str5 <- paste0("Using correlation for critical value calculation: ",
                 des1$corr_known)

  str6 <- paste0("Number of looks in phase 3: ", des1$K)

  str7 <- paste0("Max information for pairwise comparion: ",
                 round(des1$maxInformation, 2))

  str8 <- paste0("Interim adaptation look in Phase 3: ", des1$L, ", ",
                 "z-statistic value: ", paste(round(des1$zL, 3), collapse = ", "))

  str9 <- paste0("theta: ", round(des1$theta, 3))

  str10 <- paste0("Conditional type I error: ", round(des1$conditionalAlpha, 4),
                  ", conditional power: ", round(des1$conditionalPower, 3))

  str11 <- paste0("Muller & Schafer method for secondary trial: ",
                  des1$MullerSchafer)

  df1a <- data.frame(x = rep("", 12))
  colnames(df1a) <- NULL
  rownames(df1a) <- c(str1, str2, str3, str4, str5, str6, str7, str8,
                      str9, str10, str11, "")

  b <- data.frame(informationRates = des1$informationRates,
                  efficacyBounds = des1$efficacyBounds,
                  futilityBounds = des1$futilityBounds,
                  information = des1$information)

  b[1:3] <- lapply(b[1:3], formatC, format = "f", digits = 3)
  b[4] <- lapply(b[4], formatC, format = "f", digits = 2)

  if (des1$kMax > 1 && any(des1$futilityBounds[1:(des1$kMax-1)] > -8)) {
    df1b <- t(b)
    rownames(df1b) <- c("Information rate",
                        "Efficacy boundary (Z)",
                        "Futility boundary (Z)",
                        "Information")
  } else {
    df1b <- t(b[, c(1,2,4)])
    rownames(df1b) <- c("Information rate",
                        "Efficacy boundary (Z)",
                        "Information")
  }
  colnames(df1b) <- paste("Stage", seq_len(ncol(df1b)), sep=" ")

  des2 <- x$secondaryTrial
  k <- des2$kMax

  str1 <- "Secondary trial:"

  if (k>1) {
    str2 <- paste0("Group-sequential design with ", k, " stages")
  } else {
    str2 <- "Fixed design"
  }

  str3 <- paste0("Maximum information: ", round(des2$maxInformation, 2))

  str4 <- paste0("Overall power: ",
                 round(des2$overallReject, 4), ", ",
                 "overall significance level (1-sided): ",
                 round(des2$alpha, 4))

  df2a <- data.frame(x = rep("", 5))
  colnames(df2a) <- NULL
  rownames(df2a) <- c(str1, str2, str3, str4, "")


  if (k>1) {
    b <- data.frame(informationRates = des2$informationRates,
                    efficacyBounds = des2$efficacyBounds,
                    futilityBounds = des2$futilityBounds,
                    cumulativeRejection = des2$cumulativeRejection,
                    cumulativeFutility = des2$cumulativeFutility,
                    cumulativeAlphaSpent = des2$cumulativeAlphaSpent,
                    information = des2$information)

    # format number of digits after decimal for each column
    j2 <- 7
    j3 <- c(1,2,3)
    j4 <- c(4,5,6)

    b[j2] <- lapply(b[j2], formatC, format = "f", digits = 2)
    b[j3] <- lapply(b[j3], formatC, format = "f", digits = 3)
    b[j4] <- lapply(b[j4], formatC, format = "f", digits = 4)

    if (des2$typeBetaSpending != 'none' ||
        (k > 1 && any(des2$futilityBounds[1:(k-1)] > -8))) {
      df2b <- t(b)
      rownames(df2b) <- c("Information rate",
                          "Efficacy boundary (Z)",
                          "Futility boundary (Z)",
                          "Cumulative rejection",
                          "Cumulative futility",
                          "Cumulative alpha spent",
                          "Information")

    } else {
      df2b <- t(b[, c(1,2,4,6,7)])
      rownames(df2b) <- c("Information rate",
                          "Efficacy boundary (Z)",
                          "Cumulative rejection",
                          "Cumulative alpha spent",
                          "Information")
    }

    colnames(df2b) <- paste("Stage", seq_len(ncol(df2b)), sep=" ")
  } else {
    b <- data.frame(efficacyBounds = des2$efficacyBounds)

    # format number of digits after decimal for each column
    j3 <- 1

    b[j3] <- lapply(b[j3], formatC, format = "f", digits = 3)

    df2b <- t(b)

    rownames(df2b) <- "Efficacy boundary (Z)"

    colnames(df2b) <- NA
  }



  des3 <- x$integratedTrial

  str1 <- "Integrated trial:"

  str2 <- "Adaptive Phase 2/3 seamless design"

  str3 <- paste0("Total number of looks in Phase 3: ", des3$kMax - 1)
  str4 <- paste0("Maximum information for pairwise comparion: ",
                 round(des3$maxInformation, 2))

  str5 <- paste0("Interim adaptation look in Phase 3: ", des3$L, ", ",
                 "z-statistic value: ", paste(round(des3$zL, 3), collapse = ", "))

  df3a <- data.frame(x = rep("", 6))
  colnames(df3a) <- NULL
  rownames(df3a) <- c(str1, str2, str3, str4, str5, "")

  b <- data.frame(informationRates = des3$informationRates,
                  efficacyBounds = des3$efficacyBounds,
                  futilityBounds = des3$futilityBounds,
                  information = des3$information)

  # format number of digits after decimal for each column
  j2 <- 4
  j3 <- c(1,2,3)

  b[j2] <- lapply(b[j2], formatC, format = "f", digits = 2)
  b[j3] <- lapply(b[j3], formatC, format = "f", digits = 3)


  if ((des3$kMax > 1 && any(des3$futilityBounds[1:(des3$kMax-1)] > -8))) {
    df3b <- t(b)
    rownames(df3b) <- c("Information rate",
                        "Efficacy bounds (Z)",
                        "Futility bounds (Z)",
                        "Information")
  } else {
    df3b <- t(b[, c(1,2,4)])
    rownames(df3b) <- c("Information rate",
                        "Efficacy bounds (Z)",
                        "Information")
  }
  colnames(df3b) <- paste("Stage", seq_len(ncol(df3b)), sep=" ")


  print(df1a, ..., na.print = "" , quote = FALSE )
  print(df1b, ..., na.print = "" , quote = FALSE )
  print(df2a, ..., na.print = "" , quote = FALSE )
  print(df2b, ..., na.print = "" , quote = FALSE )
  print(df3a, ..., na.print = "" , quote = FALSE )
  print(df3b, ..., na.print = "" , quote = FALSE )
  invisible(x)
}


#' @title Print Simulation Results for Phase 2/3 Seamless Design
#' @description Prints the summary statistics from simulation.
#'
#' @param x The lrsim_seamless object to print.
#' @param ... Ensures that all arguments starting from "..." are named.
#'
#' @return A tabular printout of the summary statistics from simulation runs.
#'
#' @keywords internal
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @export
print.lrsim_seamless <- function(x, ...) {
  a <- x$overview
  k <- a$K + 1

  if (k>1) {
    str1 <- paste0("Phase 2/3 seamless group-sequential design")
  } else {
    str1 <- "Phase 2/3 seamless design"
  }

  if (a$rho1 != 0 || a$rho2 != 0) {
    str1 <- paste0(str1, " for weighted log-rank test, FH(",
                   a$rho1, ", ", a$rho2, ")")
  } else {
    str1 <- paste0(str1, " for log-rank test")
  }

  str2 <- paste0("Empirical power: ", round(a$overallReject, 4))

  str3 <- paste0("Number of active arms in phase 2: ", a$M)

  str4 <- paste0("Number of looks in phase 3: ", a$K)

  str5 <- paste0("Expected # events: ",
                 round(a$expectedNumberOfEvents, 1))

  str6 <- paste0("Expected # dropouts: ",
                 round(a$expectedNumberOfDropouts, 1))

  str7 <- paste0("Expected # subjects: ",
                 round(a$expectedNumberOfSubjects, 1))

  str8 <- paste0("Expected study duration: ",
                 round(a$expectedStudyDuration, 1))

  str9 <- paste0("n: ", a$n, ", ",
                 "fixed follow-up: ", a$fixedFollowup)

  str10 <- paste0("Number of simulations: ", a$numberOfIterations)

  str11 <- paste0("Efficacy bounds (z-scale): ",
                  paste(round(a$criticalValues, 3), collapse = ", "))

  df1 <- data.frame(x = rep("", 12))
  colnames(df1) <- NULL
  rownames(df1) <- c(str1, str2, str3, str4, str5, str6, str7, str8,
                     str9, str10, str11, "")

  df2 <- t(data.frame(selectAsBest = a$selectAsBest))
  rownames(df2) <- "Selected as best in phase 2"
  colnames(df2) <- paste("Arm", seq_len(ncol(df2)), sep=" ")

  df3 <- data.frame(activeArm = rep(c(seq(1, a$M), "Overall"), each = k),
                    stage = rep(seq(1, k), times = a$M + 1),
                    cumReject = c(a$cumulativeRejection),
                    nEvents = c(a$numberOfEvents),
                    nDropouts = c(a$numberOfDropouts),
                    nSubjects = c(a$numberOfSubjects),
                    analysisTime = c(a$analysisTime))

  # format number of digits after decimal for each column
  j1 <- c(4,5,6,7)
  j4 <- 3
  df3[j1] <- lapply(df3[j1], formatC, format = "f", digits = 1)
  df3[j4] <- lapply(df3[j4], formatC, format = "f", digits = 4)

  print(df1, ..., na.print = "" , quote = FALSE )
  print(df2, ..., na.print = "" , quote = FALSE )
  cat("\n")
  print(df3, ..., na.print = "" , quote = FALSE, row.names= FALSE)
  invisible(x)
}


#' @title Print Multi-Arm Multi-Stage Design
#' @description Prints the stopping boundaries and power for
#' multi-arm multi-stage design.
#'
#' @param x The design object to print.
#' @param ... Ensures that all arguments starting from "..." are named.
#'
#' @return A tabular printout of the design elements.
#'
#' @keywords internal
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @export
print.mams <- function(x, ...) {
  a <- x$overallResults
  s <- x$byStageResults
  t <- x$settings
  k <- a$kMax

  # overall results
  str1 <- "Multi-arm multi-stage design"

  str2 <- paste0("Overall power: ",
                 round(a$overallReject, 4), ", ",
                 "overall alpha (1-sided): ",
                 round(a$alpha, 4))

  if (x$settings$typeBetaSpending != 'none' ||
      (k > 1 && any(x$byStageResults$futilityBounds[1:(k-1)] > -8))) {
    str2 <- paste0(str2, ", ",
                   "attained alpha: ", round(a$attainedAlpha, 4))
  }

  str3 <- paste0("Number of active arms: ", a$M)

  str4 <- paste0("Randomization ratio of each active vs. control: ", a$r)

  str5 <- paste0("Using correlation for critical value calculation: ",
                 a$corr_known)

  str6 <- paste0("Number of looks: ", a$kMax)

  str7 <- paste0("Max information for pairwise comparion: ",
                 round(a$information, 2))

  str8 <- paste0("Max information for overall study: ",
                 round(a$informationOverall, 2))


  if (k > 1) {
    str9 <- paste0("Expected information under H1: ",
                   round(a$expectedInformationH1, 2), ", ",
                   "expected information under H0: ",
                   round(a$expectedInformationH0, 2))

    str10 <- paste0("Expected overall info under H1: ",
                    round(a$expectedInformationOverallH1, 2), ", ",
                    "expected overall info under H0: ",
                    round(a$expectedInformationOverallH0, 2))

    asf <- tolower(x$settings$typeAlphaSpending)
    asfpar <- round(x$settings$parameterAlphaSpending, 3)
    asfuser <- round(x$settings$userAlphaSpending, 4)

    bsf <- tolower(x$settings$typeBetaSpending)
    bsfpar <- round(x$settings$parameterBetaSpending, 3)
    bsfuser <- round(x$settings$userBetaSpending, 4)

    if (asf == "of") {
      str11 <- paste0("Alpha spending: O'Brien-Fleming")
    } else if (asf == "p") {
      str11 <- paste0("Alpha spending: Pocock")
    } else if (asf == "wt") {
      str11 <- paste0("Alpha spending: Wang-Tsiatis(Delta = ", asfpar, ")")
    } else if (asf == "sfof") {
      str11 <- paste0("Alpha spending: Lan-DeMets O'Brien-Fleming")
    } else if (asf == "sfp") {
      str11 <- paste0("Alpha spending: Lan-DeMets Pocock")
    } else if (asf == "sfkd") {
      str11 <- paste0("Alpha spending: KD(rho = ", asfpar, ")")
    } else if (asf == "sfhsd") {
      str11 <- paste0("Alpha spending: HSD(gamma = ", asfpar, ")")
    } else if (asf == "user") {
      str11 <- paste0("Alpha spending: User defined(",
                     paste(asfuser, collapse = ","), ")")
    } else {
      str11 <- "Alpha spending: None"
    }

    if (bsf == "of") {
      str12 <- paste0("beta spending: O'Brien-Fleming")
    } else if (bsf == "p") {
      str12 <- paste0("beta spending: Pocock")
    } else if (bsf == "wt") {
      str12 <- paste0("beta spending: Wang-Tsiatis(Delta = ", bsfpar, ")")
    } else if (bsf == "sfof") {
      str12 <- paste0("beta spending: Lan-DeMets O'Brien-Fleming")
    } else if (bsf == "sfp") {
      str12 <- paste0("beta spending: Lan-DeMets Pocock")
    } else if (bsf == "sfkd") {
      str12 <- paste0("beta spending: KD(rho = ", bsfpar, ")")
    } else if (bsf == "sfhsd") {
      str12 <- paste0("beta spending: HSD(gamma = ", bsfpar, ")")
    } else if (bsf == "user") {
      str12 <- paste0("beta spending: User defined(",
                      paste(bsfuser, collapse = ","), ")")
    } else {
      str12 <- "beta spending: None"
    }

    if (!any(is.na(x$settings$spendingTime)) &&
        !all(x$settings$spendingTime == s$informationRates)) {
      str13 <- paste0("Spending time: ",
                      paste(round(x$settings$spendingTime, 3), collapse = ","))
      df1 <- data.frame(x = rep("", 14))
      colnames(df1) <- NULL
      rownames(df1) <- c(str1, str2, str3, str4, str5, str6, str7,
                         str8, str9, str10,
                         paste(str11, str12, sep = ", "), str13, "")
    } else {
      df1 <- data.frame(x = rep("", 12))
      colnames(df1) <- NULL
      rownames(df1) <- c(str1, str2, str3, str4, str5, str6, str7,
                         str8, str9, str10,
                         paste(str11, str12, sep = ", "), "")
    }
  } else {
    df1 <- data.frame(x = rep("", 9))
    colnames(df1) <- NULL
    rownames(df1) <- c(str1, str2, str3, str4, str5, str6, str7, str8, "")
  }


  # by stage results
  if (k > 1) {
    b <- s[, c("informationRates",
               "efficacyBounds",
               "futilityBounds",
               "cumulativeRejection",
               "cumulativeFutility",
               "cumulativeAlphaSpent",
               "efficacyTheta",
               "futilityTheta",
               "efficacyP",
               "futilityP",
               "information",
               "informationOverall",
               "cumulativeRejectionH0",
               "cumulativeFutilityH0")]

    # format number of digits after decimal for each column
    j2 <- c(11,12)
    j3 <- c(1,2,3,7,8)
    j4 <- c(4,5,6,9,10,13,14)

    b[j2] <- lapply(b[j2], formatC, format = "f", digits = 2)
    b[j3] <- lapply(b[j3], formatC, format = "f", digits = 3)
    b[j4] <- lapply(b[j4], formatC, format = "f", digits = 4)

    if (x$settings$typeBetaSpending != 'none' ||
        (k > 1 && any(x$byStageResults$futilityBounds[1:(k-1)] > -8))) {
      df2 <- t(b)
      rownames(df2) <- c("Information rate",
                         "Efficacy boundary (Z)",
                         "Futility boundary (Z)",
                         "Cumulative rejection",
                         "Cumulative futility",
                         "Cumulative alpha spent",
                         "Efficacy boundary (theta)",
                         "Futility boundary (theta)",
                         "Efficacy boundary (p)",
                         "Futility boundary (p)",
                         "Information for pairwise comp",
                         "Information for overall study",
                         "Cumulative rejection under H0",
                         "Cumulative futility under H0")

    } else {
      df2 <- t(b[, c(1,2,4,6,7,9,11,12)])
      rownames(df2) <- c("Information rate",
                         "Efficacy boundary (Z)",
                         "Cumulative rejection",
                         "Cumulative alpha spent",
                         "Efficacy boundary (theta)",
                         "Efficacy boundary (p)",
                         "Information for pairwise comp",
                         "Information for overall study")
    }

    colnames(df2) <- paste("Stage", seq_len(ncol(df2)), sep=" ")
  } else {
    b <- s[, c("efficacyBounds", "efficacyTheta", "efficacyP")]

    # format number of digits after decimal for each column
    j3 <- c(1,2)
    j4 <- 3

    b[j3] <- lapply(b[j3], formatC, format = "f", digits = 3)
    b[j4] <- lapply(b[j4], formatC, format = "f", digits = 4)

    df2 <- t(b)

    rownames(df2) <- c("Efficacy boundary (Z)",
                       "Efficacy boundary (theta)",
                       "Efficacy boundary (p)")

    colnames(df2) <- NA
  }

  df3 <- x$byLevelBounds
  df3[3] <- lapply(df3[3], formatC, format = "f", digits = 3)
  colnames(df3) <- c("Level", "Stage", "Boundary (Z)")

  print(df1, ..., na.print = "" , quote = FALSE )
  print(df2, ..., na.print = "" , quote = FALSE )
  if (a$M > 1) {
    cat("\nBy level critical values\n")
    print(df3, ..., na.print = "" , quote = FALSE)
  }
  invisible(x)
}



#' @title Print Adaptive Multi-Arm Multi-Stage Design
#' @description Prints the stopping boundaries and power for adaptive
#' multi-arm multi-stage design.
#'
#' @param x The design object to print.
#' @param ... Ensures that all arguments starting from "..." are named.
#'
#' @return A tabular printout of the design elements.
#'
#' @keywords internal
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @export
print.adaptDesign_mams <- function(x, ...) {
  des1 <- x$primaryTrial

  str1 <- "Primary trial:"
  str2 <- "Multi-arm multi-stage design"
  str3 <- paste0("Number of active arms: ", des1$M)
  str4 <- paste0("Randomization ratio of each active vs. control: ", des1$r)
  str5 <- paste0("Using correlation for critical value calculation: ",
                 des1$corr_known)
  str6 <- paste0("Max information for pairwise comparion: ",
                 round(des1$maxInformation, 2))
  str7 <- paste0("Number of looks: ", des1$kMax)
  str8 <- paste0("Interim adaptation look: ", des1$L, ", ",
                 "z-statistic value: ", paste(round(des1$zL, 3), collapse = ", "))
  str9 <- paste0("theta: ", paste(round(des1$theta, 3), collapse = ", "))
  str10 <- paste0("Conditional type I error: ", round(des1$conditionalAlpha, 4),
                  ", conditional power: ", round(des1$conditionalPower, 3))
  str11 <- paste0("Muller & Schafer method for secondary trial: ",
                  des1$MullerSchafer)

  df1a <- data.frame(x = rep("", 12))
  colnames(df1a) <- NULL
  rownames(df1a) <- c(str1, str2, str3, str4, str5, str6, str7, str8,
                      str9, str10, str11, "")

  b <- data.frame(informationRates = des1$informationRates,
                  efficacyBounds = des1$efficacyBounds,
                  futilityBounds = des1$futilityBounds,
                  information = des1$information)

  b[1:3] <- lapply(b[1:3], formatC, format = "f", digits = 3)
  b[4] <- lapply(b[4], formatC, format = "f", digits = 2)

  if (des1$kMax > 1 && any(des1$futilityBounds[1:(des1$kMax-1)] > -8)) {
    df1b <- t(b)
    rownames(df1b) <- c("Information rate",
                        "Efficacy boundary (Z)",
                        "Futility boundary (Z)",
                        "Information")
  } else {
    df1b <- t(b[, c(1,2,4)])
    rownames(df1b) <- c("Information rate",
                        "Efficacy boundary (Z)",
                        "Information")
  }
  colnames(df1b) <- paste("Stage", seq_len(ncol(df1b)), sep=" ")


  df1c <- des1$byLevelBounds
  df1c[3] <- lapply(df1c[3], formatC, format = "f", digits = 3)
  colnames(df1c) <- c("Level", "Stage", "Boundary (Z)")


  des2 <- x$secondaryTrial

  # overall results
  str1 <- "Secondary trial:"
  str2 <- "Multi-arm multi-stage design"
  str3 <- paste0("Number of selected active arms: ", des2$M, ", ",
                 "selected active arms: ", paste(des2$selected, collapse = ", "))
  str4 <- paste0("Randomization ratio of each active vs. control: ", des2$r)
  str5 <- paste0("Maximum information: ", round(des2$maxInformation, 2))
  str6 <- paste0("Overall power: ",
                 round(des2$overallReject, 4), ", ",
                 "overall significance level (1-sided): ",
                 round(des2$alpha, 4))

  df2a <- data.frame(x = rep("", 7))
  colnames(df2a) <- NULL
  rownames(df2a) <- c(str1, str2, str3, str4, str5, str6, "")

  b <- data.frame(informationRates = des2$informationRates,
                  cumulativeRejection = des2$cumulativeRejection,
                  cumulativeAlphaSpent = des2$cumulativeAlphaSpent,
                  information = des2$information)

  # format number of digits after decimal for each column
  j2 <- 4
  j3 <- 1
  j4 <- c(2,3)

  b[j2] <- lapply(b[j2], formatC, format = "f", digits = 2)
  b[j3] <- lapply(b[j3], formatC, format = "f", digits = 3)
  b[j4] <- lapply(b[j4], formatC, format = "f", digits = 4)

  df2b <- t(b)
  rownames(df2b) <- c("Information rate",
                      "Cumulative rejection",
                      "Cumulative alpha spent",
                      "Information")

  colnames(df2b) <- paste("Stage", seq_len(ncol(df2b)), sep=" ")

  df2c <- des2$byHypothesisBounds
  df2c[c(3,4)] <- lapply(df2c[c(3,4)], formatC, format = "f", digits = 3)
  colnames(df2c) <- c("Hypothesis",
                      "Stage",
                      "Efficacy boundary (Z)",
                      "Futility boundary (Z)")

  des3 <- x$integratedTrial

  str1 <- "Integrated trial:"
  str2 <- paste0("Adaptive multi-arm multi-stage design")
  str3 <- paste0("Number of active arms before adaptation: ", des3$M)
  str4 <- paste0("Number of selected active arms: ", des3$MNew, ", ",
                 "selected active arms: ", paste(des3$selected, collapse = ", "))
  str5 <- paste0("Total number of looks: ", des3$kMax)
  str6 <- paste0("Interim adaptation look: ", des3$L, ", ",
                 "z-statistic value: ", paste(round(des3$zL, 3), collapse = ", "))

  df3a <- data.frame(x = rep("", 7))
  colnames(df3a) <- NULL
  rownames(df3a) <- c(str1, str2, str3, str4, str5, str6, "")

  b <- data.frame(informationRates = des3$informationRates,
                  efficacyBounds = des3$efficacyBounds,
                  futilityBounds = des3$futilityBounds,
                  information = des3$information)

  # format number of digits after decimal for each column
  j2 <- 4
  j3 <- c(1,2,3)

  b[j2] <- lapply(b[j2], formatC, format = "f", digits = 2)
  b[j3] <- lapply(b[j3], formatC, format = "f", digits = 3)

  if ((des3$kMax > 1 && any(des3$futilityBounds[1:(des3$kMax-1)] > -8))) {
    df3b <- t(b)
    rownames(df3b) <- c("Information rate",
                        "Efficacy bounds (Z)",
                        "Futility bounds (Z)",
                        "Information")
  } else {
    df3b <- t(b[, c(1,2,4)])
    rownames(df3b) <- c("Information rate",
                        "Efficacy bounds (Z)",
                        "Information")
  }
  colnames(df3b) <- paste("Stage", seq_len(ncol(df3b)), sep=" ")

  df3c <- des3$byIntersectionBounds
  df3c[3] <- lapply(df3c[3], formatC, format = "f", digits = 3)
  colnames(df3c) <- c("Interection hypothesis", "Stage", "Boundary (Z)")

  print(df1a, ..., na.print = "" , quote = FALSE )
  print(df1b, ..., na.print = "" , quote = FALSE )
  if (des1$M > 1) {
    cat("\nBy level critical values for primary trial\n")
    print(df1c, ..., na.print = "" , quote = FALSE)
  }

  print(df2a, ..., na.print = "" , quote = FALSE )
  print(df2b, ..., na.print = "" , quote = FALSE )
  cat("\nBy hypothesis critical values for secondary trial\n")
  print(df2c, ..., na.print = "" , quote = FALSE)

  print(df3a, ..., na.print = "" , quote = FALSE )
  print(df3b, ..., na.print = "" , quote = FALSE )
  if (des3$MNew > 1) {
    cat("\nBy intersection hypothesis critical values for integrated trial\n")
    print(df3c, ..., na.print = "" , quote = FALSE)
  }
  invisible(x)
}


#' @title Print Simulation Results for Multi-Arm Multi-Stage Design
#' @description Prints the summary statistics from simulation.
#'
#' @param x The lrsim_mams object to print.
#' @param ... Ensures that all arguments starting from "..." are named.
#'
#' @return A tabular printout of the summary statistics from simulation runs.
#'
#' @keywords internal
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @export
print.lrsim_mams <- function(x, ...) {
  a <- x$overview
  M <- a$M
  k <- a$kMax

  str1 <- paste0("Multi-arm multi-stage design")
  if (a$rho1 != 0 || a$rho2 != 0) {
    str1 <- paste0(str1, " for weighted log-rank test, FH(",
                   a$rho1, ", ", a$rho2, ")")
  } else {
    str1 <- paste0(str1, " for log-rank test")
  }

  str2 <- paste0("Empirical power: ", round(a$overallReject, 4))

  str3 <- paste0("Number of active arms: ", a$M)

  str4 <- paste0("Number of looks: ", a$kMax)

  df1a <- data.frame(x = rep("", 5))
  colnames(df1a) <- NULL
  rownames(df1a) <- c(str1, str2, str3, str4, "")

  df1b <- data.frame(a$criticalValues)
  rownames(df1b) <- paste("Stage", 1:k, sep=" ")
  colnames(df1b) <- paste("Level", M:1, sep=" ")
  j <- 1:ncol(df1b)
  df1b[j] <- lapply(df1b[j], formatC, format = "f", digits = 3)

  df2a <- as.data.frame(a$cumulativeRejection)
  rownames(df2a) <- paste("Stage", 1:k, sep=" ")
  colnames(df2a) <- c(paste("Active", 1:M, sep=" "), "Overall")
  j <- 1:ncol(df2a)
  df2a[j] <- lapply(df2a[j], formatC, format = "f", digits = 3)

  df2b <- as.data.frame(a$rejectByNumber)
  rownames(df2b) <- c(paste("Stage", 1:k, sep=" "), "Overall")
  colnames(df2b) <- paste("", 0:M, sep=" ")
  j <- 1:ncol(df2b)
  df2b[j] <- lapply(df2b[j], formatC, format = "f", digits = 3)

  df2c <- a$rejectBySet
  df2c[2] <- lapply(df2c[2], formatC, format = "f", digits = 3)
  colnames(df2c) <- c("Set of active arms", "Probability of rejection")

  df3 <- as.data.frame(t(matrix(c(a$expectedNumberOfEvents,
                                  a$expectedNumberOfDropouts,
                                  a$expectedNumberOfSubjects,
                                  a$expectedStudyDuration),
                                ncol = 4)))

  rownames(df3) <- paste("Expected", c("# events         ",
                                       "# dropouts       ",
                                       "# subjects       ",
                                       "study duration   "), sep=" ")
  colnames(df3) <- c(paste("Active", 1:M, sep=" "), "Control", "Total")
  j <- 1:ncol(df3)
  df3[j] <- lapply(df3[j], formatC, format = "f", digits = 1)

  df4a <- as.data.frame(a$numberOfEvents)
  df4b <- as.data.frame(a$numberOfDropouts)
  df4c <- as.data.frame(a$numberOfSubjects)
  df4d <- as.data.frame(a$analysisTime)
  df4 <- rbind(df4a, df4b, df4c, df4d)
  rownames(df4) <- paste(rep(c("Number of events  ",
                               "Number of dropouts",
                               "Number of subjects",
                               "Analysis time     "), each = k),
                         paste("Stage", rep(1:k, times = 4)))
  colnames(df4) <- c(paste("Active", 1:M, sep=" "), "Control", "Total")

  j <- 1:ncol(df4)
  df4[j] <- lapply(df4[j], formatC, format = "f", digits = 1)
  df4

  print(df1a, ..., na.print = "" , quote = FALSE )

  cat("By level critical boundaries\n")
  print(df1b, ..., na.print = "" , quote = FALSE )
  cat("\n")

  cat("Cumulative probability of rejection by treatment\n")
  print(df2a, ..., na.print = "" , quote = FALSE )
  cat("\n")

  cat("Probability of rejection by number of active arms\n")
  print(df2b, ..., na.print = "" , quote = FALSE )
  cat("\n")

  cat("Overall probability of rejection by set of active arms\n")
  print(df2c, ..., na.print = "" , quote = FALSE)
  cat("\n")

  print(df3, ..., na.print = "" , quote = FALSE, row.names= TRUE)
  cat("\n")

  print(df4, ..., na.print = "" , quote = FALSE, row.names= TRUE)
  invisible(x)
}

