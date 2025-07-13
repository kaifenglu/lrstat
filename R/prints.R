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
  a = x$overallResults
  s = x$byStageResults
  k = a$kMax

  if (k>1) {
    str1 = paste0("Group-sequential design with ", k, " stages")
  } else {
    str1 = "Fixed design"
  }

  str2 <- paste0("theta: ", round(a$theta, 3), ", ",
                 "maximum information: ", round(a$information, 2))

  str3 <- paste0("Overall power: ",
                 round(a$overallReject, 3), ", ",
                 "overall alpha (1-sided): ",
                 round(a$alpha, 4))

  if (x$settings$typeBetaSpending != 'none' ||
      (k > 1 && any(x$byStageResults$futilityBounds[1:(k-1)] > -6))) {
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

    asf = tolower(x$settings$typeAlphaSpending)
    asfpar = round(x$settings$parameterAlphaSpending, 3)
    asfuser = round(x$settings$userAlphaSpending, 4)

    bsf = tolower(x$settings$typeBetaSpending)
    bsfpar = round(x$settings$parameterBetaSpending, 3)
    bsfuser = round(x$settings$userBetaSpending, 4)

    if (asf == "of") {
      str6 = paste0("Alpha spending: O'Brien-Fleming")
    } else if (asf == "p") {
      str6 = paste0("Alpha spending: Pocock")
    } else if (asf == "wt") {
      str6 = paste0("Alpha spending: Wang-Tsiatis(Delta = ", asfpar, ")")
    } else if (asf == "sfof") {
      str6 = paste0("Alpha spending: Lan-DeMets O'Brien-Fleming")
    } else if (asf == "sfp") {
      str6 = paste0("Alpha spending: Lan-DeMets Pocock")
    } else if (asf == "sfkd") {
      str6 = paste0("Alpha spending: KD(rho = ", asfpar, ")")
    } else if (asf == "sfhsd") {
      str6 = paste0("Alpha spending: HSD(gamma = ", asfpar, ")")
    } else if (asf == "user") {
      str6 = paste0("Alpha spending: User defined(",
                    paste(asfuser, collapse = ","), ")")
    } else {
      str6 = "Alpha spending: None"
    }

    if (bsf == "of") {
      str7 = paste0("beta spending: O'Brien-Fleming")
    } else if (bsf == "p") {
      str7 = paste0("beta spending: Pocock")
    } else if (bsf == "wt") {
      str7 = paste0("beta spending: Wang-Tsiatis(Delta = ", bsfpar, ")")
    } else if (bsf == "sfof") {
      str7 = paste0("beta spending: Lan-DeMets O'Brien-Fleming")
    } else if (bsf == "sfp") {
      str7 = paste0("beta spending: Lan-DeMets Pocock")
    } else if (bsf == "sfkd") {
      str7 = paste0("beta spending: KD(rho = ", bsfpar, ")")
    } else if (bsf == "sfhsd") {
      str7 = paste0("beta spending: HSD(gamma = ", bsfpar, ")")
    } else if (bsf == "user") {
      str7 = paste0("beta spending: User defined(",
                    paste(bsfuser, collapse = ","), ")")
    } else {
      str7 = "beta spending: None"
    }

    if (!any(is.na(x$settings$spendingTime)) &&
        !all(x$settings$spendingTime == s$informationRates)) {
      str8 = paste0("Spending time: ",
                    paste(round(x$settings$spendingTime, 3), collapse = ","))
      df1 = data.frame(x = rep("", 8))
      colnames(df1) = NULL
      rownames(df1) = c(str1, str2, str3, str4, str5,
                        paste(str6, str7, sep = ", "), str8, "")
    } else {
      df1 = data.frame(x = rep("", 7))
      colnames(df1) = NULL
      rownames(df1) = c(str1, str2, str3, str4, str5,
                        paste(str6, str7, sep = ", "), "")
    }
  } else {
    df1 = data.frame(x = rep("", 5))
    colnames(df1) = NULL
    rownames(df1) = c(str1, str2, str3, str4, "")
  }


  if (k>1) {
    b <- s[, c("informationRates", "efficacyBounds", "futilityBounds",
               "cumulativeRejection", "cumulativeFutility",
               "cumulativeAlphaSpent", "efficacyTheta", "futilityTheta",
               "efficacyP", "futilityP", "information",
               "cumulativeRejectionH0", "cumulativeFutilityH0")]

    # format number of digits after decimal for each column
    j2 <- 11
    j3 <- c(1,2,3,4,7,8,13)
    j4 <- c(5,6,9,10,12)

    b[j2] <- lapply(b[j2], formatC, format = "f", digits = 2)
    b[j3] <- lapply(b[j3], formatC, format = "f", digits = 3)
    b[j4] <- lapply(b[j4], formatC, format = "f", digits = 4)

    if (x$settings$typeBetaSpending != 'none' ||
        (k > 1 && any(x$byStageResults$futilityBounds[1:(k-1)] > -6))) {
      df = t(b)
      rownames(df) = c("Information rate",
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
      df = t(b[,c(1,2,4,6,7,9,11)])
      rownames(df) = c("Information rate",
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

    df = t(b)

    rownames(df) = c("Efficacy boundary (Z)",
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
  a = x$overallResults
  s = x$byStageResults
  k = a$kMax

  if (k>1) {
    str1 = paste0("Group-sequential design with ", k, " stages")
  } else {
    str1 = "Fixed design"
  }

  str1 <- paste(str1, "for equivalence test")

  str2 <- paste0("Lower equivalence limit: ", round(a$thetaLower, 3), ", ",
                 "upper equivalence limit: ", round(a$thetaUpper, 3), ", ",
                 "parameter value: ", round(a$theta, 3))

  str3 <- paste0("Maximum information: ", round(a$information, 2))

  str4 <- paste0("Overall power: ",
                 round(a$overallReject, 3), ", ",
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

    asf = tolower(x$settings$typeAlphaSpending)
    asfpar = round(x$settings$parameterAlphaSpending, 3)
    asfuser = round(x$settings$userAlphaSpending, 4)

    if (asf == "of") {
      str6 = paste0("Alpha spending: O'Brien-Fleming")
    } else if (asf == "p") {
      str6 = paste0("Alpha spending: Pocock")
    } else if (asf == "wt") {
      str6 = paste0("Alpha spending: Wang-Tsiatis(Delta = ", asfpar, ")")
    } else if (asf == "sfof") {
      str6 = paste0("Alpha spending: Lan-DeMets O'Brien-Fleming")
    } else if (asf == "sfp") {
      str6 = paste0("Alpha spending: Lan-DeMets Pocock")
    } else if (asf == "sfkd") {
      str6 = paste0("Alpha spending: KD(rho = ", asfpar, ")")
    } else if (asf == "sfhsd") {
      str6 = paste0("Alpha spending: HSD(gamma = ", asfpar, ")")
    } else if (asf == "user") {
      str6 = paste0("Alpha spending: User defined(",
                    paste(asfuser, collapse = ","), ")")
    } else {
      str6 = "Alpha spending: None"
    }

    if (!any(is.na(x$settings$spendingTime)) &&
        !all(x$settings$spendingTime == s$informationRates)) {
      str7 = paste0("Spending time: ",
                    paste(round(x$settings$spendingTime, 3), collapse = ","))
      df1 = data.frame(x = rep("", 8))
      colnames(df1) = NULL
      rownames(df1) = c(str1, str2, str3, str4, str5, str6, str7, "")
    } else {
      df1 = data.frame(x = rep("", 7))
      colnames(df1) = NULL
      rownames(df1) = c(str1, str2, str3, str4, str5, str6, "")
    }
  } else {
    df1 = data.frame(x = rep("", 5))
    colnames(df1) = NULL
    rownames(df1) = c(str1, str2, str3, str4, "")
  }


  if (k>1) {
    b <- s[, c("informationRates", "efficacyBounds",
               "cumulativeRejection", "cumulativeAlphaSpent",
               "cumulativeAttainedAlphaH10", "cumulativeAttainedAlphaH20",
               "efficacyThetaLower", "efficacyThetaUpper",
               "efficacyP", "information")]

    # format number of digits after decimal for each column
    j2 <- 10
    j3 <- c(1,2,3,7,8)
    j4 <- c(4,5,6,9)

    b[j2] <- lapply(b[j2], formatC, format = "f", digits = 2)
    b[j3] <- lapply(b[j3], formatC, format = "f", digits = 3)
    b[j4] <- lapply(b[j4], formatC, format = "f", digits = 4)

    df = t(b)
    rownames(df) = c("Information rate",
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

    df = t(b)

    rownames(df) = c("Boundary for each 1-sided test (Z)",
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
  des1 = x$primaryTrial

  str1 = "Primary trial:"
  str2 = paste0("Group-sequential design with ", des1$kMax, " stages")
  str3 = paste0("Interim adaptation look: ",  des1$L, ", ",
                "z-statistic value: ", round(des1$zL, 3))

  str4 = paste0("Conditional type I error: ",
                round(des1$conditionalAlpha, 4))

  if (!is.na(des1$conditionalPower)) {
    str5 = paste0("Conditional power: ", round(des1$conditionalPower, 3),
                  ", ", "predictive power: ",
                  round(des1$predictivePower, 3))

    str6 = paste0("Muller & Schafer method for secondary trial: ",
                  des1$MullerSchafer)

    df1a = data.frame(x = rep("", 7))
    colnames(df1a) = NULL
    rownames(df1a) = c(str1, str2, str3, str4, str5, str6, "")
  } else {
    str5 = paste0("Muller & Schafer method for secondary trial: ",
                  des1$MullerSchafer)

    df1a = data.frame(x = rep("", 6))
    colnames(df1a) = NULL
    rownames(df1a) = c(str1, str2, str3, str4, str5, "")
  }

  b <- data.frame(informationRates = des1$informationRates,
                  efficacyBounds = des1$efficacyBounds,
                  futilityBounds = des1$futilityBounds)

  b[1:3] <- lapply(b[1:3], formatC, format = "f", digits = 3)

  if (!all(des1$futilityBounds[1:(des1$kMax-1)] == -6)) {
    df1b = t(b)
    rownames(df1b) = c("Information rate",
                       "Efficacy boundary (Z)",
                       "Futility boundary (Z)")
  } else {
    df1b = t(b[1:2])
    rownames(df1b) = c("Information rate",
                       "Efficacy boundary (Z)")
  }

  colnames(df1b) <- paste("Stage", seq_len(ncol(df1b)), sep=" ")


  des2 = x$secondaryTrial
  a = des2$overallResults
  s = des2$byStageResults
  k = a$kMax

  str1 = "Secondary trial:"

  if (k>1) {
    str2 = paste0("Group-sequential design with ", k, " stages")
  } else {
    str2 = "Fixed design"
  }

  str3 <- paste0("theta: ", round(a$theta, 3), ", ",
                 "maximum information: ", round(a$information, 2))

  str4 <- paste0("Overall power: ",
                 round(a$overallReject, 3), ", ",
                 "overall significance level (1-sided): ",
                 round(a$alpha, 4))

  str5 <- paste0("Drift parameter: ", round(a$drift, 3), ", ",
                 "inflation factor: ", round(a$inflationFactor, 3))

  if (k>1) {
    str6 <- paste0("Expected information under H1: ",
                   round(a$expectedInformationH1, 2), ", ",
                   "Expected information under H0: ",
                   round(a$expectedInformationH0, 2))
    df2a = data.frame(x = rep("", 7))
    colnames(df2a) = NULL
    rownames(df2a) = c(str1, str2, str3, str4, str5, str6, "")
  } else {
    df2a = data.frame(x = rep("", 6))
    colnames(df2a) = NULL
    rownames(df2a) = c(str1, str2, str3, str4, str5, "")
  }


  if (k>1) {
    b <- s[, c("informationRates", "efficacyBounds", "futilityBounds",
               "cumulativeRejection", "cumulativeFutility",
               "cumulativeAlphaSpent", "efficacyTheta", "futilityTheta",
               "efficacyP", "futilityP", "information")]

    # format number of digits after decimal for each column
    j2 <- 11
    j3 <- c(1,2,3,4,5,7,8)
    j4 <- c(6,9,10)

    b[j2] <- lapply(b[j2], formatC, format = "f", digits = 2)
    b[j3] <- lapply(b[j3], formatC, format = "f", digits = 3)
    b[j4] <- lapply(b[j4], formatC, format = "f", digits = 4)

    if (des2$settings$typeBetaSpending != 'none' ||
        (k > 1 && any(des2$byStageResults$futilityBounds[1:(k-1)] > -6))) {
      df2b = t(b)
      rownames(df2b) = c("Information rate",
                         "Efficacy boundary (Z)",
                         "Futility boundary (Z)",
                         "Cumulative rejection",
                         "Cumulative futility",
                         "Cumulative alpha spent",
                         "Efficacy boundary (theta)",
                         "Futility boundary (theta)",
                         "Efficacy boundary (p)",
                         "Futility boundary (p)",
                         "Information")

    } else {
      df2b = t(b[,c(1,2,4,6,7,9,11)])
      rownames(df2b) = c("Information rate",
                         "Efficacy boundary (Z)",
                         "Cumulative rejection",
                         "Cumulative alpha spent",
                         "Efficacy boundary (theta)",
                         "Efficacy boundary (p)",
                         "Information")
    }

    colnames(df2b) <- paste("Stage", seq_len(ncol(df2b)), sep=" ")
  } else {
    b <- s[, c("efficacyBounds", "efficacyTheta", "efficacyP")]

    # format number of digits after decimal for each column
    j3 <- c(1,2)
    j4 <- 3

    b[j3] <- lapply(b[j3], formatC, format = "f", digits = 3)
    b[j4] <- lapply(b[j4], formatC, format = "f", digits = 4)

    df2b = t(b)

    rownames(df2b) = c("Efficacy boundary (Z)",
                       "Efficacy boundary (theta)",
                       "Efficacy boundary (p)")

    colnames(df2b) <- NA
  }

  print(df1a, ..., na.print = "" , quote = FALSE )
  print(df1b, ..., na.print = "" , quote = FALSE )
  print(df2a, ..., na.print = "" , quote = FALSE )
  print(df2b, ..., na.print = "" , quote = FALSE )
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
  a = x$overallResults
  s = x$byStageResults
  k = a$kMax

  if (k>1) {
    str1 = paste0("Group-sequential design with ", k, " stages")
  } else {
    str1 = "Fixed design"
  }

  if (a$rho1 != 0 || a$rho2 != 0) {
    str1 <- paste0(str1, " for weighted log-rank test, FH(",
                   a$rho1, ", ", a$rho2, ")")
  } else {
    str1 <- paste0(str1, " for log-rank test")
  }

  str2 <- paste0("Overall power: ",
                 round(a$overallReject, 3), ", ",
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
    asf = tolower(x$settings$typeAlphaSpending)
    asfpar = round(x$settings$parameterAlphaSpending, 3)
    asfuser = round(x$settings$userAlphaSpending, 4)

    bsf = tolower(x$settings$typeBetaSpending)
    bsfpar = round(x$settings$parameterBetaSpending, 3)

    if (asf == "of") {
      str10 = paste0("Alpha spending: O'Brien-Fleming")
    } else if (asf == "p") {
      str10 = paste0("Alpha spending: Pocock")
    } else if (asf == "wt") {
      str10 = paste0("Alpha spending: Wang-Tsiatis(Delta = ", asfpar, ")")
    } else if (asf == "sfof") {
      str10 = paste0("Alpha spending: Lan-DeMets O'Brien-Fleming")
    } else if (asf == "sfp") {
      str10 = paste0("Alpha spending: Lan-DeMets Pocock")
    } else if (asf == "sfkd") {
      str10 = paste0("Alpha spending: KD(rho = ", asfpar, ")")
    } else if (asf == "sfhsd") {
      str10 = paste0("Alpha spending: HSD(gamma = ", asfpar, ")")
    } else if (asf == "user") {
      str10 = paste0("Alpha spending: User defined(",
                     paste(asfuser, collapse = ","), ")")
    } else {
      str10 = "Alpha spending: None"
    }

    if (bsf == "of") {
      str11 = paste0("beta spending: O'Brien-Fleming")
    } else if (bsf == "p") {
      str11 = paste0("beta spending: Pocock")
    } else if (bsf == "wt") {
      str11 = paste0("beta spending: Wang-Tsiatis(Delta = ", bsfpar, ")")
    } else if (bsf == "sfof") {
      str11 = paste0("beta spending: Lan-DeMets O'Brien-Fleming")
    } else if (bsf == "sfp") {
      str11 = paste0("beta spending: Lan-DeMets Pocock")
    } else if (bsf == "sfkd") {
      str11 = paste0("beta spending: KD(rho = ", bsfpar, ")")
    } else if (bsf == "sfhsd") {
      str11 = paste0("beta spending: HSD(gamma = ", bsfpar, ")")
    } else if (bsf == "user") {
      str11 = paste0("beta spending: User defined")
    } else {
      str11 = "beta spending: None"
    }

    if (!any(is.na(x$settings$spendingTime)) &&
        !all(x$settings$spendingTime == s$informationRates)) {
      str12 = paste0("Spending time: ",
                     paste(round(x$settings$spendingTime, 3), collapse = ","))
      df1 = data.frame(x = rep("", 12))
      colnames(df1) = NULL
      rownames(df1) = c(str1, str2, str3, str4, str5, str6, str7,
                        str8, str9, paste(str10, str11, sep = ", "),
                        str12, "")
    } else {
      df1 = data.frame(x = rep("", 11))
      colnames(df1) = NULL
      rownames(df1) = c(str1, str2, str3, str4, str5, str6, str7,
                        str8, str9, paste(str10, str11, sep = ", "), "")
    }
  } else {
    df1 = data.frame(x = rep("", 10))
    colnames(df1) = NULL
    rownames(df1) = c(str1, str2, str3, str4, str5, str6, str7,
                      str8, str9, "")
  }

  if (k>1) {
    if (x$settings$estimateHazardRatio) {
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
          (k > 1 && any(x$byStageResults$futilityBounds[1:(k-1)] > -6))) {
        df = t(b)
        rownames(df) = c("Information rate",
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
        df = t(b[,c(1,2,4,6,7,8,9,10,11,13,15,16)])
        rownames(df) = c("Information rate",
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

    } else {
      b <- s[, c("informationRates", "efficacyBounds", "futilityBounds",
                 "cumulativeRejection", "cumulativeFutility",
                 "cumulativeAlphaSpent",
                 "numberOfEvents", "numberOfDropouts", "numberOfSubjects",
                 "analysisTime", "efficacyP", "futilityP", "information")]

      # format number of digits after decimal for each column
      j1 <- c(7,8,9,10)
      j2 <- 13
      j3 <- c(1,2,3)
      j4 <- c(4,5,6,11,12)

      b[j1] <- lapply(b[j1], formatC, format = "f", digits = 1)
      b[j2] <- lapply(b[j2], formatC, format = "f", digits = 2)
      b[j3] <- lapply(b[j3], formatC, format = "f", digits = 3)
      b[j4] <- lapply(b[j4], formatC, format = "f", digits = 4)

      if (x$settings$typeBetaSpending != 'none' ||
          (k > 1 && any(x$byStageResults$futilityBounds[1:(k-1)] > -6))) {
        df = t(b)
        rownames(df) = c("Information rate",
                         "Efficacy boundary (Z)",
                         "Futility boundary (Z)",
                         "Cumulative rejection",
                         "Cumulative futility",
                         "Cumulative alpha spent",
                         "Number of events",
                         "Number of dropouts",
                         "Number of subjects",
                         "Analysis time",
                         "Efficacy boundary (p)",
                         "Futility boundary (p)",
                         "Information")

      } else {
        df = t(b[,c(1,2,4,6,7,8,9,10,11,13)])
        rownames(df) = c("Information rate",
                         "Efficacy boundary (Z)",
                         "Cumulative rejection",
                         "Cumulative alpha spent",
                         "Number of events",
                         "Number of dropouts",
                         "Number of subjects",
                         "Analysis time",
                         "Efficacy boundary (p)",
                         "Information")
      }

    }

    colnames(df) <- paste("Stage", seq_len(ncol(df)), sep=" ")
  } else {
    if (x$settings$estimateHazardRatio) {
      b <- s[, c("efficacyBounds", "efficacyHR", "efficacyP",
                 "HR")]

      # format number of digits after decimal for each column
      j3 <- c(1,2,4)
      j4 <- 3

      b[j3] <- lapply(b[j3], formatC, format = "f", digits = 3)
      b[j4] <- lapply(b[j4], formatC, format = "f", digits = 4)

      df = t(b)

      rownames(df) = c("Efficacy boundary (Z)",
                       "Efficacy boundary (HR)",
                       "Efficacy boundary (p)",
                       "HR")
    } else {
      b <- s[, c("efficacyBounds", "efficacyP")]

      # format number of digits after decimal for each column
      j3 <- 1
      j4 <- 2

      b[j3] <- lapply(b[j3], formatC, format = "f", digits = 3)
      b[j4] <- lapply(b[j4], formatC, format = "f", digits = 4)

      df = t(b)

      rownames(df) = c("Efficacy boundary (Z)",
                       "Efficacy boundary (p)")
    }

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
  a = x$overview
  k = a$kMax

  if (k>1) {
    str1 = paste0("Group-sequential design with ", k, " stages")
  } else {
    str1 = "Fixed design"
  }

  if (a$rho1 != 0 || a$rho2 != 0) {
    str1 <- paste0(str1, " for weighted log-rank test, FH(",
                   a$rho1, ", ", a$rho2, ")")
  } else {
    str1 = paste0(str1, " for log-rank test")
  }

  str2 <- paste0("Overall power: ", round(a$overallReject, 3))

  str3 <- paste0("Expected # events: ",
                 round(a$expectedNumberOfEvents, 1))

  str4 <- paste0("Expected # dropouts: ",
                 round(a$expectedNumberOfDropouts, 1))

  str5 <- paste0("Expected # subjects: ",
                 round(a$expectedNumberOfSubjects, 1))

  str6 <- paste0("Expected study duration: ",
                 round(a$expectedStudyDuration, 1))

  str7 <- paste0("Accrual duration: ",
                 round(a$accrualDuration, 1), ", ",
                 "fixed follow-up: ", a$fixedFollowup)

  df1 = data.frame(x = rep("", 8))
  colnames(df1) = NULL
  rownames(df1) = c(str1, str2, str3, str4, str5, str6, str7, "")


  if (k>1) {
    b <- data.frame(a$cumulativeRejection,
                    a$cumulativeFutility,
                    a$numberOfEvents,
                    a$numberOfDropouts,
                    a$numberOfSubjects,
                    a$analysisTime)

    # format number of digits after decimal for each column
    j1 <- c(3,4,5,6)
    j3 <- c(1,2)

    b[j1] <- lapply(b[j1], formatC, format = "f", digits = 1)
    b[j3] <- lapply(b[j3], formatC, format = "f", digits = 3)

    df = t(b)
    rownames(df) = c("Cumulative rejection",
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
  a = x$overallResults
  s = x$byStageResults
  k = a$kMax

  if (k>1) {
    str1 = paste0("Group-sequential design with ", k, " stages")
  } else {
    str1 = "Fixed design"
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
                 round(a$overallReject, 3), ", ",
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
    asf = tolower(x$settings$typeAlphaSpending)
    asfpar = round(x$settings$parameterAlphaSpending, 3)
    asfuser = round(x$settings$userAlphaSpending, 4)

    bsf = tolower(x$settings$typeBetaSpending)
    bsfpar = round(x$settings$parameterBetaSpending, 3)

    if (asf == "of") {
      str14 = paste0("Alpha spending: O'Brien-Fleming")
    } else if (asf == "p") {
      str14 = paste0("Alpha spending: Pocock")
    } else if (asf == "wt") {
      str14 = paste0("Alpha spending: Wang-Tsiatis(Delta = ", asfpar, ")")
    } else if (asf == "sfof") {
      str14 = paste0("Alpha spending: Lan-DeMets O'Brien-Fleming")
    } else if (asf == "sfp") {
      str14 = paste0("Alpha spending: Lan-DeMets Pocock")
    } else if (asf == "sfkd") {
      str14 = paste0("Alpha spending: KD(rho = ", asfpar, ")")
    } else if (asf == "sfhsd") {
      str14 = paste0("Alpha spending: HSD(gamma = ", asfpar, ")")
    } else if (asf == "user") {
      str14 = paste0("Alpha spending: User defined(",
                     paste(asfuser, collapse = ","), ")")
    } else {
      str14 = "Alpha spending: None"
    }

    if (bsf == "of") {
      str15 = paste0("beta spending: O'Brien-Fleming")
    } else if (bsf == "p") {
      str15 = paste0("beta spending: Pocock")
    } else if (bsf == "wt") {
      str15 = paste0("beta spending: Wang-Tsiatis(Delta = ", bsfpar, ")")
    } else if (bsf == "sfof") {
      str15 = paste0("beta spending: Lan-DeMets O'Brien-Fleming")
    } else if (bsf == "sfp") {
      str15 = paste0("beta spending: Lan-DeMets Pocock")
    } else if (bsf == "sfkd") {
      str15 = paste0("beta spending: KD(rho = ", bsfpar, ")")
    } else if (bsf == "sfhsd") {
      str15 = paste0("beta spending: HSD(gamma = ", bsfpar, ")")
    } else if (bsf == "user") {
      str15 = paste0("beta spending: User defined")
    } else {
      str15 = "beta spending: None"
    }

    if (!any(is.na(x$settings$spendingTime)) &&
        !all(x$settings$spendingTime == s$informationRates)) {
      str16 = paste0("Spending time: ",
                     paste(round(x$settings$spendingTime, 3), collapse = ","))

      if (length(x$settings$stratumFraction) > 1) {
        df1 = data.frame(x = rep("", 17))
        colnames(df1) = NULL
        rownames(df1) = c(str1, str2, str3a, str3, str4, str5, str6, str7,
                          str8, str9, str10, str11, str12, str13,
                          paste(str14, str15, sep = ", "), str16, "")

      } else {
        df1 = data.frame(x = rep("", 16))
        colnames(df1) = NULL
        rownames(df1) = c(str1, str2, str3, str4, str5, str6, str7,
                          str8, str9, str10, str11, str12, str13,
                          paste(str14, str15, sep = ", "), str16, "")
      }
    } else {
      if (length(x$settings$stratumFraction) > 1) {
        df1 = data.frame(x = rep("", 16))
        colnames(df1) = NULL
        rownames(df1) = c(str1, str2, str3a, str3, str4, str5, str6, str7,
                          str8, str9, str10, str11, str12, str13,
                          paste(str14, str15, sep = ", "), "")
      } else {
        df1 = data.frame(x = rep("", 15))
        colnames(df1) = NULL
        rownames(df1) = c(str1, str2, str3, str4, str5, str6, str7,
                          str8, str9, str10, str11, str12, str13,
                          paste(str14, str15, sep = ", "), "")
      }
    }
  } else {
    if (length(x$settings$stratumFraction) > 1) {
      df1 = data.frame(x = rep("", 15))
      colnames(df1) = NULL
      rownames(df1) = c(str1, str2, str3a, str3, str4, str5, str6, str7,
                        str8, str9, str10, str11, str12, str13, "")
    } else {
      df1 = data.frame(x = rep("", 14))
      colnames(df1) = NULL
      rownames(df1) = c(str1, str2, str3, str4, str5, str6, str7,
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
        (k > 1 && any(x$byStageResults$futilityBounds[1:(k-1)] > -6))) {
      df = t(b)
      rownames(df) = c("Information rate",
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
      df = t(b[,c(1,2,4,6,7,8,9,10,11,12,14,16)])
      rownames(df) = c("Information rate",
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

    df = t(b)

    rownames(df) = c("Efficacy boundary (Z)",
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
  a = x$overallResults
  s = x$byStageResults
  k = a$kMax

  if (k>1) {
    str1 = paste0("Group-sequential design with ", k, " stages")
  } else {
    str1 = "Fixed design"
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
                 round(a$overallReject, 3), ", ",
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
    asf = tolower(x$settings$typeAlphaSpending)
    asfpar = round(x$settings$parameterAlphaSpending, 3)
    asfuser = round(x$settings$userAlphaSpending, 4)

    if (asf == "of") {
      str13 = paste0("Alpha spending: O'Brien-Fleming")
    } else if (asf == "p") {
      str13 = paste0("Alpha spending: Pocock")
    } else if (asf == "wt") {
      str13 = paste0("Alpha spending: Wang-Tsiatis(Delta = ", asfpar, ")")
    } else if (asf == "sfof") {
      str13 = paste0("Alpha spending: Lan-DeMets O'Brien-Fleming")
    } else if (asf == "sfp") {
      str13 = paste0("Alpha spending: Lan-DeMets Pocock")
    } else if (asf == "sfkd") {
      str13 = paste0("Alpha spending: KD(rho = ", asfpar, ")")
    } else if (asf == "sfhsd") {
      str13 = paste0("Alpha spending: HSD(gamma = ", asfpar, ")")
    } else if (asf == "user") {
      str13 = paste0("Alpha spending: User defined(",
                     paste(asfuser, collapse = ","), ")")
    } else {
      str13 = "Alpha spending: None"
    }

    if (!any(is.na(x$settings$spendingTime)) &&
        !all(x$settings$spendingTime == s$informationRates)) {
      str14 = paste0("Spending time: ",
                     paste(round(x$settings$spendingTime, 3), collapse = ","))
      if (length(x$settings$stratumFraction) > 1) {
        df1 = data.frame(x = rep("", 16))
        colnames(df1) = NULL
        rownames(df1) = c(str1, str2, str3a, str3, str4, str5, str6, str7,
                          str8, str9, str10, str11, str12,
                          str13, str14, "")
      } else {
        df1 = data.frame(x = rep("", 15))
        colnames(df1) = NULL
        rownames(df1) = c(str1, str2, str3, str4, str5, str6, str7,
                          str8, str9, str10, str11, str12,
                          str13, str14, "")
      }
    } else {
      if (length(x$settings$stratumFraction) > 1) {
        df1 = data.frame(x = rep("", 15))
        colnames(df1) = NULL
        rownames(df1) = c(str1, str2, str3a, str3, str4, str5, str6, str7,
                          str8, str9, str10, str11, str12,
                          str13, "")
      } else {
        df1 = data.frame(x = rep("", 14))
        colnames(df1) = NULL
        rownames(df1) = c(str1, str2, str3, str4, str5, str6, str7,
                          str8, str9, str10, str11, str12,
                          str13, "")
      }
    }
  } else {
    if (length(x$settings$stratumFraction) > 1) {
      df1 = data.frame(x = rep("", 14))
      colnames(df1) = NULL
      rownames(df1) = c(str1, str2, str3a, str3, str4, str5, str6, str7,
                        str8, str9, str10, str11, str12, "")
    } else {
      df1 = data.frame(x = rep("", 13))
      colnames(df1) = NULL
      rownames(df1) = c(str1, str2, str3, str4, str5, str6, str7,
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
    j3 <- c(1,2,3,12,13)
    j4 <- c(4,5,6,14)

    b[j1] <- lapply(b[j1], formatC, format = "f", digits = 1)
    b[j2] <- lapply(b[j2], formatC, format = "f", digits = 2)
    b[j3] <- lapply(b[j3], formatC, format = "f", digits = 3)
    b[j4] <- lapply(b[j4], formatC, format = "f", digits = 4)

    df = t(b)
    rownames(df) = c("Information rate",
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

    df = t(b)

    rownames(df) = c("Boundary for each 1-sided test (Z)",
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
  a = x$overallResults
  s = x$byStageResults
  k = a$kMax

  if (k>1) {
    str1 = paste0("Group-sequential design with ", k, " stages")
  } else {
    str1 = "Fixed design"
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
                 round(a$overallReject, 3), ", ",
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
    asf = tolower(x$settings$typeAlphaSpending)
    asfpar = round(x$settings$parameterAlphaSpending, 3)
    asfuser = round(x$settings$userAlphaSpending, 4)

    bsf = tolower(x$settings$typeBetaSpending)
    bsfpar = round(x$settings$parameterBetaSpending, 3)

    if (asf == "of") {
      str12 = paste0("Alpha spending: O'Brien-Fleming")
    } else if (asf == "p") {
      str12 = paste0("Alpha spending: Pocock")
    } else if (asf == "wt") {
      str12 = paste0("Alpha spending: Wang-Tsiatis(Delta = ", asfpar, ")")
    } else if (asf == "sfof") {
      str12 = paste0("Alpha spending: Lan-DeMets O'Brien-Fleming")
    } else if (asf == "sfp") {
      str12 = paste0("Alpha spending: Lan-DeMets Pocock")
    } else if (asf == "sfkd") {
      str12 = paste0("Alpha spending: KD(rho = ", asfpar, ")")
    } else if (asf == "sfhsd") {
      str12 = paste0("Alpha spending: HSD(gamma = ", asfpar, ")")
    } else if (asf == "user") {
      str12 = paste0("Alpha spending: User defined(",
                     paste(asfuser, collapse = ","), ")")
    } else {
      str12 = "Alpha spending: None"
    }

    if (bsf == "of") {
      str13 = paste0("beta spending: O'Brien-Fleming")
    } else if (bsf == "p") {
      str13 = paste0("beta spending: Pocock")
    } else if (bsf == "wt") {
      str13 = paste0("beta spending: Wang-Tsiatis(Delta = ", bsfpar, ")")
    } else if (bsf == "sfof") {
      str13 = paste0("beta spending: Lan-DeMets O'Brien-Fleming")
    } else if (bsf == "sfp") {
      str13 = paste0("beta spending: Lan-DeMets Pocock")
    } else if (bsf == "sfkd") {
      str13 = paste0("beta spending: KD(rho = ", bsfpar, ")")
    } else if (bsf == "sfhsd") {
      str13 = paste0("beta spending: HSD(gamma = ", bsfpar, ")")
    } else if (bsf == "user") {
      str13 = paste0("beta spending: User defined")
    } else {
      str13 = "beta spending: None"
    }

    if (!any(is.na(x$settings$spendingTime)) &&
        !all(x$settings$spendingTime == s$informationRates)) {
      str14 = paste0("Spending time: ",
                     paste(round(x$settings$spendingTime, 3), collapse = ","))
      df1 = data.frame(x = rep("", 14))
      colnames(df1) = NULL
      rownames(df1) = c(str1, str2, str3, str4, str5, str6, str7,
                        str8, str9, str10, str11,
                        paste(str12, str13, sep = ", "), str14, "")
    } else {
      df1 = data.frame(x = rep("", 13))
      colnames(df1) = NULL
      rownames(df1) = c(str1, str2, str3, str4, str5, str6, str7,
                        str8, str9, str10, str11,
                        paste(str12, str13, sep = ", "), "")
    }
  } else {
    df1 = data.frame(x = rep("", 12))
    colnames(df1) = NULL
    rownames(df1) = c(str1, str2, str3, str4, str5, str6, str7,
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
    j3 <- c(1,2,3,4,5)
    j4 <- c(6,12,13,14,15)

    b[j1] <- lapply(b[j1], formatC, format = "f", digits = 1)
    b[j2] <- lapply(b[j2], formatC, format = "f", digits = 2)
    b[j3] <- lapply(b[j3], formatC, format = "f", digits = 3)
    b[j4] <- lapply(b[j4], formatC, format = "f", digits = 4)

    if (x$settings$typeBetaSpending != 'none' ||
        (k > 1 && any(x$byStageResults$futilityBounds[1:(k-1)] > -6))) {
      df = t(b)
      rownames(df) = c("Information rate",
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
      df = t(b[,c(1,2,4,6,7,8,9,10,11,12,14,16)])
      rownames(df) = c("Information rate",
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

    df = t(b)

    rownames(df) = c("Efficacy boundary (Z)",
                     "Efficacy boundary (rate)",
                     "Efficacy boundary (p)")
    colnames(df) <- NA
  }

  print(df1, ..., na.print = "" , quote = FALSE )
  print(df, ..., na.print = "" , quote = FALSE )
  invisible(x)
}


#' @title Print Power and Sample Size Results for One-Sample Mean
#' @description Prints the summary statistics from power calculation of
#' one-sample mean.
#'
#' @param x The designOneMean object to print.
#' @param ... Ensures that all arguments starting from "..." are named.
#'
#' @return A tabular printout of the design elements.
#'
#' @keywords internal
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @export
print.designOneMean <- function(x, ...) {
  a = x$overallResults
  s = x$byStageResults
  k = a$kMax

  if (k>1) {
    str1 = paste0("Group-sequential design with ", k, " stages")
  } else {
    str1 = "Fixed design"
  }

  str1 <- paste(str1, "for one-sample mean")

  str2 <- paste0("Mean under H0: ",
                 round(a$meanH0, 3), ", ",
                 "mean under H1: ",
                 round(a$mean, 3), ", ",
                 "standard deviation: ", round(a$stDev, 3))

  str3 <- paste0("Overall power: ",
                 round(a$overallReject, 3), ", ",
                 "overall alpha (1-sided): ",
                 round(a$alpha, 4))

  if (x$settings$typeBetaSpending != 'none' ||
      (k > 1 && any(x$byStageResults$futilityBounds[1:(k-1)] > -6))) {
    str3 <- paste0(str3, ", ",
                   "attained alpha: ", round(a$attainedAlpha, 4))
  }

  str4 <- paste0("Drift parameter: ", round(a$drift, 3), ", ",
                 "inflation factor: ", round(a$inflationFactor, 3))


  if (k>1) {
    str5 <- paste0("Maximum information: ",
                   round(a$information, 2), ", ",
                   "expected under H1: ",
                   round(a$expectedInformationH1, 2), ", ",
                   "expected under H0: ",
                   round(a$expectedInformationH0, 2))

    str6 <- paste0("Maximum # subjects: ",
                   round(a$numberOfSubjects, 1), ", ",
                   "expected under H1: ",
                   round(a$expectedNumberOfSubjectsH1, 1), ", ",
                   "expected under H0: ",
                   round(a$expectedNumberOfSubjectsH0, 1))

  } else {
    str5 <- paste0("Information: ",
                   round(a$information, 2))

    str6 <- paste0("Number of subjects: ",
                   round(a$numberOfSubjects, 1))
  }

  if (k > 1) {
    asf = tolower(x$settings$typeAlphaSpending)
    asfpar = round(x$settings$parameterAlphaSpending, 3)
    asfuser = round(x$settings$userAlphaSpending, 4)

    bsf = tolower(x$settings$typeBetaSpending)
    bsfpar = round(x$settings$parameterBetaSpending, 3)
    bsfuser = round(x$settings$userBetaSpending, 4)

    if (asf == "of") {
      str7 = paste0("Alpha spending: O'Brien-Fleming")
    } else if (asf == "p") {
      str7 = paste0("Alpha spending: Pocock")
    } else if (asf == "wt") {
      str7 = paste0("Alpha spending: Wang-Tsiatis(Delta = ", asfpar, ")")
    } else if (asf == "sfof") {
      str7 = paste0("Alpha spending: Lan-DeMets O'Brien-Fleming")
    } else if (asf == "sfp") {
      str7 = paste0("Alpha spending: Lan-DeMets Pocock")
    } else if (asf == "sfkd") {
      str7 = paste0("Alpha spending: KD(rho = ", asfpar, ")")
    } else if (asf == "sfhsd") {
      str7 = paste0("Alpha spending: HSD(gamma = ", asfpar, ")")
    } else if (asf == "user") {
      str7 = paste0("Alpha spending: User defined(",
                    paste(asfuser, collapse = ","), ")")
    } else {
      str7 = "Alpha spending: None"
    }

    if (bsf == "of") {
      str8 = paste0("beta spending: O'Brien-Fleming")
    } else if (bsf == "p") {
      str8 = paste0("beta spending: Pocock")
    } else if (bsf == "wt") {
      str8 = paste0("beta spending: Wang-Tsiatis(Delta = ", bsfpar, ")")
    } else if (bsf == "sfof") {
      str8 = paste0("beta spending: Lan-DeMets O'Brien-Fleming")
    } else if (bsf == "sfp") {
      str8 = paste0("beta spending: Lan-DeMets Pocock")
    } else if (bsf == "sfkd") {
      str8 = paste0("beta spending: KD(rho = ", bsfpar, ")")
    } else if (bsf == "sfhsd") {
      str8 = paste0("beta spending: HSD(gamma = ", bsfpar, ")")
    } else if (bsf == "user") {
      str8 = paste0("beta spending: User defined(",
                    paste(bsfuser, collapse = ","), ")")
    } else {
      str8 = "beta spending: None"
    }

    if (!any(is.na(x$settings$spendingTime)) &&
        !all(x$settings$spendingTime == s$informationRates)) {
      str9 = paste0("Spending time: ",
                    paste(round(x$settings$spendingTime, 3), collapse = ","))
      df1 = data.frame(x = rep("", 9))
      colnames(df1) = NULL
      rownames(df1) = c(str1, str2, str3, str4, str5, str6,
                        paste(str7, str8, sep = ", "), str9, "")
    } else {
      df1 = data.frame(x = rep("", 8))
      colnames(df1) = NULL
      rownames(df1) = c(str1, str2, str3, str4, str5, str6,
                        paste(str7, str8, sep = ", "), "")
    }
  } else {
    df1 = data.frame(x = rep("", 7))
    colnames(df1) = NULL
    rownames(df1) = c(str1, str2, str3, str4, str5, str6, "")
  }

  if (k>1) {
    b <- s[, c("informationRates", "efficacyBounds", "futilityBounds",
               "cumulativeRejection", "cumulativeFutility",
               "cumulativeAlphaSpent", "numberOfSubjects",
               "efficacyMean", "futilityMean",
               "efficacyP", "futilityP", "information",
               "cumulativeRejectionH0", "cumulativeFutilityH0")]

    # format number of digits after decimal for each column
    j1 <- 7
    j2 <- 12
    j3 <- c(1,2,3,4,8,9,14)
    j4 <- c(5,6,10,11,13)

    b[j1] <- lapply(b[j1], formatC, format = "f", digits = 1)
    b[j2] <- lapply(b[j2], formatC, format = "f", digits = 2)
    b[j3] <- lapply(b[j3], formatC, format = "f", digits = 3)
    b[j4] <- lapply(b[j4], formatC, format = "f", digits = 4)

    if (x$settings$typeBetaSpending != 'none' ||
        (k > 1 && any(x$byStageResults$futilityBounds[1:(k-1)] > -6))) {
      df = t(b)
      rownames(df) = c("Information rate",
                       "Efficacy boundary (Z)",
                       "Futility boundary (Z)",
                       "Cumulative rejection",
                       "Cumulative futility",
                       "Cumulative alpha spent",
                       "Number of subjects",
                       "Efficacy boundary (mean)",
                       "Futility boundary (mean)",
                       "Efficacy boundary (p)",
                       "Futility boundary (p)",
                       "Information",
                       "Cumulative rejection under H0",
                       "Cumulative futility under H0")

    } else {
      df = t(b[,c(1,2,4,6,7,8,10,12)])
      rownames(df) = c("Information rate",
                       "Efficacy boundary (Z)",
                       "Cumulative rejection",
                       "Cumulative alpha spent",
                       "Number of subjects",
                       "Efficacy boundary (mean)",
                       "Efficacy boundary (p)",
                       "Information")
    }

    colnames(df) <- paste("Stage", seq_len(ncol(df)), sep=" ")
  } else {
    b <- s[, c("efficacyBounds", "efficacyMean", "efficacyP")]

    # format number of digits after decimal for each column
    j3 <- c(1,2)
    j4 <- 3

    b[j3] <- lapply(b[j3], formatC, format = "f", digits = 3)
    b[j4] <- lapply(b[j4], formatC, format = "f", digits = 4)

    df = t(b)

    if (!x$settings$normalApproximation) {
      rownames(df) = c("Efficacy boundary (t)",
                       "Efficacy boundary (mean)",
                       "Efficacy boundary (p)")
    } else {
      rownames(df) = c("Efficacy boundary (Z)",
                       "Efficacy boundary (mean)",
                       "Efficacy boundary (p)")
    }

    colnames(df) <- NA
  }

  print(df1, ..., na.print = "" , quote = FALSE )
  print(df, ..., na.print = "" , quote = FALSE )
  invisible(x)
}


#' @title Print Power and Sample Size Results for Paired Mean Difference
#' @description Prints the summary statistics from power calculation of
#' paired mean difference.
#'
#' @param x The designPairedMeanDiff object to print.
#' @param ... Ensures that all arguments starting from "..." are named.
#'
#' @return A tabular printout of the design elements.
#'
#' @keywords internal
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @export
print.designPairedMeanDiff <- function(x, ...) {
  a = x$overallResults
  s = x$byStageResults
  k = a$kMax

  if (k>1) {
    str1 = paste0("Group-sequential design with ", k, " stages")
  } else {
    str1 = "Fixed design"
  }

  str1 <- paste(str1, "for paired mean difference")

  str2 <- paste0("Paired difference under H0: ",
                 round(a$pairedDiffH0, 3), ", ",
                 "paired difference under H1: ",
                 round(a$pairedDiff, 3), ", ",
                 "standard deviation: ", round(a$stDev, 3))

  str3 <- paste0("Overall power: ",
                 round(a$overallReject, 3), ", ",
                 "overall alpha (1-sided): ",
                 round(a$alpha, 4))

  if (x$settings$typeBetaSpending != 'none' ||
      (k > 1 && any(x$byStageResults$futilityBounds[1:(k-1)] > -6))) {
    str3 <- paste0(str3, ", ",
                   "attained alpha: ", round(a$attainedAlpha, 4))
  }

  str4 <- paste0("Drift parameter: ", round(a$drift, 3), ", ",
                 "inflation factor: ", round(a$inflationFactor, 3))

  if (k>1) {
    str5 <- paste0("Maximum information: ",
                   round(a$information, 2), ", ",
                   "expected under H1: ",
                   round(a$expectedInformationH1, 2), ", ",
                   "expected under H0: ",
                   round(a$expectedInformationH0, 2))

    str6 <- paste0("Maximum # subjects: ",
                   round(a$numberOfSubjects, 1), ", ",
                   "expected under H1: ",
                   round(a$expectedNumberOfSubjectsH1, 1), ", ",
                   "expected under H0: ",
                   round(a$expectedNumberOfSubjectsH0, 1))

  } else {
    str5 <- paste0("Information: ",
                   round(a$information, 2))

    str6 <- paste0("Number of subjects: ",
                   round(a$numberOfSubjects, 1))
  }

  if (k > 1) {
    asf = tolower(x$settings$typeAlphaSpending)
    asfpar = round(x$settings$parameterAlphaSpending, 3)
    asfuser = round(x$settings$userAlphaSpending, 4)

    bsf = tolower(x$settings$typeBetaSpending)
    bsfpar = round(x$settings$parameterBetaSpending, 3)
    bsfuser = round(x$settings$userBetaSpending, 4)

    if (asf == "of") {
      str7 = paste0("Alpha spending: O'Brien-Fleming")
    } else if (asf == "p") {
      str7 = paste0("Alpha spending: Pocock")
    } else if (asf == "wt") {
      str7 = paste0("Alpha spending: Wang-Tsiatis(Delta = ", asfpar, ")")
    } else if (asf == "sfof") {
      str7 = paste0("Alpha spending: Lan-DeMets O'Brien-Fleming")
    } else if (asf == "sfp") {
      str7 = paste0("Alpha spending: Lan-DeMets Pocock")
    } else if (asf == "sfkd") {
      str7 = paste0("Alpha spending: KD(rho = ", asfpar, ")")
    } else if (asf == "sfhsd") {
      str7 = paste0("Alpha spending: HSD(gamma = ", asfpar, ")")
    } else if (asf == "user") {
      str7 = paste0("Alpha spending: User defined(",
                    paste(asfuser, collapse = ","), ")")
    } else {
      str7 = "Alpha spending: None"
    }

    if (bsf == "of") {
      str8 = paste0("beta spending: O'Brien-Fleming")
    } else if (bsf == "p") {
      str8 = paste0("beta spending: Pocock")
    } else if (bsf == "wt") {
      str8 = paste0("beta spending: Wang-Tsiatis(Delta = ", bsfpar, ")")
    } else if (bsf == "sfof") {
      str8 = paste0("beta spending: Lan-DeMets O'Brien-Fleming")
    } else if (bsf == "sfp") {
      str8 = paste0("beta spending: Lan-DeMets Pocock")
    } else if (bsf == "sfkd") {
      str8 = paste0("beta spending: KD(rho = ", bsfpar, ")")
    } else if (bsf == "sfhsd") {
      str8 = paste0("beta spending: HSD(gamma = ", bsfpar, ")")
    } else if (bsf == "user") {
      str8 = paste0("beta spending: User defined(",
                    paste(bsfuser, collapse = ","), ")")
    } else {
      str8 = "beta spending: None"
    }

    if (!any(is.na(x$settings$spendingTime)) &&
        !all(x$settings$spendingTime == s$informationRates)) {
      str9 = paste0("Spending time: ",
                    paste(round(x$settings$spendingTime, 3), collapse = ","))
      df1 = data.frame(x = rep("", 9))
      colnames(df1) = NULL
      rownames(df1) = c(str1, str2, str3, str4, str5, str6,
                        paste(str7, str8, sep = ", "), str9, "")
    } else {
      df1 = data.frame(x = rep("", 8))
      colnames(df1) = NULL
      rownames(df1) = c(str1, str2, str3, str4, str5, str6,
                        paste(str7, str8, sep = ", "), "")
    }
  } else {
    df1 = data.frame(x = rep("", 7))
    colnames(df1) = NULL
    rownames(df1) = c(str1, str2, str3, str4, str5, str6, "")
  }

  if (k>1) {
    b <- s[, c("informationRates", "efficacyBounds", "futilityBounds",
               "cumulativeRejection", "cumulativeFutility",
               "cumulativeAlphaSpent", "numberOfSubjects",
               "efficacyPairedDiff", "futilityPairedDiff",
               "efficacyP", "futilityP", "information",
               "cumulativeRejectionH0", "cumulativeFutilityH0")]

    # format number of digits after decimal for each column
    j1 <- 7
    j2 <- 12
    j3 <- c(1,2,3,4,8,9,14)
    j4 <- c(5,6,10,11,13)

    b[j1] <- lapply(b[j1], formatC, format = "f", digits = 1)
    b[j2] <- lapply(b[j2], formatC, format = "f", digits = 2)
    b[j3] <- lapply(b[j3], formatC, format = "f", digits = 3)
    b[j4] <- lapply(b[j4], formatC, format = "f", digits = 4)

    if (x$settings$typeBetaSpending != 'none' ||
        (k > 1 && any(x$byStageResults$futilityBounds[1:(k-1)] > -6))) {
      df = t(b)
      rownames(df) = c("Information rate",
                       "Efficacy boundary (Z)",
                       "Futility boundary (Z)",
                       "Cumulative rejection",
                       "Cumulative futility",
                       "Cumulative alpha spent",
                       "Number of subjects",
                       "Efficacy boundary (paired diff)",
                       "Futility boundary (paired diff)",
                       "Efficacy boundary (p)",
                       "Futility boundary (p)",
                       "Information",
                       "Cumulative rejection under H0",
                       "Cumulative futility under H0")

    } else {
      df = t(b[,c(1,2,4,6,7,8,10,12)])
      rownames(df) = c("Information rate",
                       "Efficacy boundary (Z)",
                       "Cumulative rejection",
                       "Cumulative alpha spent",
                       "Number of subjects",
                       "Efficacy boundary (paired diff)",
                       "Efficacy boundary (p)",
                       "Information")
    }

    colnames(df) <- paste("Stage", seq_len(ncol(df)), sep=" ")
  } else {
    b <- s[, c("efficacyBounds", "efficacyPairedDiff", "efficacyP")]

    # format number of digits after decimal for each column
    j3 <- c(1,2)
    j4 <- 3

    b[j3] <- lapply(b[j3], formatC, format = "f", digits = 3)
    b[j4] <- lapply(b[j4], formatC, format = "f", digits = 4)

    df = t(b)

    if (!x$settings$normalApproximation) {
      rownames(df) = c("Efficacy boundary (t)",
                       "Efficacy boundary (paired diff)",
                       "Efficacy boundary (p)")
    } else {
      rownames(df) = c("Efficacy boundary (Z)",
                       "Efficacy boundary (paired diff)",
                       "Efficacy boundary (p)")
    }

    colnames(df) <- NA
  }

  print(df1, ..., na.print = "" , quote = FALSE )
  print(df, ..., na.print = "" , quote = FALSE )
  invisible(x)
}


#' @title Print Power and Sample Size Results for Paired Mean Ratio
#' @description Prints the summary statistics from power calculation of
#' paired mean ratio.
#'
#' @param x The designPairedMeanRatio object to print.
#' @param ... Ensures that all arguments starting from "..." are named.
#'
#' @return A tabular printout of the design elements.
#'
#' @keywords internal
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @export
print.designPairedMeanRatio <- function(x, ...) {
  a = x$overallResults
  s = x$byStageResults
  k = a$kMax

  if (k>1) {
    str1 = paste0("Group-sequential design with ", k, " stages")
  } else {
    str1 = "Fixed design"
  }

  str1 <- paste(str1, "for paired mean ratio")

  str2 <- paste0("Paired ratio under H0: ",
                 round(a$pairedRatioH0, 3), ", ",
                 "paired ratio under H1: ",
                 round(a$pairedRatio, 3), ", ",
                 "coefficient of variation: ", round(a$CV, 3))

  str3 <- paste0("Overall power: ",
                 round(a$overallReject, 3), ", ",
                 "overall alpha: ",
                 round(a$alpha, 4))

  if (x$settings$typeBetaSpending != 'none' ||
      (k > 1 && any(x$byStageResults$futilityBounds[1:(k-1)] > -6))) {
    str3 <- paste0(str3, ", ",
                   "attained alpha: ", round(a$attainedAlpha, 4))
  }

  str4 <- paste0("Drift parameter: ", round(a$drift, 3), ", ",
                 "inflation factor: ", round(a$inflationFactor, 3))


  if (k>1) {
    str5 <- paste0("Maximum information: ",
                   round(a$information, 2), ", ",
                   "expected under H1: ",
                   round(a$expectedInformationH1, 2), ", ",
                   "expected under H0: ",
                   round(a$expectedInformationH0, 2))

    str6 <- paste0("Maximum # subjects: ",
                   round(a$numberOfSubjects, 1), ", ",
                   "expected under H1: ",
                   round(a$expectedNumberOfSubjectsH1, 1), ", ",
                   "expected under H0: ",
                   round(a$expectedNumberOfSubjectsH0, 1))

  } else {
    str5 <- paste0("Information: ",
                   round(a$information, 2))

    str6 <- paste0("Number of subjects: ",
                   round(a$numberOfSubjects, 1))
  }

  if (k > 1) {
    asf = tolower(x$settings$typeAlphaSpending)
    asfpar = round(x$settings$parameterAlphaSpending, 3)
    asfuser = round(x$settings$userAlphaSpending, 4)

    bsf = tolower(x$settings$typeBetaSpending)
    bsfpar = round(x$settings$parameterBetaSpending, 3)
    bsfuser = round(x$settings$userBetaSpending, 4)

    if (asf == "of") {
      str7 = paste0("Alpha spending: O'Brien-Fleming")
    } else if (asf == "p") {
      str7 = paste0("Alpha spending: Pocock")
    } else if (asf == "wt") {
      str7 = paste0("Alpha spending: Wang-Tsiatis(Delta = ", asfpar, ")")
    } else if (asf == "sfof") {
      str7 = paste0("Alpha spending: Lan-DeMets O'Brien-Fleming")
    } else if (asf == "sfp") {
      str7 = paste0("Alpha spending: Lan-DeMets Pocock")
    } else if (asf == "sfkd") {
      str7 = paste0("Alpha spending: KD(rho = ", asfpar, ")")
    } else if (asf == "sfhsd") {
      str7 = paste0("Alpha spending: HSD(gamma = ", asfpar, ")")
    } else if (asf == "user") {
      str7 = paste0("Alpha spending: User defined(",
                    paste(asfuser, collapse = ","), ")")
    } else {
      str7 = "Alpha spending: None"
    }

    if (bsf == "of") {
      str8 = paste0("beta spending: O'Brien-Fleming")
    } else if (bsf == "p") {
      str8 = paste0("beta spending: Pocock")
    } else if (bsf == "wt") {
      str8 = paste0("beta spending: Wang-Tsiatis(Delta = ", bsfpar, ")")
    } else if (bsf == "sfof") {
      str8 = paste0("beta spending: Lan-DeMets O'Brien-Fleming")
    } else if (bsf == "sfp") {
      str8 = paste0("beta spending: Lan-DeMets Pocock")
    } else if (bsf == "sfkd") {
      str8 = paste0("beta spending: KD(rho = ", bsfpar, ")")
    } else if (bsf == "sfhsd") {
      str8 = paste0("beta spending: HSD(gamma = ", bsfpar, ")")
    } else if (bsf == "user") {
      str8 = paste0("beta spending: User defined(",
                    paste(bsfuser, collapse = ","), ")")
    } else {
      str8 = "beta spending: None"
    }

    if (!any(is.na(x$settings$spendingTime)) &&
        !all(x$settings$spendingTime == s$informationRates)) {
      str9 = paste0("Spending time: ",
                    paste(round(x$settings$spendingTime, 3), collapse = ","))
      df1 = data.frame(x = rep("", 9))
      colnames(df1) = NULL
      rownames(df1) = c(str1, str2, str3, str4, str5, str6,
                        paste(str7, str8, sep = ", "), str9, "")
    } else {
      df1 = data.frame(x = rep("", 8))
      colnames(df1) = NULL
      rownames(df1) = c(str1, str2, str3, str4, str5, str6,
                        paste(str7, str8, sep = ", "), "")
    }
  } else {
    df1 = data.frame(x = rep("", 7))
    colnames(df1) = NULL
    rownames(df1) = c(str1, str2, str3, str4, str5, str6, "")
  }

  if (k>1) {
    b <- s[, c("informationRates", "efficacyBounds", "futilityBounds",
               "cumulativeRejection", "cumulativeFutility",
               "cumulativeAlphaSpent", "numberOfSubjects",
               "efficacyPairedRatio", "futilityPairedRatio",
               "efficacyP", "futilityP", "information",
               "cumulativeRejectionH0", "cumulativeFutilityH0")]

    # format number of digits after decimal for each column
    j1 <- 7
    j2 <- 12
    j3 <- c(1,2,3,4,8,9,14)
    j4 <- c(5,6,10,11,13)

    b[j1] <- lapply(b[j1], formatC, format = "f", digits = 1)
    b[j2] <- lapply(b[j2], formatC, format = "f", digits = 2)
    b[j3] <- lapply(b[j3], formatC, format = "f", digits = 3)
    b[j4] <- lapply(b[j4], formatC, format = "f", digits = 4)

    if (x$settings$typeBetaSpending != 'none' ||
        (k > 1 && any(x$byStageResults$futilityBounds[1:(k-1)] > -6))) {
      df = t(b)
      rownames(df) = c("Information rate",
                       "Efficacy boundary (Z)",
                       "Futility boundary (Z)",
                       "Cumulative rejection",
                       "Cumulative futility",
                       "Cumulative alpha spent",
                       "Number of subjects",
                       "Efficacy boundary (paired ratio)",
                       "Futility boundary (paired ratio)",
                       "Efficacy boundary (p)",
                       "Futility boundary (p)",
                       "Information",
                       "Cumulative rejection under H0",
                       "Cumulative futility under H0")

    } else {
      df = t(b[,c(1,2,4,6,7,8,10,12)])
      rownames(df) = c("Information rate",
                       "Efficacy boundary (Z)",
                       "Cumulative rejection",
                       "Cumulative alpha spent",
                       "Number of subjects",
                       "Efficacy boundary (paired ratio)",
                       "Efficacy boundary (p)",
                       "Information")
    }

    colnames(df) <- paste("Stage", seq_len(ncol(df)), sep=" ")
  } else {
    b <- s[, c("efficacyBounds", "efficacyPairedRatio", "efficacyP")]

    # format number of digits after decimal for each column
    j3 <- c(1,2)
    j4 <- 3

    b[j3] <- lapply(b[j3], formatC, format = "f", digits = 3)
    b[j4] <- lapply(b[j4], formatC, format = "f", digits = 4)

    df = t(b)

    if (!x$settings$normalApproximation) {
      rownames(df) = c("Efficacy boundary (t)",
                       "Efficacy boundary (paired ratio)",
                       "Efficacy boundary (p)")
    } else {
      rownames(df) = c("Efficacy boundary (Z)",
                       "Efficacy boundary (paired ratio)",
                       "Efficacy boundary (p)")
    }

    colnames(df) <- NA
  }

  print(df1, ..., na.print = "" , quote = FALSE )
  print(df, ..., na.print = "" , quote = FALSE )
  invisible(x)
}


#' @title Print Power and Sample Size Results for Two-Sample Mean Difference
#' @description Prints the summary statistics from power calculation of
#' two-sample mean difference.
#'
#' @param x The designMeanDiff object to print.
#' @param ... Ensures that all arguments starting from "..." are named.
#'
#' @return A tabular printout of the design elements.
#'
#' @keywords internal
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @export
print.designMeanDiff <- function(x, ...) {
  a = x$overallResults
  s = x$byStageResults
  k = a$kMax

  if (k>1) {
    str1 = paste0("Group-sequential design with ", k, " stages")
  } else {
    str1 = "Fixed design"
  }

  str1 <- paste(str1, "for two-sample mean difference")

  str2 <- paste0("Mean difference under H0: ",
                 round(a$meanDiffH0, 3), ", ",
                 "mean difference under H1: ",
                 round(a$meanDiff, 3), ", ",
                 "standard deviation: ", round(a$stDev, 3))

  str3 <- paste0("Overall power: ",
                 round(a$overallReject, 3), ", ",
                 "overall alpha (1-sided): ",
                 round(a$alpha, 4))

  if (x$settings$typeBetaSpending != 'none' ||
      (k > 1 && any(x$byStageResults$futilityBounds[1:(k-1)] > -6))) {
    str3 <- paste0(str3, ", ",
                   "attained alpha: ", round(a$attainedAlpha, 4))
  }

  str4 <- paste0("Drift parameter: ", round(a$drift, 3), ", ",
                 "inflation factor: ", round(a$inflationFactor, 3))

  if (k>1) {
    str5 <- paste0("Maximum information: ",
                   round(a$information, 2), ", ",
                   "expected under H1: ",
                   round(a$expectedInformationH1, 2), ", ",
                   "expected under H0: ",
                   round(a$expectedInformationH0, 2))

    str6 <- paste0("Maximum # subjects: ",
                   round(a$numberOfSubjects, 1), ", ",
                   "expected under H1: ",
                   round(a$expectedNumberOfSubjectsH1, 1), ", ",
                   "expected under H0: ",
                   round(a$expectedNumberOfSubjectsH0, 1))

  } else {
    str5 <- paste0("Information: ",
                   round(a$information, 2))

    str6 <- paste0("Number of subjects: ",
                   round(a$numberOfSubjects, 1))
  }

  str7 <- paste0("Allocation ratio: ",
                 round(x$settings$allocationRatioPlanned, 3))

  if (k > 1) {
    asf = tolower(x$settings$typeAlphaSpending)
    asfpar = round(x$settings$parameterAlphaSpending, 3)
    asfuser = round(x$settings$userAlphaSpending, 4)

    bsf = tolower(x$settings$typeBetaSpending)
    bsfpar = round(x$settings$parameterBetaSpending, 3)
    bsfuser = round(x$settings$userBetaSpending, 4)

    if (asf == "of") {
      str8 = paste0("Alpha spending: O'Brien-Fleming")
    } else if (asf == "p") {
      str8 = paste0("Alpha spending: Pocock")
    } else if (asf == "wt") {
      str8 = paste0("Alpha spending: Wang-Tsiatis(Delta = ", asfpar, ")")
    } else if (asf == "sfof") {
      str8 = paste0("Alpha spending: Lan-DeMets O'Brien-Fleming")
    } else if (asf == "sfp") {
      str8 = paste0("Alpha spending: Lan-DeMets Pocock")
    } else if (asf == "sfkd") {
      str8 = paste0("Alpha spending: KD(rho = ", asfpar, ")")
    } else if (asf == "sfhsd") {
      str8 = paste0("Alpha spending: HSD(gamma = ", asfpar, ")")
    } else if (asf == "user") {
      str8 = paste0("Alpha spending: User defined(",
                    paste(asfuser, collapse = ","), ")")
    } else {
      str8 = "Alpha spending: None"
    }

    if (bsf == "of") {
      str9 = paste0("beta spending: O'Brien-Fleming")
    } else if (bsf == "p") {
      str9 = paste0("beta spending: Pocock")
    } else if (bsf == "wt") {
      str9 = paste0("beta spending: Wang-Tsiatis(Delta = ", bsfpar, ")")
    } else if (bsf == "sfof") {
      str9 = paste0("beta spending: Lan-DeMets O'Brien-Fleming")
    } else if (bsf == "sfp") {
      str9 = paste0("beta spending: Lan-DeMets Pocock")
    } else if (bsf == "sfkd") {
      str9 = paste0("beta spending: KD(rho = ", bsfpar, ")")
    } else if (bsf == "sfhsd") {
      str9 = paste0("beta spending: HSD(gamma = ", bsfpar, ")")
    } else if (bsf == "user") {
      str9 = paste0("beta spending: User defined(",
                    paste(bsfuser, collapse = ","), ")")
    } else {
      str9 = "beta spending: None"
    }

    if (!any(is.na(x$settings$spendingTime)) &&
        !all(x$settings$spendingTime == s$informationRates)) {
      str10 = paste0("Spending time: ",
                     paste(round(x$settings$spendingTime, 3), collapse = ","))
      df1 = data.frame(x = rep("", 10))
      colnames(df1) = NULL
      rownames(df1) = c(str1, str2, str3, str4, str5, str6, str7,
                        paste(str8, str9, sep = ", "), str10, "")
    } else {
      df1 = data.frame(x = rep("", 9))
      colnames(df1) = NULL
      rownames(df1) = c(str1, str2, str3, str4, str5, str6, str7,
                        paste(str8, str9, sep = ", "), "")
    }
  } else {
    df1 = data.frame(x = rep("", 8))
    colnames(df1) = NULL
    rownames(df1) = c(str1, str2, str3, str4, str5, str6, str7, "")
  }

  if (k>1) {
    b <- s[, c("informationRates", "efficacyBounds", "futilityBounds",
               "cumulativeRejection", "cumulativeFutility",
               "cumulativeAlphaSpent", "numberOfSubjects",
               "efficacyMeanDiff", "futilityMeanDiff",
               "efficacyP", "futilityP", "information",
               "cumulativeRejectionH0", "cumulativeFutilityH0")]

    # format number of digits after decimal for each column
    j1 <- 7
    j2 <- 12
    j3 <- c(1,2,3,4,8,9,14)
    j4 <- c(5,6,10,11,13)

    b[j1] <- lapply(b[j1], formatC, format = "f", digits = 1)
    b[j2] <- lapply(b[j2], formatC, format = "f", digits = 2)
    b[j3] <- lapply(b[j3], formatC, format = "f", digits = 3)
    b[j4] <- lapply(b[j4], formatC, format = "f", digits = 4)

    if (x$settings$typeBetaSpending != 'none' ||
        (k > 1 && any(x$byStageResults$futilityBounds[1:(k-1)] > -6))) {
      df = t(b)
      rownames(df) = c("Information rate",
                       "Efficacy boundary (Z)",
                       "Futility boundary (Z)",
                       "Cumulative rejection",
                       "Cumulative futility",
                       "Cumulative alpha spent",
                       "Number of subjects",
                       "Efficacy boundary (mean diff)",
                       "Futility boundary (mean diff)",
                       "Efficacy boundary (p)",
                       "Futility boundary (p)",
                       "Information",
                       "Cumulative rejection under H0",
                       "Cumulative futility under H0")

    } else {
      df = t(b[,c(1,2,4,6,7,8,10,12)])
      rownames(df) = c("Information rate",
                       "Efficacy boundary (Z)",
                       "Cumulative rejection",
                       "Cumulative alpha spent",
                       "Number of subjects",
                       "Efficacy boundary (mean diff)",
                       "Efficacy boundary (p)",
                       "Information")
    }

    colnames(df) <- paste("Stage", seq_len(ncol(df)), sep=" ")
  } else {
    b <- s[, c("efficacyBounds", "efficacyMeanDiff", "efficacyP")]

    # format number of digits after decimal for each column
    j3 <- c(1,2)
    j4 <- 3

    b[j3] <- lapply(b[j3], formatC, format = "f", digits = 3)
    b[j4] <- lapply(b[j4], formatC, format = "f", digits = 4)

    df = t(b)

    if (x$settings$normalApproximation) {
      rownames(df) = c("Efficacy boundary (Z)",
                       "Efficacy boundary (mean diff)",
                       "Efficacy boundary (p)")
    } else {
      rownames(df) = c("Efficacy boundary (t)",
                       "Efficacy boundary (mean diff)",
                       "Efficacy boundary (p)")
    }

    colnames(df) <- NA
  }

  print(df1, ..., na.print = "" , quote = FALSE )
  print(df, ..., na.print = "" , quote = FALSE )
  invisible(x)
}


#' @title Print Power and Sample Size Results for Two-Sample Mean Ratio
#' @description Prints the summary statistics from power calculation of
#' two-sample mean ratio.
#'
#' @param x The designMeanRatio object to print.
#' @param ... Ensures that all arguments starting from "..." are named.
#'
#' @return A tabular printout of the design elements.
#'
#' @keywords internal
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @export
print.designMeanRatio <- function(x, ...) {
  a = x$overallResults
  s = x$byStageResults
  k = a$kMax

  if (k>1) {
    str1 = paste0("Group-sequential design with ", k, " stages")
  } else {
    str1 = "Fixed design"
  }

  str1 <- paste(str1, "for two-sample mean ratio")

  str2 <- paste0("Mean ratio under H0: ",
                 round(a$meanRatioH0, 3), ", ",
                 "mean ratio under H1: ",
                 round(a$meanRatio, 3), ", ",
                 "coefficient of variation: ", round(a$CV, 3))

  str3 <- paste0("Overall power: ",
                 round(a$overallReject, 3), ", ",
                 "overall alpha (1-sided): ",
                 round(a$alpha, 4))

  if (x$settings$typeBetaSpending != 'none' ||
      (k > 1 && any(x$byStageResults$futilityBounds[1:(k-1)] > -6))) {
    str3 <- paste0(str3, ", ",
                   "attained alpha: ", round(a$attainedAlpha, 4))
  }

  str4 <- paste0("Drift parameter: ", round(a$drift, 3), ", ",
                 "inflation factor: ", round(a$inflationFactor, 3))

  if (k>1) {
    str5 <- paste0("Maximum information: ",
                   round(a$information, 2), ", ",
                   "expected under H1: ",
                   round(a$expectedInformationH1, 2), ", ",
                   "expected under H0: ",
                   round(a$expectedInformationH0, 2))

    str6 <- paste0("Maximum # subjects: ",
                   round(a$numberOfSubjects, 1), ", ",
                   "expected under H1: ",
                   round(a$expectedNumberOfSubjectsH1, 1), ", ",
                   "expected under H0: ",
                   round(a$expectedNumberOfSubjectsH0, 1))

  } else {
    str5 <- paste0("Information: ",
                   round(a$information, 2))

    str6 <- paste0("Number of subjects: ",
                   round(a$numberOfSubjects, 1))
  }

  str7 <- paste0("Allocation ratio: ",
                 round(x$settings$allocationRatioPlanned, 3))

  if (k > 1) {
    asf = tolower(x$settings$typeAlphaSpending)
    asfpar = round(x$settings$parameterAlphaSpending, 3)
    asfuser = round(x$settings$userAlphaSpending, 4)

    bsf = tolower(x$settings$typeBetaSpending)
    bsfpar = round(x$settings$parameterBetaSpending, 3)
    bsfuser = round(x$settings$userBetaSpending, 4)

    if (asf == "of") {
      str8 = paste0("Alpha spending: O'Brien-Fleming")
    } else if (asf == "p") {
      str8 = paste0("Alpha spending: Pocock")
    } else if (asf == "wt") {
      str8 = paste0("Alpha spending: Wang-Tsiatis(Delta = ", asfpar, ")")
    } else if (asf == "sfof") {
      str8 = paste0("Alpha spending: Lan-DeMets O'Brien-Fleming")
    } else if (asf == "sfp") {
      str8 = paste0("Alpha spending: Lan-DeMets Pocock")
    } else if (asf == "sfkd") {
      str8 = paste0("Alpha spending: KD(rho = ", asfpar, ")")
    } else if (asf == "sfhsd") {
      str8 = paste0("Alpha spending: HSD(gamma = ", asfpar, ")")
    } else if (asf == "user") {
      str8 = paste0("Alpha spending: User defined(",
                    paste(asfuser, collapse = ","), ")")
    } else {
      str8 = "Alpha spending: None"
    }

    if (bsf == "of") {
      str9 = paste0("beta spending: O'Brien-Fleming")
    } else if (bsf == "p") {
      str9 = paste0("beta spending: Pocock")
    } else if (bsf == "wt") {
      str9 = paste0("beta spending: Wang-Tsiatis(Delta = ", bsfpar, ")")
    } else if (bsf == "sfof") {
      str9 = paste0("beta spending: Lan-DeMets O'Brien-Fleming")
    } else if (bsf == "sfp") {
      str9 = paste0("beta spending: Lan-DeMets Pocock")
    } else if (bsf == "sfkd") {
      str9 = paste0("beta spending: KD(rho = ", bsfpar, ")")
    } else if (bsf == "sfhsd") {
      str9 = paste0("beta spending: HSD(gamma = ", bsfpar, ")")
    } else if (bsf == "user") {
      str9 = paste0("beta spending: User defined(",
                    paste(bsfuser, collapse = ","), ")")
    } else {
      str9 = "beta spending: None"
    }

    if (!any(is.na(x$settings$spendingTime)) &&
        !all(x$settings$spendingTime == s$informationRates)) {
      str10 = paste0("Spending time: ",
                     paste(round(x$settings$spendingTime, 3), collapse = ","))
      df1 = data.frame(x = rep("", 10))
      colnames(df1) = NULL
      rownames(df1) = c(str1, str2, str3, str4, str5, str6, str7,
                        paste(str8, str9, sep = ", "), str10, "")
    } else {
      df1 = data.frame(x = rep("", 9))
      colnames(df1) = NULL
      rownames(df1) = c(str1, str2, str3, str4, str5, str6, str7,
                        paste(str8, str9, sep = ", "), "")
    }
  } else {
    df1 = data.frame(x = rep("", 8))
    colnames(df1) = NULL
    rownames(df1) = c(str1, str2, str3, str4, str5, str6, str7, "")
  }

  if (k>1) {
    b <- s[, c("informationRates", "efficacyBounds", "futilityBounds",
               "cumulativeRejection", "cumulativeFutility",
               "cumulativeAlphaSpent", "numberOfSubjects",
               "efficacyMeanRatio", "futilityMeanRatio",
               "efficacyP", "futilityP", "information",
               "cumulativeRejectionH0", "cumulativeFutilityH0")]

    # format number of digits after decimal for each column
    j1 <- 7
    j2 <- 12
    j3 <- c(1,2,3,4,8,9,14)
    j4 <- c(5,6,10,11,13)

    b[j1] <- lapply(b[j1], formatC, format = "f", digits = 1)
    b[j2] <- lapply(b[j2], formatC, format = "f", digits = 2)
    b[j3] <- lapply(b[j3], formatC, format = "f", digits = 3)
    b[j4] <- lapply(b[j4], formatC, format = "f", digits = 4)

    if (x$settings$typeBetaSpending != 'none' ||
        (k > 1 && any(x$byStageResults$futilityBounds[1:(k-1)] > -6))) {
      df = t(b)
      rownames(df) = c("Information rate",
                       "Efficacy boundary (Z)",
                       "Futility boundary (Z)",
                       "Cumulative rejection",
                       "Cumulative futility",
                       "Cumulative alpha spent",
                       "Number of subjects",
                       "Efficacy boundary (mean ratio)",
                       "Futility boundary (mean ratio)",
                       "Efficacy boundary (p)",
                       "Futility boundary (p)",
                       "Information",
                       "Cumulative rejection under H0",
                       "Cumulative futility under H0")

    } else {
      df = t(b[,c(1,2,4,6,7,8,10,12)])
      rownames(df) = c("Information rate",
                       "Efficacy boundary (Z)",
                       "Cumulative rejection",
                       "Cumulative alpha spent",
                       "Number of subjects",
                       "Efficacy boundary (mean ratio)",
                       "Efficacy boundary (p)",
                       "Information")
    }

    colnames(df) <- paste("Stage", seq_len(ncol(df)), sep=" ")
  } else {
    b <- s[, c("efficacyBounds", "efficacyMeanRatio", "efficacyP")]

    # format number of digits after decimal for each column
    j3 <- c(1,2)
    j4 <- 3

    b[j3] <- lapply(b[j3], formatC, format = "f", digits = 3)
    b[j4] <- lapply(b[j4], formatC, format = "f", digits = 4)

    df = t(b)

    if (x$settings$normalApproximation) {
      rownames(df) = c("Efficacy boundary (Z)",
                       "Efficacy boundary (mean ratio)",
                       "Efficacy boundary (p)")
    } else {
      rownames(df) = c("Efficacy boundary (t)",
                       "Efficacy boundary (mean ratio)",
                       "Efficacy boundary (p)")
    }

    colnames(df) <- NA
  }

  print(df1, ..., na.print = "" , quote = FALSE )
  print(df, ..., na.print = "" , quote = FALSE )
  invisible(x)
}


#' @title Print Power and Sample Size Results for Mean Difference in 2x2
#' Crossover
#' @description Prints the summary statistics from power calculation of
#' mean difference in 2x2 crossover.
#'
#' @param x The designMeanDiffXO object to print.
#' @param ... Ensures that all arguments starting from "..." are named.
#'
#' @return A tabular printout of the design elements.
#'
#' @keywords internal
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @export
print.designMeanDiffXO <- function(x, ...) {
  a = x$overallResults
  s = x$byStageResults
  k = a$kMax

  if (k>1) {
    str1 = paste0("Group-sequential design with ", k, " stages")
  } else {
    str1 = "Fixed design"
  }

  str1 <- paste(str1, "for mean difference in 2x2 crossover")

  str2 <- paste0("Mean difference under H0: ",
                 round(a$meanDiffH0, 3), ", ",
                 "mean difference under H1: ",
                 round(a$meanDiff, 3), ", ",
                 "standard deviation: ",
                 round(a$stDev, 3))

  str3 <- paste0("Overall power: ",
                 round(a$overallReject, 3), ", ",
                 "overall alpha (1-sided): ",
                 round(a$alpha, 4))

  if (x$settings$typeBetaSpending != 'none' ||
      (k > 1 && any(x$byStageResults$futilityBounds[1:(k-1)] > -6))) {
    str3 <- paste0(str3, ", ",
                   "attained alpha: ", round(a$attainedAlpha, 4))
  }

  str4 <- paste0("Drift parameter: ", round(a$drift, 3), ", ",
                 "inflation factor: ", round(a$inflationFactor, 3))

  if (k>1) {
    str5 <- paste0("Maximum information: ",
                   round(a$information, 2), ", ",
                   "expected under H1: ",
                   round(a$expectedInformationH1, 2), ", ",
                   "expected under H0: ",
                   round(a$expectedInformationH0, 2))

    str6 <- paste0("Maximum # subjects: ",
                   round(a$numberOfSubjects, 1), ", ",
                   "expected under H1: ",
                   round(a$expectedNumberOfSubjectsH1, 1), ", ",
                   "expected under H0: ",
                   round(a$expectedNumberOfSubjectsH0, 1))

  } else {
    str5 <- paste0("Information: ",
                   round(a$information, 2))

    str6 <- paste0("Number of subjects: ",
                   round(a$numberOfSubjects, 1))
  }

  str7 <- paste0("Sequence allocation ratio: ",
                 round(x$settings$allocationRatioPlanned, 3))

  if (k > 1) {
    asf = tolower(x$settings$typeAlphaSpending)
    asfpar = round(x$settings$parameterAlphaSpending, 3)
    asfuser = round(x$settings$userAlphaSpending, 4)

    bsf = tolower(x$settings$typeBetaSpending)
    bsfpar = round(x$settings$parameterBetaSpending, 3)
    bsfuser = round(x$settings$userBetaSpending, 4)

    if (asf == "of") {
      str8 = paste0("Alpha spending: O'Brien-Fleming")
    } else if (asf == "p") {
      str8 = paste0("Alpha spending: Pocock")
    } else if (asf == "wt") {
      str8 = paste0("Alpha spending: Wang-Tsiatis(Delta = ", asfpar, ")")
    } else if (asf == "sfof") {
      str8 = paste0("Alpha spending: Lan-DeMets O'Brien-Fleming")
    } else if (asf == "sfp") {
      str8 = paste0("Alpha spending: Lan-DeMets Pocock")
    } else if (asf == "sfkd") {
      str8 = paste0("Alpha spending: KD(rho = ", asfpar, ")")
    } else if (asf == "sfhsd") {
      str8 = paste0("Alpha spending: HSD(gamma = ", asfpar, ")")
    } else if (asf == "user") {
      str8 = paste0("Alpha spending: User defined(",
                    paste(asfuser, collapse = ","), ")")
    } else {
      str8 = "Alpha spending: None"
    }

    if (bsf == "of") {
      str9 = paste0("beta spending: O'Brien-Fleming")
    } else if (bsf == "p") {
      str9 = paste0("beta spending: Pocock")
    } else if (bsf == "wt") {
      str9 = paste0("beta spending: Wang-Tsiatis(Delta = ", bsfpar, ")")
    } else if (bsf == "sfof") {
      str9 = paste0("beta spending: Lan-DeMets O'Brien-Fleming")
    } else if (bsf == "sfp") {
      str9 = paste0("beta spending: Lan-DeMets Pocock")
    } else if (bsf == "sfkd") {
      str9 = paste0("beta spending: KD(rho = ", bsfpar, ")")
    } else if (bsf == "sfhsd") {
      str9 = paste0("beta spending: HSD(gamma = ", bsfpar, ")")
    } else if (bsf == "user") {
      str9 = paste0("beta spending: User defined(",
                    paste(bsfuser, collapse = ","), ")")
    } else {
      str9 = "beta spending: None"
    }

    if (!any(is.na(x$settings$spendingTime)) &&
        !all(x$settings$spendingTime == s$informationRates)) {
      str10 = paste0("Spending time: ",
                     paste(round(x$settings$spendingTime, 3), collapse = ","))
      df1 = data.frame(x = rep("", 10))
      colnames(df1) = NULL
      rownames(df1) = c(str1, str2, str3, str4, str5, str6, str7,
                        paste(str8, str9, sep = ", "), str10, "")
    } else {
      df1 = data.frame(x = rep("", 9))
      colnames(df1) = NULL
      rownames(df1) = c(str1, str2, str3, str4, str5, str6, str7,
                        paste(str8, str9, sep = ", "), "")
    }
  } else {
    df1 = data.frame(x = rep("", 8))
    colnames(df1) = NULL
    rownames(df1) = c(str1, str2, str3, str4, str5, str6, str7, "")
  }

  if (k>1) {
    b <- s[, c("informationRates", "efficacyBounds", "futilityBounds",
               "cumulativeRejection", "cumulativeFutility",
               "cumulativeAlphaSpent", "numberOfSubjects",
               "efficacyMeanDiff", "futilityMeanDiff",
               "efficacyP", "futilityP", "information",
               "cumulativeRejectionH0", "cumulativeFutilityH0")]

    # format number of digits after decimal for each column
    j1 <- 7
    j2 <- 12
    j3 <- c(1,2,3,4,8,9,14)
    j4 <- c(5,6,10,11,13)

    b[j1] <- lapply(b[j1], formatC, format = "f", digits = 1)
    b[j2] <- lapply(b[j2], formatC, format = "f", digits = 2)
    b[j3] <- lapply(b[j3], formatC, format = "f", digits = 3)
    b[j4] <- lapply(b[j4], formatC, format = "f", digits = 4)

    if (x$settings$typeBetaSpending != 'none' ||
        (k > 1 && any(x$byStageResults$futilityBounds[1:(k-1)] > -6))) {
      df = t(b)
      rownames(df) = c("Information rate",
                       "Efficacy boundary (Z)",
                       "Futility boundary (Z)",
                       "Cumulative rejection",
                       "Cumulative futility",
                       "Cumulative alpha spent",
                       "Number of subjects",
                       "Efficacy boundary (mean diff)",
                       "Futility boundary (mean diff)",
                       "Efficacy boundary (p)",
                       "Futility boundary (p)",
                       "Information",
                       "Cumulative rejection under H0",
                       "Cumulative futility under H0")

    } else {
      df = t(b[,c(1,2,4,6,7,8,10,12)])
      rownames(df) = c("Information rate",
                       "Efficacy boundary (Z)",
                       "Cumulative rejection",
                       "Cumulative alpha spent",
                       "Number of subjects",
                       "Efficacy boundary (mean diff)",
                       "Efficacy boundary (p)",
                       "Information")
    }

    colnames(df) <- paste("Stage", seq_len(ncol(df)), sep=" ")
  } else {
    b <- s[, c("efficacyBounds", "efficacyMeanDiff", "efficacyP")]

    # format number of digits after decimal for each column
    j3 <- c(1,2)
    j4 <- 3

    b[j3] <- lapply(b[j3], formatC, format = "f", digits = 3)
    b[j4] <- lapply(b[j4], formatC, format = "f", digits = 4)

    df = t(b)

    if (x$settings$normalApproximation) {
      rownames(df) = c("Efficacy boundary (Z)",
                       "Efficacy boundary (mean diff)",
                       "Efficacy boundary (p)")
    } else {
      rownames(df) = c("Efficacy boundary (t)",
                       "Efficacy boundary (mean diff)",
                       "Efficacy boundary (p)")
    }

    colnames(df) <- NA
  }

  print(df1, ..., na.print = "" , quote = FALSE )
  print(df, ..., na.print = "" , quote = FALSE )
  invisible(x)
}


#' @title Print Power and Sample Size Results for Mean Ratio in 2x2 Crossover
#' @description Prints the summary statistics from power calculation of
#' mean ratio in 2x2 crossover.
#'
#' @param x The designMeanRatioXO object to print.
#' @param ... Ensures that all arguments starting from "..." are named.
#'
#' @return A tabular printout of the design elements.
#'
#' @keywords internal
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @export
print.designMeanRatioXO <- function(x, ...) {
  a = x$overallResults
  s = x$byStageResults
  k = a$kMax

  if (k>1) {
    str1 = paste0("Group-sequential design with ", k, " stages")
  } else {
    str1 = "Fixed design"
  }

  str1 <- paste(str1, "for mean ratio in 2x2 crossover")

  str2 <- paste0("Mean ratio under H0: ",
                 round(a$meanRatioH0, 3), ", ",
                 "mean ratio under H1: ",
                 round(a$meanRatio, 3), ", ",
                 "coefficient of variation: ", round(a$CV, 3))

  str3 <- paste0("Overall power: ",
                 round(a$overallReject, 3), ", ",
                 "overall alpha (1-sided): ",
                 round(a$alpha, 4))

  if (x$settings$typeBetaSpending != 'none' ||
      (k > 1 && any(x$byStageResults$futilityBounds[1:(k-1)] > -6))) {
    str3 <- paste0(str3, ", ",
                   "attained alpha: ", round(a$attainedAlpha, 4))
  }

  str4 <- paste0("Drift parameter: ", round(a$drift, 3), ", ",
                 "inflation factor: ", round(a$inflationFactor, 3))

  if (k>1) {
    str5 <- paste0("Maximum information: ",
                   round(a$information, 2), ", ",
                   "expected under H1: ",
                   round(a$expectedInformationH1, 2), ", ",
                   "expected under H0: ",
                   round(a$expectedInformationH0, 2))

    str6 <- paste0("Maximum # subjects: ",
                   round(a$numberOfSubjects, 1), ", ",
                   "expected under H1: ",
                   round(a$expectedNumberOfSubjectsH1, 1), ", ",
                   "expected under H0: ",
                   round(a$expectedNumberOfSubjectsH0, 1))

  } else {
    str5 <- paste0("Information: ",
                   round(a$information, 2))

    str6 <- paste0("Number of subjects: ",
                   round(a$numberOfSubjects, 1))
  }

  str7 <- paste0("Sequence allocation ratio: ",
                 round(x$settings$allocationRatioPlanned, 3))

  if (k > 1) {
    asf = tolower(x$settings$typeAlphaSpending)
    asfpar = round(x$settings$parameterAlphaSpending, 3)
    asfuser = round(x$settings$userAlphaSpending, 4)

    bsf = tolower(x$settings$typeBetaSpending)
    bsfpar = round(x$settings$parameterBetaSpending, 3)
    bsfuser = round(x$settings$userBetaSpending, 4)

    if (asf == "of") {
      str8 = paste0("Alpha spending: O'Brien-Fleming")
    } else if (asf == "p") {
      str8 = paste0("Alpha spending: Pocock")
    } else if (asf == "wt") {
      str8 = paste0("Alpha spending: Wang-Tsiatis(Delta = ", asfpar, ")")
    } else if (asf == "sfof") {
      str8 = paste0("Alpha spending: Lan-DeMets O'Brien-Fleming")
    } else if (asf == "sfp") {
      str8 = paste0("Alpha spending: Lan-DeMets Pocock")
    } else if (asf == "sfkd") {
      str8 = paste0("Alpha spending: KD(rho = ", asfpar, ")")
    } else if (asf == "sfhsd") {
      str8 = paste0("Alpha spending: HSD(gamma = ", asfpar, ")")
    } else if (asf == "user") {
      str8 = paste0("Alpha spending: User defined(",
                    paste(asfuser, collapse = ","), ")")
    } else {
      str8 = "Alpha spending: None"
    }

    if (bsf == "of") {
      str9 = paste0("beta spending: O'Brien-Fleming")
    } else if (bsf == "p") {
      str9 = paste0("beta spending: Pocock")
    } else if (bsf == "wt") {
      str9 = paste0("beta spending: Wang-Tsiatis(Delta = ", bsfpar, ")")
    } else if (bsf == "sfof") {
      str9 = paste0("beta spending: Lan-DeMets O'Brien-Fleming")
    } else if (bsf == "sfp") {
      str9 = paste0("beta spending: Lan-DeMets Pocock")
    } else if (bsf == "sfkd") {
      str9 = paste0("beta spending: KD(rho = ", bsfpar, ")")
    } else if (bsf == "sfhsd") {
      str9 = paste0("beta spending: HSD(gamma = ", bsfpar, ")")
    } else if (bsf == "user") {
      str9 = paste0("beta spending: User defined(",
                    paste(bsfuser, collapse = ","), ")")
    } else {
      str9 = "beta spending: None"
    }

    if (!any(is.na(x$settings$spendingTime)) &&
        !all(x$settings$spendingTime == s$informationRates)) {
      str10 = paste0("Spending time: ",
                     paste(round(x$settings$spendingTime, 3), collapse = ","))
      df1 = data.frame(x = rep("", 10))
      colnames(df1) = NULL
      rownames(df1) = c(str1, str2, str3, str4, str5, str6, str7,
                        paste(str8, str9, sep = ", "), str10, "")
    } else {
      df1 = data.frame(x = rep("", 9))
      colnames(df1) = NULL
      rownames(df1) = c(str1, str2, str3, str4, str5, str6, str7,
                        paste(str8, str9, sep = ", "), "")
    }
  } else {
    df1 = data.frame(x = rep("", 8))
    colnames(df1) = NULL
    rownames(df1) = c(str1, str2, str3, str4, str5, str6, str7, "")
  }

  if (k>1) {
    b <- s[, c("informationRates", "efficacyBounds", "futilityBounds",
               "cumulativeRejection", "cumulativeFutility",
               "cumulativeAlphaSpent", "numberOfSubjects",
               "efficacyMeanRatio", "futilityMeanRatio",
               "efficacyP", "futilityP", "information",
               "cumulativeRejectionH0", "cumulativeFutilityH0")]

    # format number of digits after decimal for each column
    j1 <- 7
    j2 <- 12
    j3 <- c(1,2,3,4,8,9,14)
    j4 <- c(5,6,10,11,13)

    b[j1] <- lapply(b[j1], formatC, format = "f", digits = 1)
    b[j2] <- lapply(b[j2], formatC, format = "f", digits = 2)
    b[j3] <- lapply(b[j3], formatC, format = "f", digits = 3)
    b[j4] <- lapply(b[j4], formatC, format = "f", digits = 4)

    if (x$settings$typeBetaSpending != 'none' ||
        (k > 1 && any(x$byStageResults$futilityBounds[1:(k-1)] > -6))) {
      df = t(b)
      rownames(df) = c("Information rate",
                       "Efficacy boundary (Z)",
                       "Futility boundary (Z)",
                       "Cumulative rejection",
                       "Cumulative futility",
                       "Cumulative alpha spent",
                       "Number of subjects",
                       "Efficacy boundary (mean ratio)",
                       "Futility boundary (mean ratio)",
                       "Efficacy boundary (p)",
                       "Futility boundary (p)",
                       "Information",
                       "Cumulative rejection under H0",
                       "Cumulative futility under H0")

    } else {
      df = t(b[,c(1,2,4,6,7,8,10,12)])
      rownames(df) = c("Information rate",
                       "Efficacy boundary (Z)",
                       "Cumulative rejection",
                       "Cumulative alpha spent",
                       "Number of subjects",
                       "Efficacy boundary (mean ratio)",
                       "Efficacy boundary (p)",
                       "Information")
    }

    colnames(df) <- paste("Stage", seq_len(ncol(df)), sep=" ")
  } else {
    b <- s[, c("efficacyBounds", "efficacyMeanRatio", "efficacyP")]

    # format number of digits after decimal for each column
    j3 <- c(1,2)
    j4 <- 3

    b[j3] <- lapply(b[j3], formatC, format = "f", digits = 3)
    b[j4] <- lapply(b[j4], formatC, format = "f", digits = 4)

    df = t(b)

    if (!x$settings$normalApproximation) {
      rownames(df) = c("Efficacy boundary (t)",
                       "Efficacy boundary (mean ratio)",
                       "Efficacy boundary (p)")
    } else {
      rownames(df) = c("Efficacy boundary (Z)",
                       "Efficacy boundary (mean ratio)",
                       "Efficacy boundary (p)")
    }

    colnames(df) <- NA
  }

  print(df1, ..., na.print = "" , quote = FALSE )
  print(df, ..., na.print = "" , quote = FALSE )
  invisible(x)
}


#' @title Print Power and Sample Size Results for Equivalence in Paired
#' Mean Difference
#' @description Prints the summary statistics from power calculation of
#' equivalence in paired mean difference.
#'
#' @param x The designPairedMeanDiffEquiv object to print.
#' @param ... Ensures that all arguments starting from "..." are named.
#'
#' @return A tabular printout of the design elements.
#'
#' @keywords internal
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @export
print.designPairedMeanDiffEquiv <- function(x, ...) {
  a = x$overallResults
  s = x$byStageResults
  k = a$kMax

  if (k>1) {
    str1 = paste0("Group-sequential design with ", k, " stages")
  } else {
    str1 = "Fixed design"
  }

  str1 <- paste(str1, "for equivalence in paired mean difference")

  str2 <- paste0("Lower limit for paired difference: ",
                 round(a$pairedDiffLower, 3), ", ",
                 "upper limit for paired difference: ",
                 round(a$pairedDiffUpper, 3))

  str3 <- paste0("Paired difference under H1: ",
                 round(a$pairedDiff, 3), ", ",
                 "standard deviation for paired difference: ",
                 round(a$stDev, 3))

  str4 <- paste0("Overall power: ",
                 round(a$overallReject, 3), ", ",
                 "overall alpha: ",
                 round(a$alpha, 4), ", ",
                 "attained alpha: ",
                 round(a$attainedAlpha, 4))

  if (k>1) {
    str5 <- paste0("Maximum information: ",
                   round(a$information, 2), ", ",
                   "expected under H1: ",
                   round(a$expectedInformationH1, 2), ", ",
                   "expected under H0: ",
                   round(a$expectedInformationH0, 2))

    str6 <- paste0("Maximum # subjects: ",
                   round(a$numberOfSubjects, 1), ", ",
                   "expected under H1: ",
                   round(a$expectedNumberOfSubjectsH1, 1), ", ",
                   "expected under H0: ",
                   round(a$expectedNumberOfSubjectsH0, 1))

  } else {
    str5 <- paste0("Information: ",
                   round(a$information, 2))

    str6 <- paste0("Number of subjects: ",
                   round(a$numberOfSubjects, 1))
  }

  if (k > 1) {
    asf = tolower(x$settings$typeAlphaSpending)
    asfpar = round(x$settings$parameterAlphaSpending, 3)
    asfuser = round(x$settings$userAlphaSpending, 4)

    if (asf == "of") {
      str7 = paste0("Alpha spending: O'Brien-Fleming")
    } else if (asf == "p") {
      str7 = paste0("Alpha spending: Pocock")
    } else if (asf == "wt") {
      str7 = paste0("Alpha spending: Wang-Tsiatis(Delta = ", asfpar, ")")
    } else if (asf == "sfof") {
      str7 = paste0("Alpha spending: Lan-DeMets O'Brien-Fleming")
    } else if (asf == "sfp") {
      str7 = paste0("Alpha spending: Lan-DeMets Pocock")
    } else if (asf == "sfkd") {
      str7 = paste0("Alpha spending: KD(rho = ", asfpar, ")")
    } else if (asf == "sfhsd") {
      str7 = paste0("Alpha spending: HSD(gamma = ", asfpar, ")")
    } else if (asf == "user") {
      str7 = paste0("Alpha spending: User defined(",
                    paste(asfuser, collapse = ","), ")")
    } else {
      str7 = "Alpha spending: None"
    }

    if (!any(is.na(x$settings$spendingTime)) &&
        !all(x$settings$spendingTime == s$informationRates)) {
      str8 = paste0("Spending time: ",
                    paste(round(x$settings$spendingTime, 3), collapse = ","))
      df1 = data.frame(x = rep("", 9))
      colnames(df1) = NULL
      rownames(df1) = c(str1, str2, str3, str4, str5, str6, str7, str8, "")
    } else {
      df1 = data.frame(x = rep("", 8))
      colnames(df1) = NULL
      rownames(df1) = c(str1, str2, str3, str4, str5, str6, str7, "")
    }
  } else {
    df1 = data.frame(x = rep("", 7))
    colnames(df1) = NULL
    rownames(df1) = c(str1, str2, str3, str4, str5, str6, "")
  }

  if (k>1) {
    b <- s[, c("informationRates", "efficacyBounds",
               "cumulativeRejection", "cumulativeAlphaSpent",
               "cumulativeAttainedAlpha", "numberOfSubjects",
               "efficacyPairedDiffLower", "efficacyPairedDiffUpper",
               "efficacyP", "information")]

    # format number of digits after decimal for each column
    j1 <- 6
    j2 <- 10
    j3 <- c(1,2,3,7,8)
    j4 <- c(4,5,9)

    b[j1] <- lapply(b[j1], formatC, format = "f", digits = 1)
    b[j2] <- lapply(b[j2], formatC, format = "f", digits = 2)
    b[j3] <- lapply(b[j3], formatC, format = "f", digits = 3)
    b[j4] <- lapply(b[j4], formatC, format = "f", digits = 4)

    df = t(b)
    rownames(df) = c("Information rate",
                     "Boundary for each 1-sided test (Z)",
                     "Cumulative rejection",
                     "Cumulative alpha for each 1-sided test",
                     "Cumulative alpha attained under H0",
                     "Number of subjects",
                     "Boundary for lower limit (paired diff)",
                     "Boundary for upper limit (paired diff)",
                     "Boundary for each 1-sided test (p)",
                     "Information")

    colnames(df) <- paste("Stage", seq_len(ncol(df)), sep=" ")
  } else {
    b <- s[, c("efficacyBounds", "efficacyPairedDiffLower",
               "efficacyPairedDiffUpper", "efficacyP")]

    # format number of digits after decimal for each column
    j3 <- c(1,2,3)
    j4 <- 4

    b[j3] <- lapply(b[j3], formatC, format = "f", digits = 3)
    b[j4] <- lapply(b[j4], formatC, format = "f", digits = 4)

    df = t(b)

    if (x$settings$normalApproximation) {
      rownames(df) = c("Boundary for each 1-sided test (Z)",
                       "Boundary for lower limit (paired diff)",
                       "Boundary for upper limit (paired diff)",
                       "Boundary for each 1-sided test (p)")
    } else {
      rownames(df) = c("Boundary for each 1-sided test (t)",
                       "Boundary for lower limit (paired diff)",
                       "Boundary for upper limit (paired diff)",
                       "Boundary for each 1-sided test (p)")
    }

    colnames(df) <- NA
  }

  print(df1, ..., na.print = "" , quote = FALSE )
  print(df, ..., na.print = "" , quote = FALSE )
  invisible(x)
}


#' @title Print Power and Sample Size Results for Equivalence in Paired
#' Mean Ratio
#' @description Prints the summary statistics from power calculation of
#' equivalence in paired mean ratio.
#'
#' @param x The designPairedMeanRatioEquiv object to print.
#' @param ... Ensures that all arguments starting from "..." are named.
#'
#' @return A tabular printout of the design elements.
#'
#' @keywords internal
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @export
print.designPairedMeanRatioEquiv <- function(x, ...) {
  a = x$overallResults
  s = x$byStageResults
  k = a$kMax

  if (k>1) {
    str1 = paste0("Group-sequential design with ", k, " stages")
  } else {
    str1 = "Fixed design"
  }

  str1 <- paste(str1, "for equivalence in paired mean ratio")

  str2 <- paste0("Lower limit for paired ratio: ",
                 round(a$pairedRatioLower, 3), ", ",
                 "upper limit for paired ratio: ",
                 round(a$pairedRatioUpper, 3))

  str3 <- paste0("Paired ratio under H1: ",
                 round(a$pairedRatio, 3), ", ",
                 "coefficient of variation for paired ratio: ",
                 round(a$CV, 3))

  str4 <- paste0("Overall power: ",
                 round(a$overallReject, 3), ", ",
                 "overall alpha: ",
                 round(a$alpha, 4), ", ",
                 "attained alpha: ",
                 round(a$attainedAlpha, 4))

  if (k>1) {
    str5 <- paste0("Maximum information: ",
                   round(a$information, 2), ", ",
                   "expected under H1: ",
                   round(a$expectedInformationH1, 2), ", ",
                   "expected under H0: ",
                   round(a$expectedInformationH0, 2))

    str6 <- paste0("Maximum # subjects: ",
                   round(a$numberOfSubjects, 1), ", ",
                   "expected under H1: ",
                   round(a$expectedNumberOfSubjectsH1, 1), ", ",
                   "expected under H0: ",
                   round(a$expectedNumberOfSubjectsH0, 1))

  } else {
    str5 <- paste0("Information: ",
                   round(a$information, 2))

    str6 <- paste0("Number of subjects: ",
                   round(a$numberOfSubjects, 1))
  }

  if (k > 1) {
    asf = tolower(x$settings$typeAlphaSpending)
    asfpar = round(x$settings$parameterAlphaSpending, 3)
    asfuser = round(x$settings$userAlphaSpending, 4)

    if (asf == "of") {
      str7 = paste0("Alpha spending: O'Brien-Fleming")
    } else if (asf == "p") {
      str7 = paste0("Alpha spending: Pocock")
    } else if (asf == "wt") {
      str7 = paste0("Alpha spending: Wang-Tsiatis(Delta = ", asfpar, ")")
    } else if (asf == "sfof") {
      str7 = paste0("Alpha spending: Lan-DeMets O'Brien-Fleming")
    } else if (asf == "sfp") {
      str7 = paste0("Alpha spending: Lan-DeMets Pocock")
    } else if (asf == "sfkd") {
      str7 = paste0("Alpha spending: KD(rho = ", asfpar, ")")
    } else if (asf == "sfhsd") {
      str7 = paste0("Alpha spending: HSD(gamma = ", asfpar, ")")
    } else if (asf == "user") {
      str7 = paste0("Alpha spending: User defined(",
                    paste(asfuser, collapse = ","), ")")
    } else {
      str7 = "Alpha spending: None"
    }

    if (!any(is.na(x$settings$spendingTime)) &&
        !all(x$settings$spendingTime == s$informationRates)) {
      str8 = paste0("Spending time: ",
                    paste(round(x$settings$spendingTime, 3), collapse = ","))
      df1 = data.frame(x = rep("", 9))
      colnames(df1) = NULL
      rownames(df1) = c(str1, str2, str3, str4, str5, str6, str7, str8, "")
    } else {
      df1 = data.frame(x = rep("", 8))
      colnames(df1) = NULL
      rownames(df1) = c(str1, str2, str3, str4, str5, str6, str7, "")
    }
  } else {
    df1 = data.frame(x = rep("", 7))
    colnames(df1) = NULL
    rownames(df1) = c(str1, str2, str3, str4, str5, str6, "")
  }

  if (k>1) {
    b <- s[, c("informationRates", "efficacyBounds",
               "cumulativeRejection", "cumulativeAlphaSpent",
               "cumulativeAttainedAlpha", "numberOfSubjects",
               "efficacyPairedRatioLower", "efficacyPairedRatioUpper",
               "efficacyP", "information")]

    # format number of digits after decimal for each column
    j1 <- 6
    j2 <- 10
    j3 <- c(1,2,3,7,8)
    j4 <- c(4,5,9)

    b[j1] <- lapply(b[j1], formatC, format = "f", digits = 1)
    b[j2] <- lapply(b[j2], formatC, format = "f", digits = 2)
    b[j3] <- lapply(b[j3], formatC, format = "f", digits = 3)
    b[j4] <- lapply(b[j4], formatC, format = "f", digits = 4)

    df = t(b)
    rownames(df) = c("Information rate",
                     "Boundary for each 1-sided test (Z)",
                     "Cumulative rejection",
                     "Cumulative alpha for each 1-sided test",
                     "Cumulative alpha attained under H0",
                     "Number of subjects",
                     "Boundary for lower limit (paired ratio)",
                     "Boundary for upper limit (paired ratio)",
                     "Boundary for each 1-sided test (p)",
                     "Information")

    colnames(df) <- paste("Stage", seq_len(ncol(df)), sep=" ")
  } else {
    b <- s[, c("efficacyBounds", "efficacyPairedRatioLower",
               "efficacyPairedRatioUpper", "efficacyP")]

    # format number of digits after decimal for each column
    j3 <- c(1,2,3)
    j4 <- 4

    b[j3] <- lapply(b[j3], formatC, format = "f", digits = 3)
    b[j4] <- lapply(b[j4], formatC, format = "f", digits = 4)

    df = t(b)

    if (x$settings$normalApproximation) {
      rownames(df) = c("Boundary for each 1-sided test (Z)",
                       "Boundary for lower limit (paired ratio)",
                       "Boundary for upper limit (paired ratio)",
                       "Boundary for each 1-sided test (p)")
    } else {
      rownames(df) = c("Boundary for each 1-sided test (t)",
                       "Boundary for lower limit (paired ratio)",
                       "Boundary for upper limit (paired ratio)",
                       "Boundary for each 1-sided test (p)")
    }

    colnames(df) <- NA
  }

  print(df1, ..., na.print = "" , quote = FALSE )
  print(df, ..., na.print = "" , quote = FALSE )
  invisible(x)
}


#' @title Print Power and Sample Size Results for Equivalence in
#' Two-Sample Mean Difference
#' @description Prints the summary statistics from power calculation of
#' equivalence in two-sample mean difference.
#'
#' @param x The designMeanDiffEquiv object to print.
#' @param ... Ensures that all arguments starting from "..." are named.
#'
#' @return A tabular printout of the design elements.
#'
#' @keywords internal
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @export
print.designMeanDiffEquiv <- function(x, ...) {
  a = x$overallResults
  s = x$byStageResults
  k = a$kMax

  if (k>1) {
    str1 = paste0("Group-sequential design with ", k, " stages")
  } else {
    str1 = "Fixed design"
  }

  str1 <- paste(str1, "for equivalence in two-sample mean difference")

  str2 <- paste0("Lower limit for mean difference: ",
                 round(a$meanDiffLower, 3), ", ",
                 "upper limit for mean difference: ",
                 round(a$meanDiffUpper, 3))

  str3 <- paste0("Mean difference under H1: ",
                 round(a$meanDiff, 3), ", ",
                 "standard deviation: ",
                 round(a$stDev, 3))

  str4 <- paste0("Overall power: ",
                 round(a$overallReject, 3), ", ",
                 "overall alpha: ",
                 round(a$alpha, 4), ", ",
                 "attained alpha: ",
                 round(a$attainedAlpha, 4))

  if (k>1) {
    str5 <- paste0("Maximum information: ",
                   round(a$information, 2), ", ",
                   "expected under H1: ",
                   round(a$expectedInformationH1, 2), ", ",
                   "expected under H0: ",
                   round(a$expectedInformationH0, 2))

    str6 <- paste0("Maximum # subjects: ",
                   round(a$numberOfSubjects, 1), ", ",
                   "expected under H1: ",
                   round(a$expectedNumberOfSubjectsH1, 1), ", ",
                   "expected under H0: ",
                   round(a$expectedNumberOfSubjectsH0, 1))

  } else {
    str5 <- paste0("Information: ",
                   round(a$information, 2))

    str6 <- paste0("Number of subjects: ",
                   round(a$numberOfSubjects, 1))
  }

  str7 <- paste0("Allocation ratio: ",
                 round(x$settings$allocationRatioPlanned, 3))

  if (k > 1) {
    asf = tolower(x$settings$typeAlphaSpending)
    asfpar = round(x$settings$parameterAlphaSpending, 3)
    asfuser = round(x$settings$userAlphaSpending, 4)

    if (asf == "of") {
      str8 = paste0("Alpha spending: O'Brien-Fleming")
    } else if (asf == "p") {
      str8 = paste0("Alpha spending: Pocock")
    } else if (asf == "wt") {
      str8 = paste0("Alpha spending: Wang-Tsiatis(Delta = ", asfpar, ")")
    } else if (asf == "sfof") {
      str8 = paste0("Alpha spending: Lan-DeMets O'Brien-Fleming")
    } else if (asf == "sfp") {
      str8 = paste0("Alpha spending: Lan-DeMets Pocock")
    } else if (asf == "sfkd") {
      str8 = paste0("Alpha spending: KD(rho = ", asfpar, ")")
    } else if (asf == "sfhsd") {
      str8 = paste0("Alpha spending: HSD(gamma = ", asfpar, ")")
    } else if (asf == "user") {
      str8 = paste0("Alpha spending: User defined(",
                    paste(asfuser, collapse = ","), ")")
    } else {
      str8 = "Alpha spending: None"
    }

    if (!any(is.na(x$settings$spendingTime)) &&
        !all(x$settings$spendingTime == s$informationRates)) {
      str9 = paste0("Spending time: ",
                    paste(round(x$settings$spendingTime, 3), collapse = ","))
      df1 = data.frame(x = rep("", 10))
      colnames(df1) = NULL
      rownames(df1) = c(str1, str2, str3, str4, str5, str6, str7,
                        str8, str9, "")
    } else {
      df1 = data.frame(x = rep("", 9))
      colnames(df1) = NULL
      rownames(df1) = c(str1, str2, str3, str4, str5, str6, str7, str8, "")
    }
  } else {
    df1 = data.frame(x = rep("", 8))
    colnames(df1) = NULL
    rownames(df1) = c(str1, str2, str3, str4, str5, str6, str7, "")
  }

  if (k>1) {
    b <- s[, c("informationRates", "efficacyBounds",
               "cumulativeRejection", "cumulativeAlphaSpent",
               "cumulativeAttainedAlpha", "numberOfSubjects",
               "efficacyMeanDiffLower", "efficacyMeanDiffUpper",
               "efficacyP", "information")]

    # format number of digits after decimal for each column
    j1 <- 6
    j2 <- 10
    j3 <- c(1,2,3,7,8)
    j4 <- c(4,5,9)

    b[j1] <- lapply(b[j1], formatC, format = "f", digits = 1)
    b[j2] <- lapply(b[j2], formatC, format = "f", digits = 2)
    b[j3] <- lapply(b[j3], formatC, format = "f", digits = 3)
    b[j4] <- lapply(b[j4], formatC, format = "f", digits = 4)

    df = t(b)
    rownames(df) = c("Information rate",
                     "Boundary for each 1-sided test (Z)",
                     "Cumulative rejection",
                     "Cumulative alpha for each 1-sided test",
                     "Cumulative alpha attained under H0",
                     "Number of subjects",
                     "Boundary for lower limit (mean diff)",
                     "Boundary for upper limit (mean diff)",
                     "Boundary for each 1-sided test (p)",
                     "Information")

    colnames(df) <- paste("Stage", seq_len(ncol(df)), sep=" ")
  } else {
    b <- s[, c("efficacyBounds", "efficacyMeanDiffLower",
               "efficacyMeanDiffUpper", "efficacyP")]

    # format number of digits after decimal for each column
    j3 <- c(1,2,3)
    j4 <- 4

    b[j3] <- lapply(b[j3], formatC, format = "f", digits = 3)
    b[j4] <- lapply(b[j4], formatC, format = "f", digits = 4)

    df = t(b)

    if (x$settings$normalApproximation) {
      rownames(df) = c("Boundary for each 1-sided test (Z)",
                       "Boundary for lower limit (mean diff)",
                       "Boundary for upper limit (mean diff)",
                       "Boundary for each 1-sided test (p)")
    } else {
      rownames(df) = c("Boundary for each 1-sided test (t)",
                       "Boundary for lower limit (mean diff)",
                       "Boundary for upper limit (mean diff)",
                       "Boundary for each 1-sided test (p)")
    }

    colnames(df) <- NA
  }

  print(df1, ..., na.print = "" , quote = FALSE )
  print(df, ..., na.print = "" , quote = FALSE )
  invisible(x)
}


#' @title Print Power and Sample Size Results for Equivalence in
#' Two-Sample Mean Ratio
#' @description Prints the summary statistics from power calculation of
#' equivalence in two-sample mean ratio.
#'
#' @param x The designMeanRatioEquiv object to print.
#' @param ... Ensures that all arguments starting from "..." are named.
#'
#' @return A tabular printout of the design elements.
#'
#' @keywords internal
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @export
print.designMeanRatioEquiv <- function(x, ...) {
  a = x$overallResults
  s = x$byStageResults
  k = a$kMax

  if (k>1) {
    str1 = paste0("Group-sequential design with ", k, " stages")
  } else {
    str1 = "Fixed design"
  }

  str1 <- paste(str1, "for equivalence in two-sample mean ratio")

  str2 <- paste0("Lower limit for mean ratio: ",
                 round(a$meanRatioLower, 3), ", ",
                 "upper limit for mean ratio: ",
                 round(a$meanRatioUpper, 3))

  str3 <- paste0("Mean ratio under H1: ",
                 round(a$meanRatio, 3), ", ",
                 "coefficient of variation: ",
                 round(a$CV, 3))

  str4 <- paste0("Overall power: ",
                 round(a$overallReject, 3), ", ",
                 "overall alpha (1-sided): ",
                 round(a$alpha, 4), ", ",
                 "attained alpha: ",
                 round(a$attainedAlpha, 4))

  if (k>1) {
    str5 <- paste0("Maximum information: ",
                   round(a$information, 2), ", ",
                   "expected under H1: ",
                   round(a$expectedInformationH1, 2), ", ",
                   "expected under H0: ",
                   round(a$expectedInformationH0, 2))

    str6 <- paste0("Maximum # subjects: ",
                   round(a$numberOfSubjects, 1), ", ",
                   "expected under H1: ",
                   round(a$expectedNumberOfSubjectsH1, 1), ", ",
                   "expected under H0: ",
                   round(a$expectedNumberOfSubjectsH0, 1))

  } else {
    str5 <- paste0("Information: ",
                   round(a$information, 2))

    str6 <- paste0("Number of subjects: ",
                   round(a$numberOfSubjects, 1))
  }

  str7 <- paste0("Allocation ratio: ",
                 round(x$settings$allocationRatioPlanned, 3))

  if (k > 1) {
    asf = tolower(x$settings$typeAlphaSpending)
    asfpar = round(x$settings$parameterAlphaSpending, 3)
    asfuser = round(x$settings$userAlphaSpending, 4)

    if (asf == "of") {
      str8 = paste0("Alpha spending: O'Brien-Fleming")
    } else if (asf == "p") {
      str8 = paste0("Alpha spending: Pocock")
    } else if (asf == "wt") {
      str8 = paste0("Alpha spending: Wang-Tsiatis(Delta = ", asfpar, ")")
    } else if (asf == "sfof") {
      str8 = paste0("Alpha spending: Lan-DeMets O'Brien-Fleming")
    } else if (asf == "sfp") {
      str8 = paste0("Alpha spending: Lan-DeMets Pocock")
    } else if (asf == "sfkd") {
      str8 = paste0("Alpha spending: KD(rho = ", asfpar, ")")
    } else if (asf == "sfhsd") {
      str8 = paste0("Alpha spending: HSD(gamma = ", asfpar, ")")
    } else if (asf == "user") {
      str8 = paste0("Alpha spending: User defined(",
                    paste(asfuser, collapse = ","), ")")
    } else {
      str8 = "Alpha spending: None"
    }

    if (!any(is.na(x$settings$spendingTime)) &&
        !all(x$settings$spendingTime == s$informationRates)) {
      str9 = paste0("Spending time: ",
                    paste(round(x$settings$spendingTime, 3), collapse = ","))
      df1 = data.frame(x = rep("", 10))
      colnames(df1) = NULL
      rownames(df1) = c(str1, str2, str3, str4, str5, str6, str7,
                        str8, str9, "")
    } else {
      df1 = data.frame(x = rep("", 9))
      colnames(df1) = NULL
      rownames(df1) = c(str1, str2, str3, str4, str5, str6, str7, str8, "")
    }
  } else {
    df1 = data.frame(x = rep("", 8))
    colnames(df1) = NULL
    rownames(df1) = c(str1, str2, str3, str4, str5, str6, str7, "")
  }

  if (k>1) {
    b <- s[, c("informationRates", "efficacyBounds",
               "cumulativeRejection", "cumulativeAlphaSpent",
               "cumulativeAttainedAlpha", "numberOfSubjects",
               "efficacyMeanRatioLower", "efficacyMeanRatioUpper",
               "efficacyP", "information")]

    # format number of digits after decimal for each column
    j1 <- 6
    j2 <- 10
    j3 <- c(1,2,3,7,8)
    j4 <- c(4,5,9)

    b[j1] <- lapply(b[j1], formatC, format = "f", digits = 1)
    b[j2] <- lapply(b[j2], formatC, format = "f", digits = 2)
    b[j3] <- lapply(b[j3], formatC, format = "f", digits = 3)
    b[j4] <- lapply(b[j4], formatC, format = "f", digits = 4)

    df = t(b)
    rownames(df) = c("Information rate",
                     "Boundary for each 1-sided test (Z)",
                     "Cumulative rejection",
                     "Cumulative alpha for each 1-sided test",
                     "Cumulative alpha spent",
                     "Number of subjects",
                     "Boundary for lower limit (mean ratio)",
                     "Boundary for upper limit (mean ratio)",
                     "Boundary for each 1-sided test (p)",
                     "Information")

    colnames(df) <- paste("Stage", seq_len(ncol(df)), sep=" ")
  } else {
    b <- s[, c("efficacyBounds", "efficacyMeanRatioLower",
               "efficacyMeanRatioUpper", "efficacyP")]

    # format number of digits after decimal for each column
    j3 <- c(1,2,3)
    j4 <- 4

    b[j3] <- lapply(b[j3], formatC, format = "f", digits = 3)
    b[j4] <- lapply(b[j4], formatC, format = "f", digits = 4)

    df = t(b)

    if (x$settings$normalApproximation) {
      rownames(df) = c("Boundary for each 1-sided test (Z)",
                       "Boundary for lower limit (mean ratio)",
                       "Boundary for upper limit (mean ratio)",
                       "Boundary for each 1-sided test (p)")
    } else {
      rownames(df) = c("Boundary for each 1-sided test (t)",
                       "Boundary for lower limit (mean ratio)",
                       "Boundary for upper limit (mean ratio)",
                       "Boundary for each 1-sided test (p)")
    }

    colnames(df) <- NA
  }

  print(df1, ..., na.print = "" , quote = FALSE )
  print(df, ..., na.print = "" , quote = FALSE )
  invisible(x)
}


#' @title Print Power and Sample Size Results for Equivalence in
#' Mean Difference in 2x2 Crossover
#' @description Prints the summary statistics from power calculation of
#' equivalence in mean difference in 2x2 crossover.
#'
#' @param x The designMeanDiffXOEquiv object to print.
#' @param ... Ensures that all arguments starting from "..." are named.
#'
#' @return A tabular printout of the design elements.
#'
#' @keywords internal
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @export
print.designMeanDiffXOEquiv <- function(x, ...) {
  a = x$overallResults
  s = x$byStageResults
  k = a$kMax

  if (k>1) {
    str1 = paste0("Group-sequential design with ", k, " stages")
  } else {
    str1 = "Fixed design"
  }

  str1 <- paste(str1, "for equivalence in mean difference in 2x2 crossover")

  str2 <- paste0("Lower limit for mean difference: ",
                 round(a$meanDiffLower, 3), ", ",
                 "upper limit for mean difference: ",
                 round(a$meanDiffUpper, 3))

  str3 <- paste0("Mean difference under H1: ",
                 round(a$meanDiff, 3), ", ",
                 "standard deviation: ",
                 round(a$stDev, 3))

  str4 <- paste0("Overall power: ",
                 round(a$overallReject, 3), ", ",
                 "overall alpha (1-sided): ",
                 round(a$alpha, 4), ", ",
                 "attained alpha: ",
                 round(a$attainedAlpha, 4))

  if (k>1) {
    str5 <- paste0("Maximum information: ",
                   round(a$information, 2), ", ",
                   "expected under H1: ",
                   round(a$expectedInformationH1, 2), ", ",
                   "expected under H0: ",
                   round(a$expectedInformationH0, 2))

    str6 <- paste0("Maximum # subjects: ",
                   round(a$numberOfSubjects, 1), ", ",
                   "expected under H1: ",
                   round(a$expectedNumberOfSubjectsH1, 1), ", ",
                   "expected under H0: ",
                   round(a$expectedNumberOfSubjectsH0, 1))

  } else {
    str5 <- paste0("Information: ",
                   round(a$information, 2))

    str6 <- paste0("Number of subjects: ",
                   round(a$numberOfSubjects, 1))
  }

  str7 <- paste0("Sequence allocation ratio: ",
                 round(x$settings$allocationRatioPlanned, 3))

  if (k > 1) {
    asf = tolower(x$settings$typeAlphaSpending)
    asfpar = round(x$settings$parameterAlphaSpending, 3)
    asfuser = round(x$settings$userAlphaSpending, 4)

    if (asf == "of") {
      str8 = paste0("Alpha spending: O'Brien-Fleming")
    } else if (asf == "p") {
      str8 = paste0("Alpha spending: Pocock")
    } else if (asf == "wt") {
      str8 = paste0("Alpha spending: Wang-Tsiatis(Delta = ", asfpar, ")")
    } else if (asf == "sfof") {
      str8 = paste0("Alpha spending: Lan-DeMets O'Brien-Fleming")
    } else if (asf == "sfp") {
      str8 = paste0("Alpha spending: Lan-DeMets Pocock")
    } else if (asf == "sfkd") {
      str8 = paste0("Alpha spending: KD(rho = ", asfpar, ")")
    } else if (asf == "sfhsd") {
      str8 = paste0("Alpha spending: HSD(gamma = ", asfpar, ")")
    } else if (asf == "user") {
      str8 = paste0("Alpha spending: User defined(",
                    paste(asfuser, collapse = ","), ")")
    } else {
      str8 = "Alpha spending: None"
    }

    if (!any(is.na(x$settings$spendingTime)) &&
        !all(x$settings$spendingTime == s$informationRates)) {
      str9 = paste0("Spending time: ",
                    paste(round(x$settings$spendingTime, 3), collapse = ","))
      df1 = data.frame(x = rep("", 10))
      colnames(df1) = NULL
      rownames(df1) = c(str1, str2, str3, str4, str5, str6, str7,
                        str8, str9, "")
    } else {
      df1 = data.frame(x = rep("", 9))
      colnames(df1) = NULL
      rownames(df1) = c(str1, str2, str3, str4, str5, str6, str7,
                        str8, "")
    }
  } else {
    df1 = data.frame(x = rep("", 8))
    colnames(df1) = NULL
    rownames(df1) = c(str1, str2, str3, str4, str5, str6, str7, "")
  }

  if (k>1) {
    b <- s[, c("informationRates", "efficacyBounds",
               "cumulativeRejection", "cumulativeAlphaSpent",
               "cumulativeAttainedAlpha", "numberOfSubjects",
               "efficacyMeanDiffLower", "efficacyMeanDiffUpper",
               "efficacyP", "information")]

    # format number of digits after decimal for each column
    j1 <- 6
    j2 <- 10
    j3 <- c(1,2,3,7,8)
    j4 <- c(4,5,9)

    b[j1] <- lapply(b[j1], formatC, format = "f", digits = 1)
    b[j2] <- lapply(b[j2], formatC, format = "f", digits = 2)
    b[j3] <- lapply(b[j3], formatC, format = "f", digits = 3)
    b[j4] <- lapply(b[j4], formatC, format = "f", digits = 4)

    df = t(b)
    rownames(df) = c("Information rate",
                     "Boundary for each 1-sided test (Z)",
                     "Cumulative rejection",
                     "Cumulative alpha for each 1-sided test",
                     "Cumulative alpha spent",
                     "Number of subjects",
                     "Boundary for lower limit (mean diff)",
                     "Boundary for upper limit (mean diff)",
                     "Boundary for each 1-sided test (p)",
                     "Information")

    colnames(df) <- paste("Stage", seq_len(ncol(df)), sep=" ")
  } else {
    b <- s[, c("efficacyBounds", "efficacyMeanDiffLower",
               "efficacyMeanDiffUpper", "efficacyP")]

    # format number of digits after decimal for each column
    j3 <- c(1,2,3)
    j4 <- 4

    b[j3] <- lapply(b[j3], formatC, format = "f", digits = 3)
    b[j4] <- lapply(b[j4], formatC, format = "f", digits = 4)

    df = t(b)

    if (x$settings$normalApproximation) {
      rownames(df) = c("Boundary for each 1-sided test (Z)",
                       "Boundary for lower limit (mean diff)",
                       "Boundary for upper limit (mean diff)",
                       "Boundary for each 1-sided test (p)")
    } else {
      rownames(df) = c("Boundary for each 1-sided test (t)",
                       "Boundary for lower limit (mean diff)",
                       "Boundary for upper limit (mean diff)",
                       "Boundary for each 1-sided test (p)")
    }

    colnames(df) <- NA
  }

  print(df1, ..., na.print = "" , quote = FALSE )
  print(df, ..., na.print = "" , quote = FALSE )
  invisible(x)
}


#' @title Print Power and Sample Size Results for Equivalence in
#' Mean Ratio in 2x2 Crossover
#' @description Prints the summary statistics from power calculation of
#' equivalence in mean ratio in 2x2 crossover.
#'
#' @param x The designMeanRatioXOEquiv object to print.
#' @param ... Ensures that all arguments starting from "..." are named.
#'
#' @return A tabular printout of the design elements.
#'
#' @keywords internal
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @export
print.designMeanRatioXOEquiv <- function(x, ...) {
  a = x$overallResults
  s = x$byStageResults
  k = a$kMax

  if (k>1) {
    str1 = paste0("Group-sequential design with ", k, " stages")
  } else {
    str1 = "Fixed design"
  }

  str1 <- paste(str1, "for equivalence in mean ratio in 2x2 crossover")

  str2 <- paste0("Lower limit for mean ratio: ",
                 round(a$meanRatioLower, 3), ", ",
                 "upper limit for mean ratio: ",
                 round(a$meanRatioUpper, 3))

  str3 <- paste0("Mean ratio under H1: ",
                 round(a$meanRatio, 3), ", ",
                 "coefficient of variation: ",
                 round(a$CV, 3))

  str4 <- paste0("Overall power: ",
                 round(a$overallReject, 3), ", ",
                 "overall alpha (1-sided): ",
                 round(a$alpha, 4), ", ",
                 "attained alpha: ",
                 round(a$attainedAlpha, 4))

  if (k>1) {
    str5 <- paste0("Maximum information: ",
                   round(a$information, 2), ", ",
                   "expected under H1: ",
                   round(a$expectedInformationH1, 2), ", ",
                   "expected under H0: ",
                   round(a$expectedInformationH0, 2))

    str6 <- paste0("Maximum # subjects: ",
                   round(a$numberOfSubjects, 1), ", ",
                   "expected under H1: ",
                   round(a$expectedNumberOfSubjectsH1, 1), ", ",
                   "expected under H0: ",
                   round(a$expectedNumberOfSubjectsH0, 1))

  } else {
    str5 <- paste0("Information: ",
                   round(a$information, 2))

    str6 <- paste0("Number of subjects: ",
                   round(a$numberOfSubjects, 1))
  }

  str7 <- paste0("Sequence allocation ratio: ",
                 round(x$settings$allocationRatioPlanned, 3))

  if (k > 1) {
    asf = tolower(x$settings$typeAlphaSpending)
    asfpar = round(x$settings$parameterAlphaSpending, 3)
    asfuser = round(x$settings$userAlphaSpending, 4)

    if (asf == "of") {
      str8 = paste0("Alpha spending: O'Brien-Fleming")
    } else if (asf == "p") {
      str8 = paste0("Alpha spending: Pocock")
    } else if (asf == "wt") {
      str8 = paste0("Alpha spending: Wang-Tsiatis(Delta = ", asfpar, ")")
    } else if (asf == "sfof") {
      str8 = paste0("Alpha spending: Lan-DeMets O'Brien-Fleming")
    } else if (asf == "sfp") {
      str8 = paste0("Alpha spending: Lan-DeMets Pocock")
    } else if (asf == "sfkd") {
      str8 = paste0("Alpha spending: KD(rho = ", asfpar, ")")
    } else if (asf == "sfhsd") {
      str8 = paste0("Alpha spending: HSD(gamma = ", asfpar, ")")
    } else if (asf == "user") {
      str8 = paste0("Alpha spending: User defined(",
                    paste(asfuser, collapse = ","), ")")
    } else {
      str8 = "Alpha spending: None"
    }

    if (!any(is.na(x$settings$spendingTime)) &&
        !all(x$settings$spendingTime == s$informationRates)) {
      str9 = paste0("Spending time: ",
                    paste(round(x$settings$spendingTime, 3), collapse = ","))
      df1 = data.frame(x = rep("", 10))
      colnames(df1) = NULL
      rownames(df1) = c(str1, str2, str3, str4, str5, str6, str7,
                        str8, str9, "")
    } else {
      df1 = data.frame(x = rep("", 9))
      colnames(df1) = NULL
      rownames(df1) = c(str1, str2, str3, str4, str5, str6, str7, str8, "")
    }
  } else {
    df1 = data.frame(x = rep("", 8))
    colnames(df1) = NULL
    rownames(df1) = c(str1, str2, str3, str4, str5, str6, str7, "")
  }

  if (k>1) {
    b <- s[, c("informationRates", "efficacyBounds",
               "cumulativeRejection", "cumulativeAlphaSpent",
               "cumulativeAttainedAlpha", "numberOfSubjects",
               "efficacyMeanRatioLower", "efficacyMeanRatioUpper",
               "efficacyP", "information")]

    # format number of digits after decimal for each column
    j1 <- 6
    j2 <- 10
    j3 <- c(1,2,3,7,8)
    j4 <- c(4,5,9)

    b[j1] <- lapply(b[j1], formatC, format = "f", digits = 1)
    b[j2] <- lapply(b[j2], formatC, format = "f", digits = 2)
    b[j3] <- lapply(b[j3], formatC, format = "f", digits = 3)
    b[j4] <- lapply(b[j4], formatC, format = "f", digits = 4)

    df = t(b)
    rownames(df) = c("Information rate",
                     "Boundary for each 1-sided test (Z)",
                     "Cumulative rejection",
                     "Cumulative alpha for each 1-sided test",
                     "Cumulative alpha spent",
                     "Number of subjects",
                     "Boundary for lower limit (mean ratio)",
                     "Boundary for upper limit (mean ratio)",
                     "Boundary for each 1-sided test (p)",
                     "Information")

    colnames(df) <- paste("Stage", seq_len(ncol(df)), sep=" ")
  } else {
    b <- s[, c("efficacyBounds", "efficacyMeanRatioLower",
               "efficacyMeanRatioUpper", "efficacyP")]

    # format number of digits after decimal for each column
    j3 <- c(1,2,3)
    j4 <- 4

    b[j3] <- lapply(b[j3], formatC, format = "f", digits = 3)
    b[j4] <- lapply(b[j4], formatC, format = "f", digits = 4)

    df = t(b)

    if (x$settings$normalApproximation) {
      rownames(df) = c("Boundary for each 1-sided test (Z)",
                       "Boundary for lower limit (mean ratio)",
                       "Boundary for upper limit (mean ratio)",
                       "Boundary for each 1-sided test (p)")
    } else {
      rownames(df) = c("Boundary for each 1-sided test (t)",
                       "Boundary for lower limit (mean ratio)",
                       "Boundary for upper limit (mean ratio)",
                       "Boundary for each 1-sided test (p)")
    }

    colnames(df) <- NA
  }

  print(df1, ..., na.print = "" , quote = FALSE )
  print(df, ..., na.print = "" , quote = FALSE )
  invisible(x)
}


#' @title Print Power and Sample Size Results for Two-Sample Wilcoxon Test
#' @description Prints the summary statistics from power calculation of
#' two-sample Wilcoxon test.
#'
#' @param x The designWilcoxon object to print.
#' @param ... Ensures that all arguments starting from "..." are named.
#'
#' @return A tabular printout of the design elements.
#'
#' @keywords internal
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @export
print.designWilcoxon <- function(x, ...) {
  a = x$overallResults
  s = x$byStageResults
  k = a$kMax

  if (k>1) {
    str1 = paste0("Group-sequential design with ", k, " stages")
  } else {
    str1 = "Fixed design"
  }

  str1 <- paste(str1, "for two-sample Wilcoxon test")

  str2 <- paste0("Probability of observations in treatment ",
                 "larger than those in control: ",
                 round(a$pLarger, 3))

  str3 <- paste0("Overall power: ",
                 round(a$overallReject, 3), ", ",
                 "overall alpha (1-sided): ",
                 round(a$alpha, 4))

  if (x$settings$typeBetaSpending != 'none' ||
      (k > 1 && any(x$byStageResults$futilityBounds[1:(k-1)] > -6))) {
    str3 <- paste0(str3, ", ",
                   "attained alpha: ", round(a$attainedAlpha, 4))
  }

  str4 <- paste0("Drift parameter: ", round(a$drift, 3), ", ",
                 "inflation factor: ", round(a$inflationFactor, 3))

  if (k>1) {
    str5 <- paste0("Maximum information: ",
                   round(a$information, 2), ", ",
                   "expected under H1: ",
                   round(a$expectedInformationH1, 2), ", ",
                   "expected under H0: ",
                   round(a$expectedInformationH0, 2))

    str6 <- paste0("Maximum # subjects: ",
                   round(a$numberOfSubjects, 1), ", ",
                   "expected under H1: ",
                   round(a$expectedNumberOfSubjectsH1, 1), ", ",
                   "expected under H0: ",
                   round(a$expectedNumberOfSubjectsH0, 1))

  } else {
    str5 <- paste0("Information: ",
                   round(a$information, 2))

    str6 <- paste0("Number of subjects: ",
                   round(a$numberOfSubjects, 1))
  }

  str7 <- paste0("Allocation ratio: ",
                 round(x$settings$allocationRatioPlanned, 3))

  if (k > 1) {
    asf = tolower(x$settings$typeAlphaSpending)
    asfpar = round(x$settings$parameterAlphaSpending, 3)
    asfuser = round(x$settings$userAlphaSpending, 4)

    bsf = tolower(x$settings$typeBetaSpending)
    bsfpar = round(x$settings$parameterBetaSpending, 3)
    bsfuser = round(x$settings$userBetaSpending, 4)

    if (asf == "of") {
      str8 = paste0("Alpha spending: O'Brien-Fleming")
    } else if (asf == "p") {
      str8 = paste0("Alpha spending: Pocock")
    } else if (asf == "wt") {
      str8 = paste0("Alpha spending: Wang-Tsiatis(Delta = ", asfpar, ")")
    } else if (asf == "sfof") {
      str8 = paste0("Alpha spending: Lan-DeMets O'Brien-Fleming")
    } else if (asf == "sfp") {
      str8 = paste0("Alpha spending: Lan-DeMets Pocock")
    } else if (asf == "sfkd") {
      str8 = paste0("Alpha spending: KD(rho = ", asfpar, ")")
    } else if (asf == "sfhsd") {
      str8 = paste0("Alpha spending: HSD(gamma = ", asfpar, ")")
    } else if (asf == "user") {
      str8 = paste0("Alpha spending: User defined(",
                    paste(asfuser, collapse = ","), ")")
    } else {
      str8 = "Alpha spending: None"
    }

    if (bsf == "of") {
      str9 = paste0("beta spending: O'Brien-Fleming")
    } else if (bsf == "p") {
      str9 = paste0("beta spending: Pocock")
    } else if (bsf == "wt") {
      str9 = paste0("beta spending: Wang-Tsiatis(Delta = ", bsfpar, ")")
    } else if (bsf == "sfof") {
      str9 = paste0("beta spending: Lan-DeMets O'Brien-Fleming")
    } else if (bsf == "sfp") {
      str9 = paste0("beta spending: Lan-DeMets Pocock")
    } else if (bsf == "sfkd") {
      str9 = paste0("beta spending: KD(rho = ", bsfpar, ")")
    } else if (bsf == "sfhsd") {
      str9 = paste0("beta spending: HSD(gamma = ", bsfpar, ")")
    } else if (bsf == "user") {
      str9 = paste0("beta spending: User defined(",
                    paste(bsfuser, collapse = ","), ")")
    } else {
      str9 = "beta spending: None"
    }

    if (!any(is.na(x$settings$spendingTime)) &&
        !all(x$settings$spendingTime == s$informationRates)) {
      str10 = paste0("Spending time: ",
                     paste(round(x$settings$spendingTime, 3), collapse = ","))
      df1 = data.frame(x = rep("", 10))
      colnames(df1) = NULL
      rownames(df1) = c(str1, str2, str3, str4, str5, str6, str7,
                        paste(str8, str9, sep = ", "), str10, "")
    } else {
      df1 = data.frame(x = rep("", 9))
      colnames(df1) = NULL
      rownames(df1) = c(str1, str2, str3, str4, str5, str6, str7,
                        paste(str8, str9, sep = ", "), "")
    }
  } else {
    df1 = data.frame(x = rep("", 8))
    colnames(df1) = NULL
    rownames(df1) = c(str1, str2, str3, str4, str5, str6, str7, "")
  }

  if (k>1) {
    b <- s[, c("informationRates", "efficacyBounds", "futilityBounds",
               "cumulativeRejection", "cumulativeFutility",
               "cumulativeAlphaSpent", "numberOfSubjects",
               "efficacyPLarger", "futilityPLarger",
               "efficacyP", "futilityP", "information",
               "cumulativeRejectionH0", "cumulativeFutilityH0")]

    # format number of digits after decimal for each column
    j1 <- 7
    j2 <- 12
    j3 <- c(1,2,3,4,14)
    j4 <- c(5,6,8,9,10,11,13)

    b[j1] <- lapply(b[j1], formatC, format = "f", digits = 1)
    b[j2] <- lapply(b[j2], formatC, format = "f", digits = 2)
    b[j3] <- lapply(b[j3], formatC, format = "f", digits = 3)
    b[j4] <- lapply(b[j4], formatC, format = "f", digits = 4)

    if (x$settings$typeBetaSpending != 'none' ||
        (k > 1 && any(x$byStageResults$futilityBounds[1:(k-1)] > -6))) {
      df = t(b)
      rownames(df) = c("Information rate",
                       "Efficacy boundary (Z)",
                       "Futility boundary (Z)",
                       "Cumulative rejection",
                       "Cumulative futility",
                       "Cumulative alpha spent",
                       "Number of subjects",
                       "Efficacy boundary (pLarger)",
                       "Futility boundary (pLarger)",
                       "Efficacy boundary (p)",
                       "Futility boundary (p)",
                       "Information",
                       "Cumulative rejection under H0",
                       "Cumulative futility under H0")

    } else {
      df = t(b[,c(1,2,4,6,7,8,10,12)])
      rownames(df) = c("Information rate",
                       "Efficacy boundary (Z)",
                       "Cumulative rejection",
                       "Cumulative alpha spent",
                       "Number of subjects",
                       "Efficacy boundary (pLarger)",
                       "Efficacy boundary (p)",
                       "Information")
    }

    colnames(df) <- paste("Stage", seq_len(ncol(df)), sep=" ")
  } else {
    b <- s[, c("efficacyBounds", "efficacyPLarger", "efficacyP")]

    # format number of digits after decimal for each column
    j3 <- 1
    j4 <- c(2,3)

    b[j3] <- lapply(b[j3], formatC, format = "f", digits = 3)
    b[j4] <- lapply(b[j4], formatC, format = "f", digits = 4)

    df = t(b)

    rownames(df) = c("Efficacy boundary (Z)",
                     "Efficacy boundary (pLarger)",
                     "Efficacy boundary (p)")
    colnames(df) <- NA
  }

  print(df1, ..., na.print = "" , quote = FALSE )
  print(df, ..., na.print = "" , quote = FALSE )
  invisible(x)
}


#' @title Print Power and Sample Size Results for Two-Sample Mean Difference
#' at the Last Time Point From the MMRM Model
#' @description Prints the summary statistics from power calculation of
#' two-sample mean difference at the last time point from the MMRM model.
#'
#' @param x The designMeanDiffMMRM object to print.
#' @param ... Ensures that all arguments starting from "..." are named.
#'
#' @return A tabular printout of the design elements.
#'
#' @keywords internal
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @export
print.designMeanDiffMMRM <- function(x, ...) {
  a = x$overallResults
  s = x$byStageResults
  kMax = a$kMax

  if (kMax>1) {
    str1 = paste0("Group-sequential design with ", kMax, " stages")
  } else {
    str1 = "Fixed design"
  }

  str1 <- paste(str1, "for two-sample mean difference",
                "from the MMRM model")

  str2 <- paste0("Mean difference under H0: ",
                 round(a$meanDiffH0, 3), ", ",
                 "mean difference under H1: ",
                 round(a$meanDiff, 3))

  str3 <- paste0("Standard deviation for treatment: ",
                 round(sqrt(x$settings$covar1[x$settings$k,x$settings$k]),
                       3), ", ",
                 "standard deviation for control: ",
                 round(sqrt(x$settings$covar2[x$settings$k,x$settings$k]),
                       3))

  str4 <- paste0("Overall power: ",
                 round(a$overallReject, 3), ", ",
                 "overall alpha (1-sided): ",
                 round(a$alpha, 4))

  if (x$settings$typeBetaSpending != 'none' ||
      (kMax > 1 && any(s$futilityBounds[1:(kMax-1)] > -6))) {
    str4 <- paste0(str4, ", ",
                   "attained alpha: ", round(a$attainedAlpha, 4))
  }

  str5 <- paste0("Drift parameter: ", round(a$drift, 3), ", ",
                 "inflation factor: ", round(a$inflationFactor, 3))

  if (kMax>1) {
    str6 <- paste0("Maximum information: ",
                   round(a$information, 2), ", ",
                   "expected under H1: ",
                   round(a$expectedInformationH1, 2), ", ",
                   "expected under H0: ",
                   round(a$expectedInformationH0, 2))

    str7 <- paste0("Maximum # subjects: ",
                   round(a$numberOfSubjects, 1), ", ",
                   "expected under H1: ",
                   round(a$expectedNumberOfSubjectsH1, 1), ", ",
                   "expected under H0: ",
                   round(a$expectedNumberOfSubjectsH0, 1))

    str8 <- paste0("Total study duration: ",
                   round(a$studyDuration, 1), ", ",
                   "expected under H1: ",
                   round(a$expectedStudyDurationH1, 1), ", ",
                   "expected under H0: ",
                   round(a$expectedStudyDurationH0, 1))
  } else {
    str6 <- paste0("Information: ",
                   round(a$information, 2))

    str7 <- paste0("Number of subjects: ",
                   round(a$numberOfSubjects, 1))

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

  if (kMax > 1) {
    asf = tolower(x$settings$typeAlphaSpending)
    asfpar = round(x$settings$parameterAlphaSpending, 3)
    asfuser = round(x$settings$userAlphaSpending, 4)

    bsf = tolower(x$settings$typeBetaSpending)
    bsfpar = round(x$settings$parameterBetaSpending, 3)
    bsfuser = round(x$settings$userBetaSpending, 4)

    if (asf == "of") {
      str11 = paste0("Alpha spending: O'Brien-Fleming")
    } else if (asf == "p") {
      str11 = paste0("Alpha spending: Pocock")
    } else if (asf == "wt") {
      str11 = paste0("Alpha spending: Wang-Tsiatis(Delta = ", asfpar, ")")
    } else if (asf == "sfof") {
      str11 = paste0("Alpha spending: Lan-DeMets O'Brien-Fleming")
    } else if (asf == "sfp") {
      str11 = paste0("Alpha spending: Lan-DeMets Pocock")
    } else if (asf == "sfkd") {
      str11 = paste0("Alpha spending: KD(rho = ", asfpar, ")")
    } else if (asf == "sfhsd") {
      str11 = paste0("Alpha spending: HSD(gamma = ", asfpar, ")")
    } else if (asf == "user") {
      str11 = paste0("Alpha spending: User defined(",
                     paste(asfuser, collapse = ","), ")")
    } else {
      str11 = "Alpha spending: None"
    }

    if (bsf == "of") {
      str12 = paste0("beta spending: O'Brien-Fleming")
    } else if (bsf == "p") {
      str12 = paste0("beta spending: Pocock")
    } else if (bsf == "wt") {
      str12 = paste0("beta spending: Wang-Tsiatis(Delta = ", bsfpar, ")")
    } else if (bsf == "sfof") {
      str12 = paste0("beta spending: Lan-DeMets O'Brien-Fleming")
    } else if (bsf == "sfp") {
      str12 = paste0("beta spending: Lan-DeMets Pocock")
    } else if (bsf == "sfkd") {
      str12 = paste0("beta spending: KD(rho = ", bsfpar, ")")
    } else if (bsf == "sfhsd") {
      str12 = paste0("beta spending: HSD(gamma = ", bsfpar, ")")
    } else if (bsf == "user") {
      str12 = paste0("beta spending: User defined(",
                     paste(bsfuser, collapse = ","), ")")
    } else {
      str12 = "beta spending: None"
    }

    if (!any(is.na(x$settings$spendingTime)) &&
        !all(x$settings$spendingTime == s$informationRates)) {
      str13 = paste0("Spending time: ",
                     paste(round(x$settings$spendingTime, 3), collapse = ","))
      df1 = data.frame(x = rep("", 13))
      colnames(df1) = NULL
      rownames(df1) = c(str1, str2, str3, str4, str5, str6, str7, str8,
                        str9, str10, paste(str11, str12, sep = ", "),
                        str13, "")
    } else {
      df1 = data.frame(x = rep("", 12))
      colnames(df1) = NULL
      rownames(df1) = c(str1, str2, str3, str4, str5, str6, str7, str8,
                        str9, str10, paste(str11, str12, sep = ", "),
                        "")
    }
  } else {
    df1 = data.frame(x = rep("", 11))
    colnames(df1) = NULL
    rownames(df1) = c(str1, str2, str3, str4, str5, str6, str7, str8,
                      str9, str10, "")
  }

  if (kMax>1) {
    b <- s[, c("informationRates", "efficacyBounds", "futilityBounds",
               "cumulativeRejection", "cumulativeFutility",
               "cumulativeAlphaSpent", "numberOfSubjects",
               "numberOfCompleters", "analysisTime",
               "efficacyMeanDiff", "futilityMeanDiff",
               "efficacyP", "futilityP", "information",
               "cumulativeRejectionH0", "cumulativeFutilityH0")]

    # format number of digits after decimal for each column
    j1 <- c(7,8,9)
    j2 <- 14
    j3 <- c(1,2,3,4,10,11,16)
    j4 <- c(5,6,12,13,15)

    b[j1] <- lapply(b[j1], formatC, format = "f", digits = 1)
    b[j2] <- lapply(b[j2], formatC, format = "f", digits = 2)
    b[j3] <- lapply(b[j3], formatC, format = "f", digits = 3)
    b[j4] <- lapply(b[j4], formatC, format = "f", digits = 4)

    if (x$settings$typeBetaSpending != 'none' ||
        (kMax > 1 && any(s$futilityBounds[1:(kMax-1)] > -6))) {
      df = t(b)
      rownames(df) = c("Information rate",
                       "Efficacy boundary (Z)",
                       "Futility boundary (Z)",
                       "Cumulative rejection",
                       "Cumulative futility",
                       "Cumulative alpha spent",
                       "Number of subjects",
                       "Number of completers",
                       "Analysis time",
                       "Efficacy boundary (mean diff)",
                       "Futility boundary (mean diff)",
                       "Efficacy boundary (p)",
                       "Futility boundary (p)",
                       "Information",
                       "Cumulative rejection under H0",
                       "Cumulative futility under H0")

    } else {
      df = t(b[,c(1,2,4,6,7,8,9,10,12,14)])
      rownames(df) = c("Information rate",
                       "Efficacy boundary (Z)",
                       "Cumulative rejection",
                       "Cumulative alpha spent",
                       "Number of subjects",
                       "Number of completers",
                       "Analysis time",
                       "Efficacy boundary (mean diff)",
                       "Efficacy boundary (p)",
                       "Information")
    }

    colnames(df) <- paste("Stage", seq_len(ncol(df)), sep=" ")
  } else {
    b <- s[, c("efficacyBounds", "efficacyMeanDiff", "efficacyP")]

    # format number of digits after decimal for each column
    j3 <- c(1,2)
    j4 <- 3

    b[j3] <- lapply(b[j3], formatC, format = "f", digits = 3)
    b[j4] <- lapply(b[j4], formatC, format = "f", digits = 4)

    df = t(b)

    if (x$settings$normalApproximation) {
      rownames(df) = c("Efficacy boundary (Z)",
                       "Efficacy boundary (mean diff)",
                       "Efficacy boundary (p)")
    } else {
      rownames(df) = c("Efficacy boundary (t)",
                       "Efficacy boundary (mean diff)",
                       "Efficacy boundary (p)")
    }

    colnames(df) <- NA
  }

  print(df1, ..., na.print = "" , quote = FALSE )
  print(df, ..., na.print = "" , quote = FALSE )
  invisible(x)
}


#' @title Print Power and Sample Size Results for Direct Treatment Effects
#' in Crossover Trials With or Without Accounting for Carryover Effects
#' @description Prints the summary statistics from power calculation of
#' direct treatment effects in crossover trials with or without accounting
#' for carryover effects.
#'
#' @param x The designMeanDiffCarryover object to print.
#' @param ... Ensures that all arguments starting from "..." are named.
#'
#' @return A tabular printout of the design elements.
#'
#' @keywords internal
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @export
print.designMeanDiffCarryover <- function(x, ...) {
  carryover = x$carryover

  if (carryover) {
    str1 = paste("Testing direct treatment effects accounting for",
                 "carryover effects")
  } else {
    str1 = paste("Testing direct treatment effects without accounting for",
                 "carryover effects")
  }

  df1 = data.frame(x = "")
  colnames(df1) = NULL
  rownames(df1) = str1

  df2 <- as.data.frame(x$design)

  str2 <- paste0("Treatment comparison of interest: ",
                 x$trtpair[1], " - ", x$trtpair[2])

  str3 <- paste0("Mean difference under H0: ",
                 round(x$meanDiffH0, 3), ", ",
                 "mean difference under H1: ",
                 round(x$meanDiff, 3))

  str4 <- paste0("Within-subject standard deviation: ",
                 round(x$stDev, 3), ", ",
                 "intra-subject correlation: ",
                 round(x$corr, 3))

  str5 <- paste0("Power: ", round(x$power, 3), ", ",
                 "alpha (1-sided): ", round(x$alpha, 4), ", ",
                 "CI half-width: ", round(x$half_width, 3))

  str6 <- paste0("Cumulative dropout rates over periods: ",
                 paste(round(x$cumdrop, 3), collapse = ", "))

  str7 <- paste0("Without accounting for carryover effects, ",
                 "variance for direct effect: ",
                 round(x$v_direct_only, 3))

  if (carryover) {
    str8 <- paste0("Accounting for carryover, ",
                   "variance for direct effect: ",
                   round(x$v_direct, 3), ", ",
                   "for carryover: ",
                   round(x$v_carry, 3))

    str9 <- paste0("Relative efficiency for direct effect: ",
                   round(x$releff_direct, 3), ", ",
                   "for carryover effect: ",
                   round(x$releff_carry, 3))
  }

  str10 <- paste0("Number of subjects: ",
                 round(x$numberOfSubjects, 1))

  str11 <- paste0("Sequence allocation ratio: ",
                  paste(round(x$allocationRatioPlanned), collapse = " "))

  str12 <- paste0("Test statistic: ",
                  ifelse(x$normalApproximation, "z-test", "t-test"))

  if (carryover) {
    df3 = data.frame(x = rep("", 11))
    colnames(df3) = NULL
    rownames(df3) = c(str2, str3, str4, str5, str6, str7, str8,
                      str9, str10, str11, str12)

  } else {
    df3 = data.frame(x = rep("", 9))
    colnames(df3) = NULL
    rownames(df3) = c(str2, str3, str4, str5, str6, str7,
                      str10, str11, str12)
  }

  print(df1, ..., na.print = "" , quote = FALSE )
  print(df2, ..., na.print = "" , quote = FALSE )
  print(df3, ..., na.print = "" , quote = FALSE )
  invisible(x)
}


#' @title Print Power and Sample Size Results for Equivalence in Direct
#' Treatment Effects in Crossover Trials With or Without Accounting for
#' Carryover Effects
#' @description Prints the summary statistics from power calculation of
#' equivalence in direct treatment effects in crossover trials with or
#' without accounting for carryover effects.
#'
#' @param x The designMeanDiffCarryoverEquiv object to print.
#' @param ... Ensures that all arguments starting from "..." are named.
#'
#' @return A tabular printout of the design elements.
#'
#' @keywords internal
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @export
print.designMeanDiffCarryoverEquiv <- function(x, ...) {
  carryover = x$carryover

  if (carryover) {
    str1 = paste("Testing equivalence in direct treatment effects",
                 "accounting for carryover effects")
  } else {
    str1 = paste("Testing equivalence in direct treatment effects",
                 "without accounting for carryover effects")
  }

  df1 = data.frame(x = "")
  colnames(df1) = NULL
  rownames(df1) = str1

  df2 <- as.data.frame(x$design)

  str2 <- paste0("Treatment comparison of interest: ",
                 x$trtpair[1], " - ", x$trtpair[2])

  str3 <- paste0("Lower equivalence limit: ",
                 round(x$meanDiffLower, 3), ", ",
                 "upper equivalence limit: ",
                 round(x$meanDiffUpper, 3))

  str4 <- paste0("Mean difference under H1: ",
                 round(x$meanDiff, 3))

  str5 <- paste0("Within-subject standard deviation: ",
                 round(x$stDev, 3), ", ",
                 "intra-subject correlation: ",
                 round(x$corr, 3))

  str6 <- paste0("Power: ", round(x$power, 3), ", ",
                 "alpha (1-sided): ", round(x$alpha, 4), ", ",
                 "CI half-width: ", round(x$half_width, 3))

  str7 <- paste0("Cumulative dropout rates over periods: ",
                 paste(round(x$cumdrop, 3), collapse = ", "))

  str8 <- paste0("Without accounting for carryover effects, ",
                 "variance for direct effect: ",
                 round(x$v_direct_only, 3))

  if (carryover) {
    str9 <- paste0("Accounting for carryover, ",
                   "variance for direct effect: ",
                   round(x$v_direct, 3), ", ",
                   "for carryover: ",
                   round(x$v_carry, 3))

    str10 <- paste0("Relative efficiency for direct effect: ",
                    round(x$releff_direct, 3), ", ",
                    "for carryover effect: ",
                    round(x$releff_carry, 3))
  }

  str11 <- paste0("Number of subjects: ",
                  round(x$numberOfSubjects, 1))

  str12 <- paste0("Sequence allocation ratio: ",
                  paste(round(x$allocationRatioPlanned), collapse = " "))

  str13 <- paste0("Test statistic: ",
                  ifelse(x$normalApproximation, "z-test", "t-test"))

  if (carryover) {
    df3 = data.frame(x = rep("", 12))
    colnames(df3) = NULL
    rownames(df3) = c(str2, str3, str4, str5, str6, str7, str8,
                      str9, str10, str11, str12, str13)

  } else {
    df3 = data.frame(x = rep("", 10))
    colnames(df3) = NULL
    rownames(df3) = c(str2, str3, str4, str5, str6, str7, str8,
                      str11, str12, str13)
  }

  print(df1, ..., na.print = "" , quote = FALSE )
  print(df2, ..., na.print = "" , quote = FALSE )
  print(df3, ..., na.print = "" , quote = FALSE )
  invisible(x)
}


#' @title Print Power and Sample Size Results for One-Way ANOVA
#' @description Prints the power and sample size for one-way analysis
#' of variance.
#'
#' @param x The designANOVA object to print.
#' @param ... Ensures that all arguments starting from "..." are named.
#'
#' @return A tabular printout of the design elements.
#'
#' @keywords internal
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @export
print.designANOVA <- function(x, ...) {
  df1 = data.frame(alpha = x$alpha, power = x$power,
                   n = x$n, ngroups = x$ngroups,
                   stDev = x$stDev, effectsize = x$effectsize)
  rownames(df1) = NULL
  print(df1, ..., na.print = "" , quote = FALSE )
  invisible(x)
}


#' @title Print Power and Sample Size Results for a Single Contrast in
#' One-Way ANOVA
#' @description Prints the power and sample size for a single contrast in
#' one-way analysis of variance.
#'
#' @param x The designANOVAContrast object to print.
#' @param ... Ensures that all arguments starting from "..." are named.
#'
#' @return A tabular printout of the design elements.
#'
#' @keywords internal
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @export
print.designANOVAContrast <- function(x, ...) {
  df1 = data.frame(alpha = x$alpha, power = x$power,
                   n = x$n, ngroups = x$ngroups, stDev = x$stDev,
                   meanContrastH0 = x$meanContrastH0,
                   meanContrast = x$meanContrast,
                   effectsize = x$effectsize)
  rownames(df1) = NULL
  print(df1, ..., na.print = "" , quote = FALSE )
  invisible(x)
}


#' @title Print Power and Sample Size Results for One-Way Repeated
#' Measures ANOVA
#' @description Prints the power and sample size for one-way repeated
#' measures analysis of variance.
#'
#' @param x The designRepeatedANOVA object to print.
#' @param ... Ensures that all arguments starting from "..." are named.
#'
#' @return A tabular printout of the design elements.
#'
#' @keywords internal
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @export
print.designRepeatedANOVA <- function(x, ...) {
  df1 = data.frame(alpha = x$alpha, power = x$power,
                   n = x$n, ngroups = x$ngroups,
                   stDev = x$stDev, corr = x$corr,
                   effectsize = x$effectsize)
  rownames(df1) = NULL
  print(df1, ..., na.print = "" , quote = FALSE )
  invisible(x)
}


#' @title Print Power and Sample Size Results for a Single Contrast in
#' One-Way Repeated Measures ANOVA
#' @description Prints the power and sample size for a single contrast in
#' one-way repeated measures analysis of variance.
#'
#' @param x The designRepeatedANOVAContrast object to print.
#' @param ... Ensures that all arguments starting from "..." are named.
#'
#' @return A tabular printout of the design elements.
#'
#' @keywords internal
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @export
print.designRepeatedANOVAContrast <- function(x, ...) {
  df1 = data.frame(alpha = x$alpha, power = x$power,
                   n = x$n, ngroups = x$ngroups, stDev = x$stDev,
                   corr = x$corr,
                   meanContrastH0 = x$meanContrastH0,
                   meanContrast = x$meanContrast,
                   effectsize = x$effectsize)
  rownames(df1) = NULL
  print(df1, ..., na.print = "" , quote = FALSE )
  invisible(x)
}


#' @title Print Power and Sample Size Results for Two-Way ANOVA
#' @description Prints the power and sample size for two-way analysis
#' of variance.
#'
#' @param x The designTwoWayANOVA object to print.
#' @param ... Ensures that all arguments starting from "..." are named.
#'
#' @return A tabular printout of the design elements.
#'
#' @keywords internal
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @export
print.designTwoWayANOVA <- function(x, ...) {
  df1 = x$powerdf
  df1$alpha = x$alpha
  df1$nlevelsA = x$nlevelsA
  df1$nlevelsB = x$nlevelsB
  df1$stDev = x$stDev
  df1$effectsizeA = x$effectsizeA
  df1$effectsizeB = x$effectsizeB
  df1$effectsizeAB = x$effectsizeAB
  rownames(df1) = NULL
  print(df1, ..., na.print = "" , quote = FALSE )
  invisible(x)
}


#' @title Print Power and Sample Size Results for One-Sample Slope
#' @description Prints the summary statistics from power calculation of
#' one-sample slope.
#'
#' @param x The designOneSlope object to print.
#' @param ... Ensures that all arguments starting from "..." are named.
#'
#' @return A tabular printout of the design elements.
#'
#' @keywords internal
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @export
print.designOneSlope <- function(x, ...) {
  a = x$overallResults
  s = x$byStageResults
  k = a$kMax

  if (k>1) {
    str1 = paste0("Group-sequential design with ", k, " stages")
  } else {
    str1 = "Fixed design"
  }

  str1 <- paste(str1, "for one-sample slope")

  str2 <- paste0("Slope under H0: ",
                 round(a$slopeH0, 3), ", ",
                 "slope under H1: ",
                 round(a$slope, 3))

  str3 <- paste0("Standard deviation of residual: ",
                 round(a$stDev, 3), ", ",
                 "standard deviation of covariate: ",
                 round(a$stDevCovariate, 3))

  str4 <- paste0("Overall power: ",
                 round(a$overallReject, 3), ", ",
                 "overall alpha (1-sided): ",
                 round(a$alpha, 4))

  if (x$settings$typeBetaSpending != 'none' ||
      (k > 1 && any(x$byStageResults$futilityBounds[1:(k-1)] > -6))) {
    str4 <- paste0(str4, ", ",
                   "attained alpha: ", round(a$attainedAlpha, 4))
  }

  str5 <- paste0("Drift parameter: ", round(a$drift, 3), ", ",
                 "inflation factor: ", round(a$inflationFactor, 3))


  if (k>1) {
    str6 <- paste0("Maximum information: ",
                   round(a$information, 2), ", ",
                   "expected under H1: ",
                   round(a$expectedInformationH1, 2), ", ",
                   "expected under H0: ",
                   round(a$expectedInformationH0, 2))

    str7 <- paste0("Maximum # subjects: ",
                   round(a$numberOfSubjects, 1), ", ",
                   "expected under H1: ",
                   round(a$expectedNumberOfSubjectsH1, 1), ", ",
                   "expected under H0: ",
                   round(a$expectedNumberOfSubjectsH0, 1))

  } else {
    str6 <- paste0("Information: ",
                   round(a$information, 2))

    str7 <- paste0("Number of subjects: ",
                   round(a$numberOfSubjects, 1))
  }

  if (k > 1) {
    asf = tolower(x$settings$typeAlphaSpending)
    asfpar = round(x$settings$parameterAlphaSpending, 3)
    asfuser = round(x$settings$userAlphaSpending, 4)

    bsf = tolower(x$settings$typeBetaSpending)
    bsfpar = round(x$settings$parameterBetaSpending, 3)
    bsfuser = round(x$settings$userBetaSpending, 4)

    if (asf == "of") {
      str8 = paste0("Alpha spending: O'Brien-Fleming")
    } else if (asf == "p") {
      str8 = paste0("Alpha spending: Pocock")
    } else if (asf == "wt") {
      str8 = paste0("Alpha spending: Wang-Tsiatis(Delta = ", asfpar, ")")
    } else if (asf == "sfof") {
      str8 = paste0("Alpha spending: Lan-DeMets O'Brien-Fleming")
    } else if (asf == "sfp") {
      str8 = paste0("Alpha spending: Lan-DeMets Pocock")
    } else if (asf == "sfkd") {
      str8 = paste0("Alpha spending: KD(rho = ", asfpar, ")")
    } else if (asf == "sfhsd") {
      str8 = paste0("Alpha spending: HSD(gamma = ", asfpar, ")")
    } else if (asf == "user") {
      str8 = paste0("Alpha spending: User defined(",
                    paste(asfuser, collapse = ","), ")")
    } else {
      str8 = "Alpha spending: None"
    }

    if (bsf == "of") {
      str9 = paste0("beta spending: O'Brien-Fleming")
    } else if (bsf == "p") {
      str9 = paste0("beta spending: Pocock")
    } else if (bsf == "wt") {
      str9 = paste0("beta spending: Wang-Tsiatis(Delta = ", bsfpar, ")")
    } else if (bsf == "sfof") {
      str9 = paste0("beta spending: Lan-DeMets O'Brien-Fleming")
    } else if (bsf == "sfp") {
      str9 = paste0("beta spending: Lan-DeMets Pocock")
    } else if (bsf == "sfkd") {
      str9 = paste0("beta spending: KD(rho = ", bsfpar, ")")
    } else if (bsf == "sfhsd") {
      str9 = paste0("beta spending: HSD(gamma = ", bsfpar, ")")
    } else if (bsf == "user") {
      str9 = paste0("beta spending: User defined(",
                    paste(bsfuser, collapse = ","), ")")
    } else {
      str9 = "beta spending: None"
    }

    if (!any(is.na(x$settings$spendingTime)) &&
        !all(x$settings$spendingTime == s$informationRates)) {
      str10 = paste0("Spending time: ",
                     paste(round(x$settings$spendingTime, 3), collapse = ","))
      df1 = data.frame(x = rep("", 10))
      colnames(df1) = NULL
      rownames(df1) = c(str1, str2, str3, str4, str5, str6, str7,
                        paste(str8, str9, sep = ", "), str10, "")
    } else {
      df1 = data.frame(x = rep("", 9))
      colnames(df1) = NULL
      rownames(df1) = c(str1, str2, str3, str4, str5, str6, str7,
                        paste(str8, str9, sep = ", "), "")
    }
  } else {
    df1 = data.frame(x = rep("", 8))
    colnames(df1) = NULL
    rownames(df1) = c(str1, str2, str3, str4, str5, str6, str7, "")
  }

  if (k>1) {
    b <- s[, c("informationRates", "efficacyBounds", "futilityBounds",
               "cumulativeRejection", "cumulativeFutility",
               "cumulativeAlphaSpent", "numberOfSubjects",
               "efficacySlope", "futilitySlope",
               "efficacyP", "futilityP", "information",
               "cumulativeRejectionH0", "cumulativeFutilityH0")]

    # format number of digits after decimal for each column
    j1 <- 7
    j2 <- 12
    j3 <- c(1,2,3,4,8,9,14)
    j4 <- c(5,6,10,11,13)

    b[j1] <- lapply(b[j1], formatC, format = "f", digits = 1)
    b[j2] <- lapply(b[j2], formatC, format = "f", digits = 2)
    b[j3] <- lapply(b[j3], formatC, format = "f", digits = 3)
    b[j4] <- lapply(b[j4], formatC, format = "f", digits = 4)

    if (x$settings$typeBetaSpending != 'none' ||
        (k > 1 && any(x$byStageResults$futilityBounds[1:(k-1)] > -6))) {
      df = t(b)
      rownames(df) = c("Information rate",
                       "Efficacy boundary (Z)",
                       "Futility boundary (Z)",
                       "Cumulative rejection",
                       "Cumulative futility",
                       "Cumulative alpha spent",
                       "Number of subjects",
                       "Efficacy boundary (slope)",
                       "Futility boundary (slope)",
                       "Efficacy boundary (p)",
                       "Futility boundary (p)",
                       "Information",
                       "Cumulative rejection under H0",
                       "Cumulative futility under H0")

    } else {
      df = t(b[,c(1,2,4,6,7,8,10,12)])
      rownames(df) = c("Information rate",
                       "Efficacy boundary (Z)",
                       "Cumulative rejection",
                       "Cumulative alpha spent",
                       "Number of subjects",
                       "Efficacy boundary (slope)",
                       "Efficacy boundary (p)",
                       "Information")
    }

    colnames(df) <- paste("Stage", seq_len(ncol(df)), sep=" ")
  } else {
    b <- s[, c("efficacyBounds", "efficacySlope", "efficacyP")]

    # format number of digits after decimal for each column
    j3 <- c(1,2)
    j4 <- 3

    b[j3] <- lapply(b[j3], formatC, format = "f", digits = 3)
    b[j4] <- lapply(b[j4], formatC, format = "f", digits = 4)

    df = t(b)

    if (!x$settings$normalApproximation) {
      rownames(df) = c("Efficacy boundary (t)",
                       "Efficacy boundary (slope)",
                       "Efficacy boundary (p)")
    } else {
      rownames(df) = c("Efficacy boundary (Z)",
                       "Efficacy boundary (slope)",
                       "Efficacy boundary (p)")
    }

    colnames(df) <- NA
  }

  print(df1, ..., na.print = "" , quote = FALSE )
  print(df, ..., na.print = "" , quote = FALSE )
  invisible(x)
}


#' @title Print Power and Sample Size Results for Two-Sample Slope Difference
#' @description Prints the summary statistics from power calculation of
#' two-sample slope difference.
#'
#' @param x The designSlopeDiff object to print.
#' @param ... Ensures that all arguments starting from "..." are named.
#'
#' @return A tabular printout of the design elements.
#'
#' @keywords internal
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @export
print.designSlopeDiff <- function(x, ...) {
  a = x$overallResults
  s = x$byStageResults
  k = a$kMax

  if (k>1) {
    str1 = paste0("Group-sequential design with ", k, " stages")
  } else {
    str1 = "Fixed design"
  }

  str1 <- paste(str1, "for two-sample slope difference")

  str2 <- paste0("Slope difference under H0: ",
                 round(a$slopeDiffH0, 3), ", ",
                 "slope difference under H1: ",
                 round(a$slopeDiff, 3))

  str3 <- paste0("Standard deviation of residual: ",
                 round(a$stDev, 3), ", ",
                 "standard deviation of covariate: ",
                 round(a$stDevCovariate, 3))

  str4 <- paste0("Overall power: ",
                 round(a$overallReject, 3), ", ",
                 "overall alpha (1-sided): ",
                 round(a$alpha, 4))

  if (x$settings$typeBetaSpending != 'none' ||
      (k > 1 && any(x$byStageResults$futilityBounds[1:(k-1)] > -6))) {
    str4 <- paste0(str4, ", ",
                   "attained alpha: ", round(a$attainedAlpha, 4))
  }

  str5 <- paste0("Drift parameter: ", round(a$drift, 3), ", ",
                 "inflation factor: ", round(a$inflationFactor, 3))

  if (k>1) {
    str6 <- paste0("Maximum information: ",
                   round(a$information, 2), ", ",
                   "expected under H1: ",
                   round(a$expectedInformationH1, 2), ", ",
                   "expected under H0: ",
                   round(a$expectedInformationH0, 2))

    str7 <- paste0("Maximum # subjects: ",
                   round(a$numberOfSubjects, 1), ", ",
                   "expected under H1: ",
                   round(a$expectedNumberOfSubjectsH1, 1), ", ",
                   "expected under H0: ",
                   round(a$expectedNumberOfSubjectsH0, 1))

  } else {
    str6 <- paste0("Information: ",
                   round(a$information, 2))

    str7 <- paste0("Number of subjects: ",
                   round(a$numberOfSubjects, 1))
  }

  str8 <- paste0("Allocation ratio: ",
                 round(x$settings$allocationRatioPlanned, 3))

  if (k > 1) {
    asf = tolower(x$settings$typeAlphaSpending)
    asfpar = round(x$settings$parameterAlphaSpending, 3)
    asfuser = round(x$settings$userAlphaSpending, 4)

    bsf = tolower(x$settings$typeBetaSpending)
    bsfpar = round(x$settings$parameterBetaSpending, 3)
    bsfuser = round(x$settings$userBetaSpending, 4)

    if (asf == "of") {
      str9 = paste0("Alpha spending: O'Brien-Fleming")
    } else if (asf == "p") {
      str9 = paste0("Alpha spending: Pocock")
    } else if (asf == "wt") {
      str9 = paste0("Alpha spending: Wang-Tsiatis(Delta = ", asfpar, ")")
    } else if (asf == "sfof") {
      str9 = paste0("Alpha spending: Lan-DeMets O'Brien-Fleming")
    } else if (asf == "sfp") {
      str9 = paste0("Alpha spending: Lan-DeMets Pocock")
    } else if (asf == "sfkd") {
      str9 = paste0("Alpha spending: KD(rho = ", asfpar, ")")
    } else if (asf == "sfhsd") {
      str9 = paste0("Alpha spending: HSD(gamma = ", asfpar, ")")
    } else if (asf == "user") {
      str9 = paste0("Alpha spending: User defined(",
                    paste(asfuser, collapse = ","), ")")
    } else {
      str9 = "Alpha spending: None"
    }

    if (bsf == "of") {
      str10 = paste0("beta spending: O'Brien-Fleming")
    } else if (bsf == "p") {
      str10 = paste0("beta spending: Pocock")
    } else if (bsf == "wt") {
      str10 = paste0("beta spending: Wang-Tsiatis(Delta = ", bsfpar, ")")
    } else if (bsf == "sfof") {
      str10 = paste0("beta spending: Lan-DeMets O'Brien-Fleming")
    } else if (bsf == "sfp") {
      str10 = paste0("beta spending: Lan-DeMets Pocock")
    } else if (bsf == "sfkd") {
      str10 = paste0("beta spending: KD(rho = ", bsfpar, ")")
    } else if (bsf == "sfhsd") {
      str10 = paste0("beta spending: HSD(gamma = ", bsfpar, ")")
    } else if (bsf == "user") {
      str10 = paste0("beta spending: User defined(",
                     paste(bsfuser, collapse = ","), ")")
    } else {
      str10 = "beta spending: None"
    }

    if (!any(is.na(x$settings$spendingTime)) &&
        !all(x$settings$spendingTime == s$informationRates)) {
      str11 = paste0("Spending time: ",
                     paste(round(x$settings$spendingTime, 3), collapse = ","))
      df1 = data.frame(x = rep("", 11))
      colnames(df1) = NULL
      rownames(df1) = c(str1, str2, str3, str4, str5, str6, str7, str8,
                        paste(str9, str10, sep = ", "), str11, "")
    } else {
      df1 = data.frame(x = rep("", 10))
      colnames(df1) = NULL
      rownames(df1) = c(str1, str2, str3, str4, str5, str6, str7, str8,
                        paste(str9, str10, sep = ", "), "")
    }
  } else {
    df1 = data.frame(x = rep("", 9))
    colnames(df1) = NULL
    rownames(df1) = c(str1, str2, str3, str4, str5, str6, str7, str8, "")
  }

  if (k>1) {
    b <- s[, c("informationRates", "efficacyBounds", "futilityBounds",
               "cumulativeRejection", "cumulativeFutility",
               "cumulativeAlphaSpent", "numberOfSubjects",
               "efficacySlopeDiff", "futilitySlopeDiff",
               "efficacyP", "futilityP", "information",
               "cumulativeRejectionH0", "cumulativeFutilityH0")]

    # format number of digits after decimal for each column
    j1 <- 7
    j2 <- 12
    j3 <- c(1,2,3,4,8,9,14)
    j4 <- c(5,6,10,11,13)

    b[j1] <- lapply(b[j1], formatC, format = "f", digits = 1)
    b[j2] <- lapply(b[j2], formatC, format = "f", digits = 2)
    b[j3] <- lapply(b[j3], formatC, format = "f", digits = 3)
    b[j4] <- lapply(b[j4], formatC, format = "f", digits = 4)

    if (x$settings$typeBetaSpending != 'none' ||
        (k > 1 && any(x$byStageResults$futilityBounds[1:(k-1)] > -6))) {
      df = t(b)
      rownames(df) = c("Information rate",
                       "Efficacy boundary (Z)",
                       "Futility boundary (Z)",
                       "Cumulative rejection",
                       "Cumulative futility",
                       "Cumulative alpha spent",
                       "Number of subjects",
                       "Efficacy boundary (slope diff)",
                       "Futility boundary (slope diff)",
                       "Efficacy boundary (p)",
                       "Futility boundary (p)",
                       "Information",
                       "Cumulative rejection under H0",
                       "Cumulative futility under H0")

    } else {
      df = t(b[,c(1,2,4,6,7,8,10,12)])
      rownames(df) = c("Information rate",
                       "Efficacy boundary (Z)",
                       "Cumulative rejection",
                       "Cumulative alpha spent",
                       "Number of subjects",
                       "Efficacy boundary (slope diff)",
                       "Efficacy boundary (p)",
                       "Information")
    }

    colnames(df) <- paste("Stage", seq_len(ncol(df)), sep=" ")
  } else {
    b <- s[, c("efficacyBounds", "efficacySlopeDiff", "efficacyP")]

    # format number of digits after decimal for each column
    j3 <- c(1,2)
    j4 <- 3

    b[j3] <- lapply(b[j3], formatC, format = "f", digits = 3)
    b[j4] <- lapply(b[j4], formatC, format = "f", digits = 4)

    df = t(b)

    if (!x$settings$normalApproximation) {
      rownames(df) = c("Efficacy boundary (t)",
                       "Efficacy boundary (slope diff)",
                       "Efficacy boundary (p)")
    } else {
      rownames(df) = c("Efficacy boundary (Z)",
                       "Efficacy boundary (slope diff)",
                       "Efficacy boundary (p)")
    }

    colnames(df) <- NA
  }

  print(df1, ..., na.print = "" , quote = FALSE )
  print(df, ..., na.print = "" , quote = FALSE )
  invisible(x)
}


#' @title Print Power and Sample Size Results for Two-Sample Slope Difference
#' at the Last Time Point From the MMRM Model
#' @description Prints the summary statistics from power calculation of
#' two-sample slope difference at the last time point from the MMRM model.
#'
#' @param x The designSlopeDiffMMRM object to print.
#' @param ... Ensures that all arguments starting from "..." are named.
#'
#' @return A tabular printout of the design elements.
#'
#' @keywords internal
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @export
print.designSlopeDiffMMRM <- function(x, ...) {
  a = x$overallResults
  s = x$byStageResults
  kMax = a$kMax

  if (kMax>1) {
    str1 = paste0("Group-sequential design with ", kMax, " stages")
  } else {
    str1 = "Fixed design"
  }

  str1 <- paste(str1, "for two-sample slope difference",
                "from the MMRM model")

  str2 <- paste0("Slope difference under H0: ",
                 round(a$slopeDiffH0, 3), ", ",
                 "slope difference under H1: ",
                 round(a$slopeDiff, 3))

  str3 <- paste0("Standard deviation (SD) of within-subject residual: ",
                 round(x$settings$stDev, 3))

  str4 <- paste0("SD of random intercept: ",
                 round(sqrt(x$settings$G[1,1]), 3), ", ",
                 "SD of random slope: ",
                 round(sqrt(x$settings$G[2,2]), 3), ", ",
                 "correlation: ",
                 round(x$settings$G[1,2]/sqrt(x$settings$G[1,1]*
                                                x$settings$G[2,2]), 3))

  str5 <- paste0("Overall power: ",
                 round(a$overallReject, 3), ", ",
                 "overall alpha (1-sided): ",
                 round(a$alpha, 4))

  if (x$settings$typeBetaSpending != 'none' ||
      (kMax > 1 && any(s$futilityBounds[1:(kMax-1)] > -6))) {
    str5 <- paste0(str5, ", ",
                   "attained alpha: ", round(a$attainedAlpha, 4))
  }

  str6 <- paste0("Drift parameter: ", round(a$drift, 3), ", ",
                 "inflation factor: ", round(a$inflationFactor, 3))

  if (kMax>1) {
    str7 <- paste0("Maximum information: ",
                   round(a$information, 2), ", ",
                   "expected under H1: ",
                   round(a$expectedInformationH1, 2), ", ",
                   "expected under H0: ",
                   round(a$expectedInformationH0, 2))

    str8 <- paste0("Maximum # subjects: ",
                   round(a$numberOfSubjects, 1), ", ",
                   "expected under H1: ",
                   round(a$expectedNumberOfSubjectsH1, 1), ", ",
                   "expected under H0: ",
                   round(a$expectedNumberOfSubjectsH0, 1))

    str9 <- paste0("Total study duration: ",
                   round(a$studyDuration, 1), ", ",
                   "expected under H1: ",
                   round(a$expectedStudyDurationH1, 1), ", ",
                   "expected under H0: ",
                   round(a$expectedStudyDurationH0, 1))
  } else {
    str7 <- paste0("Information: ",
                   round(a$information, 2))

    str8 <- paste0("Number of subjects: ",
                   round(a$numberOfSubjects, 1))

    str9 <- paste0("Study duration: ",
                   round(a$studyDuration, 1))
  }

  str10 <- paste0("Accrual duration: ",
                  round(a$accrualDuration, 1), ", ",
                  "follow-up duration: ",
                  round(a$followupTime, 1), ", ",
                  "fixed follow-up: ", a$fixedFollowup)

  str11 <- paste0("Allocation ratio: ",
                  round(x$settings$allocationRatioPlanned, 3))

  if (kMax > 1) {
    asf = tolower(x$settings$typeAlphaSpending)
    asfpar = round(x$settings$parameterAlphaSpending, 3)
    asfuser = round(x$settings$userAlphaSpending, 4)

    bsf = tolower(x$settings$typeBetaSpending)
    bsfpar = round(x$settings$parameterBetaSpending, 3)
    bsfuser = round(x$settings$userBetaSpending, 4)

    if (asf == "of") {
      str12 = paste0("Alpha spending: O'Brien-Fleming")
    } else if (asf == "p") {
      str12 = paste0("Alpha spending: Pocock")
    } else if (asf == "wt") {
      str12 = paste0("Alpha spending: Wang-Tsiatis(Delta = ", asfpar, ")")
    } else if (asf == "sfof") {
      str12 = paste0("Alpha spending: Lan-DeMets O'Brien-Fleming")
    } else if (asf == "sfp") {
      str12 = paste0("Alpha spending: Lan-DeMets Pocock")
    } else if (asf == "sfkd") {
      str12 = paste0("Alpha spending: KD(rho = ", asfpar, ")")
    } else if (asf == "sfhsd") {
      str12 = paste0("Alpha spending: HSD(gamma = ", asfpar, ")")
    } else if (asf == "user") {
      str12 = paste0("Alpha spending: User defined(",
                     paste(asfuser, collapse = ","), ")")
    } else {
      str12 = "Alpha spending: None"
    }

    if (bsf == "of") {
      str13 = paste0("beta spending: O'Brien-Fleming")
    } else if (bsf == "p") {
      str13 = paste0("beta spending: Pocock")
    } else if (bsf == "wt") {
      str13 = paste0("beta spending: Wang-Tsiatis(Delta = ", bsfpar, ")")
    } else if (bsf == "sfof") {
      str13 = paste0("beta spending: Lan-DeMets O'Brien-Fleming")
    } else if (bsf == "sfp") {
      str13 = paste0("beta spending: Lan-DeMets Pocock")
    } else if (bsf == "sfkd") {
      str13 = paste0("beta spending: KD(rho = ", bsfpar, ")")
    } else if (bsf == "sfhsd") {
      str13 = paste0("beta spending: HSD(gamma = ", bsfpar, ")")
    } else if (bsf == "user") {
      str13 = paste0("beta spending: User defined(",
                     paste(bsfuser, collapse = ","), ")")
    } else {
      str13 = "beta spending: None"
    }

    if (!any(is.na(x$settings$spendingTime)) &&
        !all(x$settings$spendingTime == s$informationRates)) {
      str14 = paste0("Spending time: ",
                     paste(round(x$settings$spendingTime, 3), collapse = ","))
      df1 = data.frame(x = rep("", 14))
      colnames(df1) = NULL
      rownames(df1) = c(str1, str2, str3, str4, str5, str6, str7, str8,
                        str9, str10, str11, paste(str12, str13, sep = ", "),
                        str14, "")
    } else {
      df1 = data.frame(x = rep("", 13))
      colnames(df1) = NULL
      rownames(df1) = c(str1, str2, str3, str4, str5, str6, str7, str8,
                        str9, str10, str11, paste(str12, str13, sep = ", "),
                        "")
    }
  } else {
    df1 = data.frame(x = rep("", 12))
    colnames(df1) = NULL
    rownames(df1) = c(str1, str2, str3, str4, str5, str6, str7, str8,
                      str9, str10, str11, "")
  }

  if (kMax>1) {
    if (!a$fixedFollowup) {
      b <- s[, c("informationRates", "efficacyBounds", "futilityBounds",
                 "cumulativeRejection", "cumulativeFutility",
                 "cumulativeAlphaSpent", "numberOfSubjects",
                 "analysisTime",
                 "efficacySlopeDiff", "futilitySlopeDiff",
                 "efficacyP", "futilityP", "information",
                 "cumulativeRejectionH0", "cumulativeFutilityH0")]

      # format number of digits after decimal for each column
      j1 <- c(7,8)
      j2 <- 13
      j3 <- c(1,2,3,4,9,10,15)
      j4 <- c(5,6,11,12,14)

      b[j1] <- lapply(b[j1], formatC, format = "f", digits = 1)
      b[j2] <- lapply(b[j2], formatC, format = "f", digits = 2)
      b[j3] <- lapply(b[j3], formatC, format = "f", digits = 3)
      b[j4] <- lapply(b[j4], formatC, format = "f", digits = 4)

      if (x$settings$typeBetaSpending != 'none' ||
          (kMax > 1 && any(s$futilityBounds[1:(kMax-1)] > -6))) {
        df = t(b)
        rownames(df) = c("Information rate",
                         "Efficacy boundary (Z)",
                         "Futility boundary (Z)",
                         "Cumulative rejection",
                         "Cumulative futility",
                         "Cumulative alpha spent",
                         "Number of subjects",
                         "Analysis time",
                         "Efficacy boundary (slope diff)",
                         "Futility boundary (slope diff)",
                         "Efficacy boundary (p)",
                         "Futility boundary (p)",
                         "Information",
                         "Cumulative rejection under H0",
                         "Cumulative futility under H0")

      } else {
        df = t(b[,c(1,2,4,6,7,8,9,11,13)])
        rownames(df) = c("Information rate",
                         "Efficacy boundary (Z)",
                         "Cumulative rejection",
                         "Cumulative alpha spent",
                         "Number of subjects",
                         "Analysis time",
                         "Efficacy boundary (slope diff)",
                         "Efficacy boundary (p)",
                         "Information")
      }
    } else {
      b <- s[, c("informationRates", "efficacyBounds", "futilityBounds",
                 "cumulativeRejection", "cumulativeFutility",
                 "cumulativeAlphaSpent", "numberOfSubjects",
                 "numberOfCompleters", "analysisTime",
                 "efficacySlopeDiff", "futilitySlopeDiff",
                 "efficacyP", "futilityP", "information",
                 "cumulativeRejectionH0", "cumulativeFutilityH0")]

      # format number of digits after decimal for each column
      j1 <- c(7,8,9)
      j2 <- 13
      j3 <- c(1,2,3,4,10,11,16)
      j4 <- c(5,6,12,13,15)

      b[j1] <- lapply(b[j1], formatC, format = "f", digits = 1)
      b[j2] <- lapply(b[j2], formatC, format = "f", digits = 2)
      b[j3] <- lapply(b[j3], formatC, format = "f", digits = 3)
      b[j4] <- lapply(b[j4], formatC, format = "f", digits = 4)

      if (x$settings$typeBetaSpending != 'none' ||
          (kMax > 1 && any(s$futilityBounds[1:(kMax-1)] > -6))) {
        df = t(b)
        rownames(df) = c("Information rate",
                         "Efficacy boundary (Z)",
                         "Futility boundary (Z)",
                         "Cumulative rejection",
                         "Cumulative futility",
                         "Cumulative alpha spent",
                         "Number of subjects",
                         "Number of completers",
                         "Analysis time",
                         "Efficacy boundary (slope diff)",
                         "Futility boundary (slope diff)",
                         "Efficacy boundary (p)",
                         "Futility boundary (p)",
                         "Information",
                         "Cumulative rejection under H0",
                         "Cumulative futility under H0")

      } else {
        df = t(b[,c(1,2,4,6,7,8,9,10,12,14)])
        rownames(df) = c("Information rate",
                         "Efficacy boundary (Z)",
                         "Cumulative rejection",
                         "Cumulative alpha spent",
                         "Number of subjects",
                         "Number of completers",
                         "Analysis time",
                         "Efficacy boundary (slope diff)",
                         "Efficacy boundary (p)",
                         "Information")
      }
    }

    colnames(df) <- paste("Stage", seq_len(ncol(df)), sep=" ")
  } else {
    b <- s[, c("efficacyBounds", "efficacySlopeDiff", "efficacyP")]

    # format number of digits after decimal for each column
    j3 <- c(1,2)
    j4 <- 3

    b[j3] <- lapply(b[j3], formatC, format = "f", digits = 3)
    b[j4] <- lapply(b[j4], formatC, format = "f", digits = 4)

    df = t(b)

    if (x$settings$normalApproximation) {
      rownames(df) = c("Efficacy boundary (Z)",
                       "Efficacy boundary (slope diff)",
                       "Efficacy boundary (p)")
    } else {
      rownames(df) = c("Efficacy boundary (t)",
                       "Efficacy boundary (slope diff)",
                       "Efficacy boundary (p)")
    }

    colnames(df) <- NA
  }

  print(df1, ..., na.print = "" , quote = FALSE )
  print(df, ..., na.print = "" , quote = FALSE )
  invisible(x)
}


#' @title Print Power and Sample Size Results for One-Sample Proportion
#' @description Prints the summary statistics from power calculation of
#' one-sample proportion.
#'
#' @param x The designOneProportion object to print.
#' @param ... Ensures that all arguments starting from "..." are named.
#'
#' @return A tabular printout of the design elements.
#'
#' @keywords internal
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @export
print.designOneProportion <- function(x, ...) {
  a = x$overallResults
  s = x$byStageResults
  k = a$kMax

  if (k>1) {
    str1 = paste0("Group-sequential design with ", k, " stages")
  } else {
    str1 = "Fixed design"
  }

  str1 <- paste(str1, "for one-sample proportion")

  str2 <- paste0("Response probability under H0: ",
                 round(a$piH0, 3), ", ",
                 "response probability under H1: ",
                 round(a$pi, 3))

  str3 <- paste0("Overall power: ",
                 round(a$overallReject, 3), ", ",
                 "overall alpha (1-sided): ",
                 round(a$alpha, 4))

  if (x$settings$typeBetaSpending != 'none' ||
      (k > 1 && any(x$byStageResults$futilityBounds[1:(k-1)] > -6)) ||
      (k == 1 && !x$settings$normalApproximation)) {
    str3 <- paste0(str3, ", ",
                   "attained alpha: ", round(a$attainedAlpha, 4))
  }

  str4 <- paste0("Drift parameter: ", round(a$drift, 3), ", ",
                 "inflation factor: ", round(a$inflationFactor, 3))

  if (k>1) {
    str5 <- paste0("Maximum information: ",
                   round(a$information, 2), ", ",
                   "expected under H1: ",
                   round(a$expectedInformationH1, 2), ", ",
                   "expected under H0: ",
                   round(a$expectedInformationH0, 2))

    str6 <- paste0("Maximum # subjects: ",
                   round(a$numberOfSubjects, 1), ", ",
                   "expected under H1: ",
                   round(a$expectedNumberOfSubjectsH1, 1), ", ",
                   "expected under H0: ",
                   round(a$expectedNumberOfSubjectsH0, 1))

  } else {
    str5 <- paste0("Information: ",
                   round(a$information, 2))

    str6 <- paste0("Number of subjects: ",
                   round(a$numberOfSubjects, 1))
  }

  str7 <- paste0("Test statistic: ",
                 ifelse(x$settings$normalApproximation, "z-test",
                        "exact test"))

  if (k > 1) {
    asf = tolower(x$settings$typeAlphaSpending)
    asfpar = round(x$settings$parameterAlphaSpending, 3)
    asfuser = round(x$settings$userAlphaSpending, 4)

    bsf = tolower(x$settings$typeBetaSpending)
    bsfpar = round(x$settings$parameterBetaSpending, 3)
    bsfuser = round(x$settings$userBetaSpending, 4)

    if (asf == "of") {
      str8 = paste0("Alpha spending: O'Brien-Fleming")
    } else if (asf == "p") {
      str8 = paste0("Alpha spending: Pocock")
    } else if (asf == "wt") {
      str8 = paste0("Alpha spending: Wang-Tsiatis(Delta = ", asfpar, ")")
    } else if (asf == "sfof") {
      str8 = paste0("Alpha spending: Lan-DeMets O'Brien-Fleming")
    } else if (asf == "sfp") {
      str8 = paste0("Alpha spending: Lan-DeMets Pocock")
    } else if (asf == "sfkd") {
      str8 = paste0("Alpha spending: KD(rho = ", asfpar, ")")
    } else if (asf == "sfhsd") {
      str8 = paste0("Alpha spending: HSD(gamma = ", asfpar, ")")
    } else if (asf == "user") {
      str8 = paste0("Alpha spending: User defined(",
                    paste(asfuser, collapse = ","), ")")
    } else {
      str8 = "Alpha spending: None"
    }

    if (bsf == "of") {
      str9 = paste0("beta spending: O'Brien-Fleming")
    } else if (bsf == "p") {
      str9 = paste0("beta spending: Pocock")
    } else if (bsf == "wt") {
      str9 = paste0("beta spending: Wang-Tsiatis(Delta = ", bsfpar, ")")
    } else if (bsf == "sfof") {
      str9 = paste0("beta spending: Lan-DeMets O'Brien-Fleming")
    } else if (bsf == "sfp") {
      str9 = paste0("beta spending: Lan-DeMets Pocock")
    } else if (bsf == "sfkd") {
      str9 = paste0("beta spending: KD(rho = ", bsfpar, ")")
    } else if (bsf == "sfhsd") {
      str9 = paste0("beta spending: HSD(gamma = ", bsfpar, ")")
    } else if (bsf == "user") {
      str9 = paste0("beta spending: User defined(",
                    paste(bsfuser, collapse = ","), ")")
    } else {
      str9 = "beta spending: None"
    }

    if (!any(is.na(x$settings$spendingTime)) &&
        !all(x$settings$spendingTime == s$informationRates)) {
      str10 = paste0("Spending time: ",
                     paste(round(x$settings$spendingTime, 3), collapse = ","))
      df1 = data.frame(x = rep("", 10))
      colnames(df1) = NULL
      rownames(df1) = c(str1, str2, str3, str4, str5, str6, str7,
                        paste(str8, str9, sep = ", "), str10, "")
    } else {
      df1 = data.frame(x = rep("", 9))
      colnames(df1) = NULL
      rownames(df1) = c(str1, str2, str3, str4, str5, str6, str7,
                        paste(str8, str9, sep = ", "), "")
    }
  } else {
    df1 = data.frame(x = rep("", 8))
    colnames(df1) = NULL
    rownames(df1) = c(str1, str2, str3, str4, str5, str6, str7, "")
  }

  if (k>1) {
    b <- s[, c("informationRates", "efficacyBounds", "futilityBounds",
               "cumulativeRejection", "cumulativeFutility",
               "cumulativeAlphaSpent", "numberOfSubjects",
               "efficacyResponses", "futilityResponses",
               "efficacyP", "futilityP", "information",
               "cumulativeRejectionH0", "cumulativeFutilityH0")]

    # format number of digits after decimal for each column
    j1 <- c(7,8,9)
    j2 <- 12
    j3 <- c(1,2,3,4,14)
    j4 <- c(5,6,10,11,13)

    b[j1] <- lapply(b[j1], formatC, format = "f", digits = 1)
    b[j2] <- lapply(b[j2], formatC, format = "f", digits = 2)
    b[j3] <- lapply(b[j3], formatC, format = "f", digits = 3)
    b[j4] <- lapply(b[j4], formatC, format = "f", digits = 4)

    if (x$settings$typeBetaSpending != 'none' ||
        (k > 1 && any(x$byStageResults$futilityBounds[1:(k-1)] > -6))) {
      df = t(b)
      rownames(df) = c("Information rate",
                       "Efficacy boundary (Z)",
                       "Futility boundary (Z)",
                       "Cumulative rejection",
                       "Cumulative futility",
                       "Cumulative alpha spent",
                       "Number of subjects",
                       "Efficacy boundary (# responses)",
                       "Futility boundary (# responses)",
                       "Efficacy boundary (p)",
                       "Futility boundary (p)",
                       "Information",
                       "Cumulative rejection under H0",
                       "Cumulative futility under H0")

    } else {
      df = t(b[,c(1,2,4,6,7,8,10,12)])
      rownames(df) = c("Information rate",
                       "Efficacy boundary (Z)",
                       "Cumulative rejection",
                       "Cumulative alpha spent",
                       "Number of subjects",
                       "Efficacy boundary (# responses)",
                       "Efficacy boundary (p)",
                       "Information")
    }

    colnames(df) <- paste("Stage", seq_len(ncol(df)), sep=" ")
  } else {
    b <- s[, c("efficacyBounds", "efficacyResponses", "efficacyP")]

    # format number of digits after decimal for each column
    j1 <- 2
    j3 <- 1
    j4 <- 3

    b[j1] <- lapply(b[j1], formatC, format = "f", digits = 1)
    b[j3] <- lapply(b[j3], formatC, format = "f", digits = 3)
    b[j4] <- lapply(b[j4], formatC, format = "f", digits = 4)

    df = t(b)

    if (!x$settings$normalApproximation) {
      rownames(df) = c("Efficacy boundary (Z)",
                       "Efficacy boundary (# responses)",
                       "Efficacy boundary (p)")
    } else {
      rownames(df) = c("Efficacy boundary (Z)",
                       "Efficacy boundary (# responses)",
                       "Efficacy boundary (p)")
    }

    colnames(df) <- NA
  }

  print(df1, ..., na.print = "" , quote = FALSE )
  print(df, ..., na.print = "" , quote = FALSE )
  invisible(x)
}


#' @title Print Power and Sample Size Results for McNemar's Test for Paired
#' Proportions
#' @description Prints the summary statistics from power calculation of
#' McNemar's test for paired proportions.
#'
#' @param x The designPairedPropMcNemar object to print.
#' @param ... Ensures that all arguments starting from "..." are named.
#'
#' @return A tabular printout of the design elements.
#'
#' @keywords internal
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @export
print.designPairedPropMcNemar <- function(x, ...) {
  a = x$overallResults
  s = x$byStageResults
  k = a$kMax

  if (k>1) {
    str1 = paste0("Group-sequential design with ", k, " stages")
  } else {
    str1 = "Fixed design"
  }

  str1 <- paste(str1, "for McNemar's test")

  str2 <- paste0("Proportion of discordant pairs: ",
                 round(a$pDiscordant, 3), ", ",
                 "risk difference: ", round(a$riskDiff, 3))

  str3 <- paste0("Overall power: ",
                 round(a$overallReject, 3), ", ",
                 "overall alpha (1-sided): ",
                 round(a$alpha, 4))

  if (x$settings$typeBetaSpending != 'none' ||
      (k > 1 && any(x$byStageResults$futilityBounds[1:(k-1)] > -6))) {
    str3 <- paste0(str3, ", ",
                   "attained alpha: ", round(a$attainedAlpha, 4))
  }

  str4 <- paste0("Drift parameter: ", round(a$drift, 3), ", ",
                 "inflation factor: ", round(a$inflationFactor, 3))

  if (k>1) {
    str5 <- paste0("Maximum information: ",
                   round(a$information, 2), ", ",
                   "expected under H1: ",
                   round(a$expectedInformationH1, 2), ", ",
                   "expected under H0: ",
                   round(a$expectedInformationH0, 2))

    str6 <- paste0("Maximum # subjects: ",
                   round(a$numberOfSubjects, 1), ", ",
                   "expected under H1: ",
                   round(a$expectedNumberOfSubjectsH1, 1), ", ",
                   "expected under H0: ",
                   round(a$expectedNumberOfSubjectsH0, 1))

  } else {
    str5 <- paste0("Information: ",
                   round(a$information, 2))

    str6 <- paste0("Number of subjects: ",
                   round(a$numberOfSubjects, 1))
  }

  str7 <- paste0("Variance of standardized test statistic: ",
                 ifelse(x$settings$nullVariance, "under H0", "under H1"))

  if (k > 1) {
    asf = tolower(x$settings$typeAlphaSpending)
    asfpar = round(x$settings$parameterAlphaSpending, 3)
    asfuser = round(x$settings$userAlphaSpending, 4)

    bsf = tolower(x$settings$typeBetaSpending)
    bsfpar = round(x$settings$parameterBetaSpending, 3)
    bsfuser = round(x$settings$userBetaSpending, 4)

    if (asf == "of") {
      str8 = paste0("Alpha spending: O'Brien-Fleming")
    } else if (asf == "p") {
      str8 = paste0("Alpha spending: Pocock")
    } else if (asf == "wt") {
      str8 = paste0("Alpha spending: Wang-Tsiatis(Delta = ", asfpar, ")")
    } else if (asf == "sfof") {
      str8 = paste0("Alpha spending: Lan-DeMets O'Brien-Fleming")
    } else if (asf == "sfp") {
      str8 = paste0("Alpha spending: Lan-DeMets Pocock")
    } else if (asf == "sfkd") {
      str8 = paste0("Alpha spending: KD(rho = ", asfpar, ")")
    } else if (asf == "sfhsd") {
      str8 = paste0("Alpha spending: HSD(gamma = ", asfpar, ")")
    } else if (asf == "user") {
      str8 = paste0("Alpha spending: User defined(",
                    paste(asfuser, collapse = ","), ")")
    } else {
      str8 = "Alpha spending: None"
    }

    if (bsf == "of") {
      str9 = paste0("beta spending: O'Brien-Fleming")
    } else if (bsf == "p") {
      str9 = paste0("beta spending: Pocock")
    } else if (bsf == "wt") {
      str9 = paste0("beta spending: Wang-Tsiatis(Delta = ", bsfpar, ")")
    } else if (bsf == "sfof") {
      str9 = paste0("beta spending: Lan-DeMets O'Brien-Fleming")
    } else if (bsf == "sfp") {
      str9 = paste0("beta spending: Lan-DeMets Pocock")
    } else if (bsf == "sfkd") {
      str9 = paste0("beta spending: KD(rho = ", bsfpar, ")")
    } else if (bsf == "sfhsd") {
      str9 = paste0("beta spending: HSD(gamma = ", bsfpar, ")")
    } else if (bsf == "user") {
      str9 = paste0("beta spending: User defined(",
                    paste(bsfuser, collapse = ","), ")")
    } else {
      str9 = "beta spending: None"
    }

    if (!any(is.na(x$settings$spendingTime)) &&
        !all(x$settings$spendingTime == s$informationRates)) {
      str10 = paste0("Spending time: ",
                     paste(round(x$settings$spendingTime, 3), collapse = ","))
      df1 = data.frame(x = rep("", 10))
      colnames(df1) = NULL
      rownames(df1) = c(str1, str2, str3, str4, str5, str6, str7,
                        paste(str8, str9, sep = ", "), str10, "")
    } else {
      df1 = data.frame(x = rep("", 9))
      colnames(df1) = NULL
      rownames(df1) = c(str1, str2, str3, str4, str5, str6, str7,
                        paste(str8, str8, sep = ", "), "")
    }
  } else {
    df1 = data.frame(x = rep("", 8))
    colnames(df1) = NULL
    rownames(df1) = c(str1, str2, str3, str4, str5, str6, str7, "")
  }

  if (k>1) {
    b <- s[, c("informationRates", "efficacyBounds", "futilityBounds",
               "cumulativeRejection", "cumulativeFutility",
               "cumulativeAlphaSpent", "numberOfSubjects",
               "efficacyRiskDiff", "futilityRiskDiff",
               "efficacyP", "futilityP", "information",
               "cumulativeRejectionH0", "cumulativeFutilityH0")]

    # format number of digits after decimal for each column
    j1 <- 7
    j2 <- 12
    j3 <- c(1,2,3,4,14)
    j4 <- c(5,6,8,9,10,11,13)

    b[j1] <- lapply(b[j1], formatC, format = "f", digits = 1)
    b[j2] <- lapply(b[j2], formatC, format = "f", digits = 2)
    b[j3] <- lapply(b[j3], formatC, format = "f", digits = 3)
    b[j4] <- lapply(b[j4], formatC, format = "f", digits = 4)

    if (x$settings$typeBetaSpending != 'none' ||
        (k > 1 && any(x$byStageResults$futilityBounds[1:(k-1)] > -6))) {
      df = t(b)
      rownames(df) = c("Information rate",
                       "Efficacy boundary (Z)",
                       "Futility boundary (Z)",
                       "Cumulative rejection",
                       "Cumulative futility",
                       "Cumulative alpha spent",
                       "Number of subjects",
                       "Efficacy boundary (risk diff)",
                       "Futility boundary (risk diff)",
                       "Efficacy boundary (p)",
                       "Futility boundary (p)",
                       "Information",
                       "Cumulative rejection under H0",
                       "Cumulative futility under H0")

    } else {
      df = t(b[,c(1,2,4,6,7,8,10,12)])
      rownames(df) = c("Information rate",
                       "Efficacy boundary (Z)",
                       "Cumulative rejection",
                       "Cumulative alpha spent",
                       "Number of subjects",
                       "Efficacy boundary (risk diff)",
                       "Efficacy boundary (p)",
                       "Information")
    }

    colnames(df) <- paste("Stage", seq_len(ncol(df)), sep=" ")
  } else {
    b <- s[, c("efficacyBounds", "efficacyRiskDiff", "efficacyP")]

    # format number of digits after decimal for each column
    j3 <- 1
    j4 <- c(2,3)

    b[j3] <- lapply(b[j3], formatC, format = "f", digits = 3)
    b[j4] <- lapply(b[j4], formatC, format = "f", digits = 4)

    df = t(b)

    rownames(df) = c("Efficacy boundary (Z)",
                     "Efficacy boundary (risk diff)",
                     "Efficacy boundary (p)")
    colnames(df) <- NA
  }

  print(df1, ..., na.print = "" , quote = FALSE )
  print(df, ..., na.print = "" , quote = FALSE )
  invisible(x)
}


#' @title Print Power and Sample Size Results for Two-Sample Risk Difference
#' @description Prints the summary statistics from power calculation of
#' two-sample risk difference.
#'
#' @param x The designRiskDiff object to print.
#' @param ... Ensures that all arguments starting from "..." are named.
#'
#' @return A tabular printout of the design elements.
#'
#' @keywords internal
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @export
print.designRiskDiff <- function(x, ...) {
  a = x$overallResults
  s = x$byStageResults
  k = a$kMax

  if (k>1) {
    str1 = paste0("Group-sequential design with ", k, " stages")
  } else {
    str1 = "Fixed design"
  }

  str1 <- paste(str1,
                "for two-sample risk difference")

  str2 <- paste0("Risk difference under H0: ", round(a$riskDiffH0, 3), ", ",
                 "proportion on treatment: ", round(a$pi1, 3), ", ",
                 "proportion on control: ", round(a$pi2, 3))

  str3 <- paste0("Overall power: ",
                 round(a$overallReject, 3), ", ",
                 "overall alpha (1-sided): ",
                 round(a$alpha, 4))

  if (x$settings$typeBetaSpending != 'none' ||
      (k > 1 && any(x$byStageResults$futilityBounds[1:(k-1)] > -6))) {
    str3 <- paste0(str3, ", ",
                   "attained alpha: ", round(a$attainedAlpha, 4))
  }

  str4 <- paste0("Drift parameter: ", round(a$drift, 3), ", ",
                 "inflation factor: ", round(a$inflationFactor, 3))

  if (k>1) {
    str5 <- paste0("Maximum information: ",
                   round(a$information, 2), ", ",
                   "expected under H1: ",
                   round(a$expectedInformationH1, 2), ", ",
                   "expected under H0: ",
                   round(a$expectedInformationH0, 2))

    str6 <- paste0("Maximum # subjects: ",
                   round(a$numberOfSubjects, 1), ", ",
                   "expected under H1: ",
                   round(a$expectedNumberOfSubjectsH1, 1), ", ",
                   "expected under H0: ",
                   round(a$expectedNumberOfSubjectsH0, 1))

  } else {
    str5 <- paste0("Information: ",
                   round(a$information, 2))

    str6 <- paste0("Number of subjects: ",
                   round(a$numberOfSubjects, 1))
  }

  str7 <- paste0("Allocation ratio: ",
                 round(x$settings$allocationRatioPlanned, 3), ", ",
                 "variance of standardized test statistic: ",
                 ifelse(x$settings$nullVariance, "under H0", "under H1"))

  if (k > 1) {
    asf = tolower(x$settings$typeAlphaSpending)
    asfpar = round(x$settings$parameterAlphaSpending, 3)
    asfuser = round(x$settings$userAlphaSpending, 4)

    bsf = tolower(x$settings$typeBetaSpending)
    bsfpar = round(x$settings$parameterBetaSpending, 3)
    bsfuser = round(x$settings$userBetaSpending, 4)

    if (asf == "of") {
      str8 = paste0("Alpha spending: O'Brien-Fleming")
    } else if (asf == "p") {
      str8 = paste0("Alpha spending: Pocock")
    } else if (asf == "wt") {
      str8 = paste0("Alpha spending: Wang-Tsiatis(Delta = ", asfpar, ")")
    } else if (asf == "sfof") {
      str8 = paste0("Alpha spending: Lan-DeMets O'Brien-Fleming")
    } else if (asf == "sfp") {
      str8 = paste0("Alpha spending: Lan-DeMets Pocock")
    } else if (asf == "sfkd") {
      str8 = paste0("Alpha spending: KD(rho = ", asfpar, ")")
    } else if (asf == "sfhsd") {
      str8 = paste0("Alpha spending: HSD(gamma = ", asfpar, ")")
    } else if (asf == "user") {
      str8 = paste0("Alpha spending: User defined(",
                    paste(asfuser, collapse = ","), ")")
    } else {
      str8 = "Alpha spending: None"
    }

    if (bsf == "of") {
      str9 = paste0("beta spending: O'Brien-Fleming")
    } else if (bsf == "p") {
      str9 = paste0("beta spending: Pocock")
    } else if (bsf == "wt") {
      str9 = paste0("beta spending: Wang-Tsiatis(Delta = ", bsfpar, ")")
    } else if (bsf == "sfof") {
      str9 = paste0("beta spending: Lan-DeMets O'Brien-Fleming")
    } else if (bsf == "sfp") {
      str9 = paste0("beta spending: Lan-DeMets Pocock")
    } else if (bsf == "sfkd") {
      str9 = paste0("beta spending: KD(rho = ", bsfpar, ")")
    } else if (bsf == "sfhsd") {
      str9 = paste0("beta spending: HSD(gamma = ", bsfpar, ")")
    } else if (bsf == "user") {
      str9 = paste0("beta spending: User defined(",
                    paste(bsfuser, collapse = ","), ")")
    } else {
      str9 = "beta spending: None"
    }

    if (!any(is.na(x$settings$spendingTime)) &&
        !all(x$settings$spendingTime == s$informationRates)) {
      str10 = paste0("Spending time: ",
                     paste(round(x$settings$spendingTime, 3), collapse = ","))
      df1 = data.frame(x = rep("", 10))
      colnames(df1) = NULL
      rownames(df1) = c(str1, str2, str3, str4, str5, str6, str7,
                        paste(str8, str9, sep = ", "), str10, "")
    } else {
      df1 = data.frame(x = rep("", 9))
      colnames(df1) = NULL
      rownames(df1) = c(str1, str2, str3, str4, str5, str6, str7,
                        paste(str8, str9, sep = ", "), "")
    }
  } else {
    df1 = data.frame(x = rep("", 8))
    colnames(df1) = NULL
    rownames(df1) = c(str1, str2, str3, str4, str5, str6, str7, "")
  }

  if (k>1) {
    b <- s[, c("informationRates", "efficacyBounds", "futilityBounds",
               "cumulativeRejection", "cumulativeFutility",
               "cumulativeAlphaSpent", "numberOfSubjects",
               "efficacyRiskDiff", "futilityRiskDiff",
               "efficacyP", "futilityP", "information",
               "cumulativeRejectionH0", "cumulativeFutilityH0")]

    # format number of digits after decimal for each column
    j1 <- 7
    j2 <- 12
    j3 <- c(1,2,3,4,8,9,14)
    j4 <- c(5,6,10,11,13)

    b[j1] <- lapply(b[j1], formatC, format = "f", digits = 1)
    b[j2] <- lapply(b[j2], formatC, format = "f", digits = 2)
    b[j3] <- lapply(b[j3], formatC, format = "f", digits = 3)
    b[j4] <- lapply(b[j4], formatC, format = "f", digits = 4)

    if (x$settings$typeBetaSpending != 'none' ||
        (k > 1 && any(x$byStageResults$futilityBounds[1:(k-1)] > -6))) {
      df = t(b)
      rownames(df) = c("Information rate",
                       "Efficacy boundary (Z)",
                       "Futility boundary (Z)",
                       "Cumulative rejection",
                       "Cumulative futility",
                       "Cumulative alpha spent",
                       "Number of subjects",
                       "Efficacy boundary (risk diff)",
                       "Futility boundary (risk diff)",
                       "Efficacy boundary (p)",
                       "Futility boundary (p)",
                       "Information",
                       "Cumulative rejection under H0",
                       "Cumulative futility under H0")

    } else {
      df = t(b[,c(1,2,4,6,7,8,10,12)])
      rownames(df) = c("Information rate",
                       "Efficacy boundary (Z)",
                       "Cumulative rejection",
                       "Cumulative alpha spent",
                       "Number of subjects",
                       "Efficacy boundary (risk diff)",
                       "Efficacy boundary (p)",
                       "Information")
    }

    colnames(df) <- paste("Stage", seq_len(ncol(df)), sep=" ")
  } else {
    b <- s[, c("efficacyBounds", "efficacyRiskDiff", "efficacyP")]

    # format number of digits after decimal for each column
    j3 <- c(1,2)
    j4 <- 3

    b[j3] <- lapply(b[j3], formatC, format = "f", digits = 3)
    b[j4] <- lapply(b[j4], formatC, format = "f", digits = 4)

    df = t(b)

    rownames(df) = c("Efficacy boundary (Z)",
                     "Efficacy boundary (risk diff)",
                     "Efficacy boundary (p)")
    colnames(df) <- NA
  }

  print(df1, ..., na.print = "" , quote = FALSE )
  print(df, ..., na.print = "" , quote = FALSE )
  invisible(x)
}


#' @title Print Power and Sample Size Results for Two-Sample Risk Ratio
#' @description Prints the summary statistics from power calculation of
#' two-sample risk ratio.
#'
#' @param x The designRiskRatio object to print.
#' @param ... Ensures that all arguments starting from "..." are named.
#'
#' @return A tabular printout of the design elements.
#'
#' @keywords internal
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @export
print.designRiskRatio <- function(x, ...) {
  a = x$overallResults
  s = x$byStageResults
  k = a$kMax

  if (k>1) {
    str1 = paste0("Group-sequential design with ", k, " stages")
  } else {
    str1 = "Fixed design"
  }

  str1 <- paste(str1, "for two-sample risk ratio")

  str2 <- paste0("Risk ratio under H0: ", round(a$riskRatioH0, 3), ", ",
                 "proportion on treatment: ", round(a$pi1, 3), ", ",
                 "proportion on control: ", round(a$pi2, 3))

  str3 <- paste0("Overall power: ",
                 round(a$overallReject, 3), ", ",
                 "overall alpha (1-sided): ",
                 round(a$alpha, 4))

  if (x$settings$typeBetaSpending != 'none' ||
      (k > 1 && any(x$byStageResults$futilityBounds[1:(k-1)] > -6))) {
    str3 <- paste0(str3, ", ",
                   "attained alpha: ", round(a$attainedAlpha, 4))
  }

  str4 <- paste0("Drift parameter: ", round(a$drift, 3), ", ",
                 "inflation factor: ", round(a$inflationFactor, 3))

  if (k>1) {
    str5 <- paste0("Maximum information: ",
                   round(a$information, 2), ", ",
                   "expected under H1: ",
                   round(a$expectedInformationH1, 2), ", ",
                   "expected under H0: ",
                   round(a$expectedInformationH0, 2))

    str6 <- paste0("Maximum # subjects: ",
                   round(a$numberOfSubjects, 1), ", ",
                   "expected under H1: ",
                   round(a$expectedNumberOfSubjectsH1, 1), ", ",
                   "expected under H0: ",
                   round(a$expectedNumberOfSubjectsH0, 1))

  } else {
    str5 <- paste0("Information: ",
                   round(a$information, 2))

    str6 <- paste0("Number of subjects: ",
                   round(a$numberOfSubjects, 1))
  }

  str7 <- paste0("Allocation ratio: ",
                 round(x$settings$allocationRatioPlanned, 3), ", ",
                 "variance of standardized test statistic: ",
                 ifelse(x$settings$nullVariance, "under H0", "under H1"))

  if (k > 1) {
    asf = tolower(x$settings$typeAlphaSpending)
    asfpar = round(x$settings$parameterAlphaSpending, 3)
    asfuser = round(x$settings$userAlphaSpending, 4)

    bsf = tolower(x$settings$typeBetaSpending)
    bsfpar = round(x$settings$parameterBetaSpending, 3)
    bsfuser = round(x$settings$userBetaSpending, 4)

    if (asf == "of") {
      str8 = paste0("Alpha spending: O'Brien-Fleming")
    } else if (asf == "p") {
      str8 = paste0("Alpha spending: Pocock")
    } else if (asf == "wt") {
      str8 = paste0("Alpha spending: Wang-Tsiatis(Delta = ", asfpar, ")")
    } else if (asf == "sfof") {
      str8 = paste0("Alpha spending: Lan-DeMets O'Brien-Fleming")
    } else if (asf == "sfp") {
      str8 = paste0("Alpha spending: Lan-DeMets Pocock")
    } else if (asf == "sfkd") {
      str8 = paste0("Alpha spending: KD(rho = ", asfpar, ")")
    } else if (asf == "sfhsd") {
      str8 = paste0("Alpha spending: HSD(gamma = ", asfpar, ")")
    } else if (asf == "user") {
      str8 = paste0("Alpha spending: User defined(",
                    paste(asfuser, collapse = ","), ")")
    } else {
      str8 = "Alpha spending: None"
    }

    if (bsf == "of") {
      str9 = paste0("beta spending: O'Brien-Fleming")
    } else if (bsf == "p") {
      str9 = paste0("beta spending: Pocock")
    } else if (bsf == "wt") {
      str9 = paste0("beta spending: Wang-Tsiatis(Delta = ", bsfpar, ")")
    } else if (bsf == "sfof") {
      str9 = paste0("beta spending: Lan-DeMets O'Brien-Fleming")
    } else if (bsf == "sfp") {
      str9 = paste0("beta spending: Lan-DeMets Pocock")
    } else if (bsf == "sfkd") {
      str9 = paste0("beta spending: KD(rho = ", bsfpar, ")")
    } else if (bsf == "sfhsd") {
      str9 = paste0("beta spending: HSD(gamma = ", bsfpar, ")")
    } else if (bsf == "user") {
      str9 = paste0("beta spending: User defined(",
                    paste(bsfuser, collapse = ","), ")")
    } else {
      str9 = "beta spending: None"
    }

    if (!any(is.na(x$settings$spendingTime)) &&
        !all(x$settings$spendingTime == s$informationRates)) {
      str10 = paste0("Spending time: ",
                     paste(round(x$settings$spendingTime, 3), collapse = ","))
      df1 = data.frame(x = rep("", 10))
      colnames(df1) = NULL
      rownames(df1) = c(str1, str2, str3, str4, str5, str6, str7,
                        paste(str8, str9, sep = ", "), str10, "")
    } else {
      df1 = data.frame(x = rep("", 9))
      colnames(df1) = NULL
      rownames(df1) = c(str1, str2, str3, str4, str5, str6, str7,
                        paste(str8, str9, sep = ", "), "")
    }
  } else {
    df1 = data.frame(x = rep("", 8))
    colnames(df1) = NULL
    rownames(df1) = c(str1, str2, str3, str4, str5, str6, str7, "")
  }

  if (k>1) {
    b <- s[, c("informationRates", "efficacyBounds", "futilityBounds",
               "cumulativeRejection", "cumulativeFutility",
               "cumulativeAlphaSpent", "numberOfSubjects",
               "efficacyRiskRatio", "futilityRiskRatio",
               "efficacyP", "futilityP", "information",
               "cumulativeRejectionH0", "cumulativeFutilityH0")]

    # format number of digits after decimal for each column
    j1 <- 7
    j2 <- 12
    j3 <- c(1,2,3,4,8,9,14)
    j4 <- c(5,6,10,11,13)

    b[j1] <- lapply(b[j1], formatC, format = "f", digits = 1)
    b[j2] <- lapply(b[j2], formatC, format = "f", digits = 2)
    b[j3] <- lapply(b[j3], formatC, format = "f", digits = 3)
    b[j4] <- lapply(b[j4], formatC, format = "f", digits = 4)

    if (x$settings$typeBetaSpending != 'none' ||
        (k > 1 && any(x$byStageResults$futilityBounds[1:(k-1)] > -6))) {
      df = t(b)
      rownames(df) = c("Information rate",
                       "Efficacy boundary (Z)",
                       "Futility boundary (Z)",
                       "Cumulative rejection",
                       "Cumulative futility",
                       "Cumulative alpha spent",
                       "Number of subjects",
                       "Efficacy boundary (risk ratio)",
                       "Futility boundary (risk ratio)",
                       "Efficacy boundary (p)",
                       "Futility boundary (p)",
                       "Information",
                       "Cumulative rejection under H0",
                       "Cumulative futility under H0")

    } else {
      df = t(b[,c(1,2,4,6,7,8,10,12)])
      rownames(df) = c("Information rate",
                       "Efficacy boundary (Z)",
                       "Cumulative rejection",
                       "Cumulative alpha spent",
                       "Number of subjects",
                       "Efficacy boundary (risk ratio)",
                       "Efficacy boundary (p)",
                       "Information")
    }

    colnames(df) <- paste("Stage", seq_len(ncol(df)), sep=" ")
  } else {
    b <- s[, c("efficacyBounds", "efficacyRiskRatio", "efficacyP")]

    # format number of digits after decimal for each column
    j3 <- c(1,2)
    j4 <- 3

    b[j3] <- lapply(b[j3], formatC, format = "f", digits = 3)
    b[j4] <- lapply(b[j4], formatC, format = "f", digits = 4)

    df = t(b)

    rownames(df) = c("Efficacy boundary (Z)",
                     "Efficacy boundary (risk ratio)",
                     "Efficacy boundary (p)")
    colnames(df) <- NA
  }

  print(df1, ..., na.print = "" , quote = FALSE )
  print(df, ..., na.print = "" , quote = FALSE )
  invisible(x)
}


#' @title Print Power and Sample Size Results for Two-Sample Risk Ratio
#' Based on the Farrington-Manning Score Test
#' @description Prints the summary statistics from power calculation of
#' two-sample risk ratio based on the Farrington-Manning score test.
#'
#' @param x The designRiskRatioFM object to print.
#' @param ... Ensures that all arguments starting from "..." are named.
#'
#' @return A tabular printout of the design elements.
#'
#' @keywords internal
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @export
print.designRiskRatioFM <- function(x, ...) {
  a = x$overallResults
  s = x$byStageResults
  k = a$kMax

  if (k>1) {
    str1 = paste0("Group-sequential design with ", k, " stages")
  } else {
    str1 = "Fixed design"
  }

  str1 <- paste(str1, "for two-sample risk ratio based on the score test")

  str2 <- paste0("Risk ratio under H0: ", round(a$riskRatioH0, 3), ", ",
                 "proportion on treatment: ", round(a$pi1, 3), ", ",
                 "proportion on control: ", round(a$pi2, 3))

  str3 <- paste0("Overall power: ",
                 round(a$overallReject, 3), ", ",
                 "overall alpha (1-sided): ",
                 round(a$alpha, 4))

  if (x$settings$typeBetaSpending != 'none' ||
      (k > 1 && any(x$byStageResults$futilityBounds[1:(k-1)] > -6))) {
    str3 <- paste0(str3, ", ",
                   "attained alpha: ", round(a$attainedAlpha, 4))
  }

  str4 <- paste0("Drift parameter: ", round(a$drift, 3), ", ",
                 "inflation factor: ", round(a$inflationFactor, 3))

  if (k>1) {
    str5 <- paste0("Maximum information: ",
                   round(a$information, 2), ", ",
                   "expected under H1: ",
                   round(a$expectedInformationH1, 2), ", ",
                   "expected under H0: ",
                   round(a$expectedInformationH0, 2))

    str6 <- paste0("Maximum # subjects: ",
                   round(a$numberOfSubjects, 1), ", ",
                   "expected under H1: ",
                   round(a$expectedNumberOfSubjectsH1, 1), ", ",
                   "expected under H0: ",
                   round(a$expectedNumberOfSubjectsH0, 1))

  } else {
    str5 <- paste0("Information: ",
                   round(a$information, 2))

    str6 <- paste0("Number of subjects: ",
                   round(a$numberOfSubjects, 1))
  }

  str7 <- paste0("Allocation ratio: ",
                 round(x$settings$allocationRatioPlanned, 3), ", ",
                 "variance of standardized test statistic: ",
                 ifelse(x$settings$nullVariance, "under H0", "under H1"))

  if (k > 1) {
    asf = tolower(x$settings$typeAlphaSpending)
    asfpar = round(x$settings$parameterAlphaSpending, 3)
    asfuser = round(x$settings$userAlphaSpending, 4)

    bsf = tolower(x$settings$typeBetaSpending)
    bsfpar = round(x$settings$parameterBetaSpending, 3)
    bsfuser = round(x$settings$userBetaSpending, 4)

    if (asf == "of") {
      str8 = paste0("Alpha spending: O'Brien-Fleming")
    } else if (asf == "p") {
      str8 = paste0("Alpha spending: Pocock")
    } else if (asf == "wt") {
      str8 = paste0("Alpha spending: Wang-Tsiatis(Delta = ", asfpar, ")")
    } else if (asf == "sfof") {
      str8 = paste0("Alpha spending: Lan-DeMets O'Brien-Fleming")
    } else if (asf == "sfp") {
      str8 = paste0("Alpha spending: Lan-DeMets Pocock")
    } else if (asf == "sfkd") {
      str8 = paste0("Alpha spending: KD(rho = ", asfpar, ")")
    } else if (asf == "sfhsd") {
      str8 = paste0("Alpha spending: HSD(gamma = ", asfpar, ")")
    } else if (asf == "user") {
      str8 = paste0("Alpha spending: User defined(",
                    paste(asfuser, collapse = ","), ")")
    } else {
      str8 = "Alpha spending: None"
    }

    if (bsf == "of") {
      str9 = paste0("beta spending: O'Brien-Fleming")
    } else if (bsf == "p") {
      str9 = paste0("beta spending: Pocock")
    } else if (bsf == "wt") {
      str9 = paste0("beta spending: Wang-Tsiatis(Delta = ", bsfpar, ")")
    } else if (bsf == "sfof") {
      str9 = paste0("beta spending: Lan-DeMets O'Brien-Fleming")
    } else if (bsf == "sfp") {
      str9 = paste0("beta spending: Lan-DeMets Pocock")
    } else if (bsf == "sfkd") {
      str9 = paste0("beta spending: KD(rho = ", bsfpar, ")")
    } else if (bsf == "sfhsd") {
      str9 = paste0("beta spending: HSD(gamma = ", bsfpar, ")")
    } else if (bsf == "user") {
      str9 = paste0("beta spending: User defined(",
                    paste(bsfuser, collapse = ","), ")")
    } else {
      str9 = "beta spending: None"
    }

    if (!any(is.na(x$settings$spendingTime)) &&
        !all(x$settings$spendingTime == s$informationRates)) {
      str10 = paste0("Spending time: ",
                     paste(round(x$settings$spendingTime, 3), collapse = ","))
      df1 = data.frame(x = rep("", 10))
      colnames(df1) = NULL
      rownames(df1) = c(str1, str2, str3, str4, str5, str6, str7,
                        paste(str8, str9, sep = ", "), str10, "")
    } else {
      df1 = data.frame(x = rep("", 9))
      colnames(df1) = NULL
      rownames(df1) = c(str1, str2, str3, str4, str5, str6, str7,
                        paste(str8, str9, sep = ", "), "")
    }
  } else {
    df1 = data.frame(x = rep("", 8))
    colnames(df1) = NULL
    rownames(df1) = c(str1, str2, str3, str4, str5, str6, str7, "")
  }

  if (k>1) {
    b <- s[, c("informationRates", "efficacyBounds", "futilityBounds",
               "cumulativeRejection", "cumulativeFutility",
               "cumulativeAlphaSpent", "numberOfSubjects",
               "efficacyRiskRatioScore", "futilityRiskRatioScore",
               "efficacyP", "futilityP", "information",
               "cumulativeRejectionH0", "cumulativeFutilityH0")]

    # format number of digits after decimal for each column
    j1 <- 7
    j2 <- 12
    j3 <- c(1,2,3,4,8,9,14)
    j4 <- c(5,6,10,11,13)

    b[j1] <- lapply(b[j1], formatC, format = "f", digits = 1)
    b[j2] <- lapply(b[j2], formatC, format = "f", digits = 2)
    b[j3] <- lapply(b[j3], formatC, format = "f", digits = 3)
    b[j4] <- lapply(b[j4], formatC, format = "f", digits = 4)

    if (x$settings$typeBetaSpending != 'none' ||
        (k > 1 && any(x$byStageResults$futilityBounds[1:(k-1)] > -6))) {
      df = t(b)
      rownames(df) = c("Information rate",
                       "Efficacy boundary (Z)",
                       "Futility boundary (Z)",
                       "Cumulative rejection",
                       "Cumulative futility",
                       "Cumulative alpha spent",
                       "Number of subjects",
                       "Efficacy boundary (risk ratio score)",
                       "Futility boundary (risk ratio score)",
                       "Efficacy boundary (p)",
                       "Futility boundary (p)",
                       "Information",
                       "Cumulative rejection under H0",
                       "Cumulative futility under H0")

    } else {
      df = t(b[,c(1,2,4,6,7,8,10,12)])
      rownames(df) = c("Information rate",
                       "Efficacy boundary (Z)",
                       "Cumulative rejection",
                       "Cumulative alpha spent",
                       "Number of subjects",
                       "Efficacy boundary (risk ratio score)",
                       "Efficacy boundary (p)",
                       "Information")
    }

    colnames(df) <- paste("Stage", seq_len(ncol(df)), sep=" ")
  } else {
    b <- s[, c("efficacyBounds", "efficacyRiskRatioScore", "efficacyP")]

    # format number of digits after decimal for each column
    j3 <- c(1,2)
    j4 <- 3

    b[j3] <- lapply(b[j3], formatC, format = "f", digits = 3)
    b[j4] <- lapply(b[j4], formatC, format = "f", digits = 4)

    df = t(b)

    rownames(df) = c("Efficacy boundary (Z)",
                     "Efficacy boundary (risk ratio score)",
                     "Efficacy boundary (p)")
    colnames(df) <- NA
  }

  print(df1, ..., na.print = "" , quote = FALSE )
  print(df, ..., na.print = "" , quote = FALSE )
  invisible(x)
}


#' @title Print Power and Sample Size Results for Two-Sample Odds Ratio
#' @description Prints the summary statistics from power calculation of
#' two-sample odds ratio.
#'
#' @param x The designOddsRatio object to print.
#' @param ... Ensures that all arguments starting from "..." are named.
#'
#' @return A tabular printout of the design elements.
#'
#' @keywords internal
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @export
print.designOddsRatio <- function(x, ...) {
  a = x$overallResults
  s = x$byStageResults
  k = a$kMax

  if (k>1) {
    str1 = paste0("Group-sequential design with ", k, " stages")
  } else {
    str1 = "Fixed design"
  }

  str1 <- paste(str1, "for two-sample odds ratio")

  str2 <- paste0("Odds ratio under H0: ", round(a$oddsRatioH0, 3), ", ",
                 "proportion on treatment: ", round(a$pi1, 3), ", ",
                 "proportion on control: ", round(a$pi2, 3))

  str3 <- paste0("Overall power: ",
                 round(a$overallReject, 3), ", ",
                 "overall alpha (1-sided): ",
                 round(a$alpha, 4))

  if (x$settings$typeBetaSpending != 'none' ||
      (k > 1 && any(x$byStageResults$futilityBounds[1:(k-1)] > -6))) {
    str3 <- paste0(str3, ", ",
                   "attained alpha: ", round(a$attainedAlpha, 4))
  }

  str4 <- paste0("Drift parameter: ", round(a$drift, 3), ", ",
                 "inflation factor: ", round(a$inflationFactor, 3))

  if (k>1) {
    str5 <- paste0("Maximum information: ",
                   round(a$information, 2), ", ",
                   "expected under H1: ",
                   round(a$expectedInformationH1, 2), ", ",
                   "expected under H0: ",
                   round(a$expectedInformationH0, 2))

    str6 <- paste0("Maximum # subjects: ",
                   round(a$numberOfSubjects, 1), ", ",
                   "expected under H1: ",
                   round(a$expectedNumberOfSubjectsH1, 1), ", ",
                   "expected under H0: ",
                   round(a$expectedNumberOfSubjectsH0, 1))

  } else {
    str5 <- paste0("Information: ",
                   round(a$information, 2))

    str6 <- paste0("Number of subjects: ",
                   round(a$numberOfSubjects, 1))
  }

  str7 <- paste0("Allocation ratio: ",
                 round(x$settings$allocationRatioPlanned, 3), ", ",
                 "variance of standardized test statistic: ",
                 ifelse(x$settings$nullVariance, "under H0", "under H1"))

  if (k > 1) {
    asf = tolower(x$settings$typeAlphaSpending)
    asfpar = round(x$settings$parameterAlphaSpending, 3)
    asfuser = round(x$settings$userAlphaSpending, 4)

    bsf = tolower(x$settings$typeBetaSpending)
    bsfpar = round(x$settings$parameterBetaSpending, 3)
    bsfuser = round(x$settings$userBetaSpending, 4)

    if (asf == "of") {
      str8 = paste0("Alpha spending: O'Brien-Fleming")
    } else if (asf == "p") {
      str8 = paste0("Alpha spending: Pocock")
    } else if (asf == "wt") {
      str8 = paste0("Alpha spending: Wang-Tsiatis(Delta = ", asfpar, ")")
    } else if (asf == "sfof") {
      str8 = paste0("Alpha spending: Lan-DeMets O'Brien-Fleming")
    } else if (asf == "sfp") {
      str8 = paste0("Alpha spending: Lan-DeMets Pocock")
    } else if (asf == "sfkd") {
      str8 = paste0("Alpha spending: KD(rho = ", asfpar, ")")
    } else if (asf == "sfhsd") {
      str8 = paste0("Alpha spending: HSD(gamma = ", asfpar, ")")
    } else if (asf == "user") {
      str8 = paste0("Alpha spending: User defined(",
                    paste(asfuser, collapse = ","), ")")
    } else {
      str8 = "Alpha spending: None"
    }

    if (bsf == "of") {
      str9 = paste0("beta spending: O'Brien-Fleming")
    } else if (bsf == "p") {
      str9 = paste0("beta spending: Pocock")
    } else if (bsf == "wt") {
      str9 = paste0("beta spending: Wang-Tsiatis(Delta = ", bsfpar, ")")
    } else if (bsf == "sfof") {
      str9 = paste0("beta spending: Lan-DeMets O'Brien-Fleming")
    } else if (bsf == "sfp") {
      str9 = paste0("beta spending: Lan-DeMets Pocock")
    } else if (bsf == "sfkd") {
      str9 = paste0("beta spending: KD(rho = ", bsfpar, ")")
    } else if (bsf == "sfhsd") {
      str9 = paste0("beta spending: HSD(gamma = ", bsfpar, ")")
    } else if (bsf == "user") {
      str9 = paste0("beta spending: User defined(",
                    paste(bsfuser, collapse = ","), ")")
    } else {
      str9 = "beta spending: None"
    }

    if (!any(is.na(x$settings$spendingTime)) &&
        !all(x$settings$spendingTime == s$informationRates)) {
      str10 = paste0("Spending time: ",
                     paste(round(x$settings$spendingTime, 3), collapse = ","))
      df1 = data.frame(x = rep("", 10))
      colnames(df1) = NULL
      rownames(df1) = c(str1, str2, str3, str4, str5, str6, str7,
                        paste(str8, str9, sep = ", "), str10, "")
    } else {
      df1 = data.frame(x = rep("", 9))
      colnames(df1) = NULL
      rownames(df1) = c(str1, str2, str3, str4, str5, str6, str7,
                        paste(str8, str9, sep = ", "), "")
    }
  } else {
    df1 = data.frame(x = rep("", 8))
    colnames(df1) = NULL
    rownames(df1) = c(str1, str2, str3, str4, str5, str6, str7, "")
  }

  if (k>1) {
    b <- s[, c("informationRates", "efficacyBounds", "futilityBounds",
               "cumulativeRejection", "cumulativeFutility",
               "cumulativeAlphaSpent", "numberOfSubjects",
               "efficacyOddsRatio", "futilityOddsRatio",
               "efficacyP", "futilityP", "information",
               "cumulativeRejectionH0", "cumulativeFutilityH0")]

    # format number of digits after decimal for each column
    j1 <- 7
    j2 <- 12
    j3 <- c(1,2,3,4,8,9,14)
    j4 <- c(5,6,10,11,13)

    b[j1] <- lapply(b[j1], formatC, format = "f", digits = 1)
    b[j2] <- lapply(b[j2], formatC, format = "f", digits = 2)
    b[j3] <- lapply(b[j3], formatC, format = "f", digits = 3)
    b[j4] <- lapply(b[j4], formatC, format = "f", digits = 4)

    if (x$settings$typeBetaSpending != 'none' ||
        (k > 1 && any(x$byStageResults$futilityBounds[1:(k-1)] > -6))) {
      df = t(b)
      rownames(df) = c("Information rate",
                       "Efficacy boundary (Z)",
                       "Futility boundary (Z)",
                       "Cumulative rejection",
                       "Cumulative futility",
                       "Cumulative alpha spent",
                       "Number of subjects",
                       "Efficacy boundary (odds ratio)",
                       "Futility boundary (odds ratio)",
                       "Efficacy boundary (p)",
                       "Futility boundary (p)",
                       "Information",
                       "Cumulative rejection under H0",
                       "Cumulative futility under H0")

    } else {
      df = t(b[,c(1,2,4,6,7,8,10,12)])
      rownames(df) = c("Information rate",
                       "Efficacy boundary (Z)",
                       "Cumulative rejection",
                       "Cumulative alpha spent",
                       "Number of subjects",
                       "Efficacy boundary (odds ratio)",
                       "Efficacy boundary (p)",
                       "Information")
    }

    colnames(df) <- paste("Stage", seq_len(ncol(df)), sep=" ")
  } else {
    b <- s[, c("efficacyBounds", "efficacyOddsRatio", "efficacyP")]

    # format number of digits after decimal for each column
    j3 <- c(1,2)
    j4 <- 3

    b[j3] <- lapply(b[j3], formatC, format = "f", digits = 3)
    b[j4] <- lapply(b[j4], formatC, format = "f", digits = 4)

    df = t(b)

    rownames(df) = c("Efficacy boundary (Z)",
                     "Efficacy boundary (odds ratio)",
                     "Efficacy boundary (p)")
    colnames(df) <- NA
  }

  print(df1, ..., na.print = "" , quote = FALSE )
  print(df, ..., na.print = "" , quote = FALSE )
  invisible(x)
}


#' @title Print Power and Sample Size Results for Equivalence in
#' Two-Sample Risk Difference
#' @description Prints the summary statistics from power calculation of
#' equivalence in two-sample risk difference.
#'
#' @param x The designRiskDiffEquiv object to print.
#' @param ... Ensures that all arguments starting from "..." are named.
#'
#' @return A tabular printout of the design elements.
#'
#' @keywords internal
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @export
print.designRiskDiffEquiv <- function(x, ...) {
  a = x$overallResults
  s = x$byStageResults
  k = a$kMax

  if (k>1) {
    str1 = paste0("Group-sequential design with ", k, " stages")
  } else {
    str1 = "Fixed design"
  }

  str1 <- paste(str1, "for equivalence in two-sample risk difference")

  str2 <- paste0("Lower limit for risk difference: ",
                 round(a$riskDiffLower, 3), ", ",
                 "upper limit for risk difference: ",
                 round(a$riskDiffUpper, 3))

  str3 <- paste0("Proportion on treatment: ", round(a$pi1, 3), ", ",
                 "proportion on control: ", round(a$pi2, 3), ", ",
                 "risk difference: ", round(a$riskDiff, 3))

  str4 <- paste0("Overall power: ",
                 round(a$overallReject, 3), ", ",
                 "overall alpha: ",
                 round(a$alpha, 4), ", ",
                 "attained under H10: ",
                 round(a$attainedAlphaH10, 4), ", ",
                 "under H20: ",
                 round(a$attainedAlphaH20, 4))

  if (k>1) {
    str5 <- paste0("Max information: ",
                   round(a$information, 2), ", ",
                   "expected under H1: ",
                   round(a$expectedInformationH1, 2), ", ",
                   "under H10: ",
                   round(a$expectedInformationH10, 2), ", ",
                   "under H20: ",
                   round(a$expectedInformationH20, 2))

    str6 <- paste0("Max # subjects: ",
                   round(a$numberOfSubjects, 1), ", ",
                   "expected under H1: ",
                   round(a$expectedNumberOfSubjectsH1, 1), ", ",
                   "under H10: ",
                   round(a$expectedNumberOfSubjectsH10, 1), ", ",
                   "under H20: ",
                   round(a$expectedNumberOfSubjectsH20, 1))

  } else {
    str5 <- paste0("Information: ",
                   round(a$information, 2))

    str6 <- paste0("Number of subjects: ",
                   round(a$numberOfSubjects, 1))
  }

  str7 <- paste0("Allocation ratio: ",
                 round(x$settings$allocationRatioPlanned, 3), ", ",
                 "variance of standardized test statistic: ",
                 ifelse(x$settings$nullVariance, "under H0", "under H1"))

  if (k > 1) {
    asf = tolower(x$settings$typeAlphaSpending)
    asfpar = round(x$settings$parameterAlphaSpending, 3)
    asfuser = round(x$settings$userAlphaSpending, 4)

    if (asf == "of") {
      str8 = paste0("Alpha spending: O'Brien-Fleming")
    } else if (asf == "p") {
      str8 = paste0("Alpha spending: Pocock")
    } else if (asf == "wt") {
      str8 = paste0("Alpha spending: Wang-Tsiatis(Delta = ", asfpar, ")")
    } else if (asf == "sfof") {
      str8 = paste0("Alpha spending: Lan-DeMets O'Brien-Fleming")
    } else if (asf == "sfp") {
      str8 = paste0("Alpha spending: Lan-DeMets Pocock")
    } else if (asf == "sfkd") {
      str8 = paste0("Alpha spending: KD(rho = ", asfpar, ")")
    } else if (asf == "sfhsd") {
      str8 = paste0("Alpha spending: HSD(gamma = ", asfpar, ")")
    } else if (asf == "user") {
      str8 = paste0("Alpha spending: User defined(",
                    paste(asfuser, collapse = ","), ")")
    } else {
      str8 = "Alpha spending: None"
    }

    if (!any(is.na(x$settings$spendingTime)) &&
        !all(x$settings$spendingTime == s$informationRates)) {
      str9 = paste0("Spending time: ",
                    paste(round(x$settings$spendingTime, 3), collapse = ","))
      df1 = data.frame(x = rep("", 10))
      colnames(df1) = NULL
      rownames(df1) = c(str1, str2, str3, str4, str5, str6, str7,
                        str8, str9, "")
    } else {
      df1 = data.frame(x = rep("", 9))
      colnames(df1) = NULL
      rownames(df1) = c(str1, str2, str3, str4, str5, str6, str7, str8, "")
    }
  } else {
    df1 = data.frame(x = rep("", 8))
    colnames(df1) = NULL
    rownames(df1) = c(str1, str2, str3, str4, str5, str6, str7, "")
  }

  if (k>1) {
    b <- s[, c("informationRates", "efficacyBounds",
               "cumulativeRejection", "cumulativeAlphaSpent",
               "cumulativeAttainedAlphaH10",
               "cumulativeAttainedAlphaH20", "numberOfSubjects",
               "efficacyRiskDiffLower", "efficacyRiskDiffUpper",
               "efficacyP", "information")]

    # format number of digits after decimal for each column
    j1 <- 7
    j2 <- 11
    j3 <- c(1,2,3,8,9)
    j4 <- c(4,5,6,10)

    b[j1] <- lapply(b[j1], formatC, format = "f", digits = 1)
    b[j2] <- lapply(b[j2], formatC, format = "f", digits = 2)
    b[j3] <- lapply(b[j3], formatC, format = "f", digits = 3)
    b[j4] <- lapply(b[j4], formatC, format = "f", digits = 4)

    df = t(b)
    rownames(df) = c("Information rate",
                     "Boundary for each 1-sided test (Z)",
                     "Cumulative rejection",
                     "Cumulative alpha for each 1-sided test",
                     "Cumulative alpha attained under H10",
                     "Cumulative alpha attained under H20",
                     "Number of subjects",
                     "Boundary for lower limit (risk diff)",
                     "Boundary for upper limit (risk diff)",
                     "Boundary for each 1-sided test (p)",
                     "Information")

    colnames(df) <- paste("Stage", seq_len(ncol(df)), sep=" ")
  } else {
    b <- s[, c("efficacyBounds", "efficacyRiskDiffLower",
               "efficacyRiskDiffUpper", "efficacyP")]

    # format number of digits after decimal for each column
    j3 <- c(1,2,3)
    j4 <- 4

    b[j3] <- lapply(b[j3], formatC, format = "f", digits = 3)
    b[j4] <- lapply(b[j4], formatC, format = "f", digits = 4)

    df = t(b)

    rownames(df) = c("Boundary for each 1-sided test (Z)",
                     "Boundary for lower limit (risk diff)",
                     "Boundary for upper limit (risk diff)",
                     "Boundary for each 1-sided test (p)")
    colnames(df) <- NA
  }

  print(df1, ..., na.print = "" , quote = FALSE )
  print(df, ..., na.print = "" , quote = FALSE )
  invisible(x)
}


#' @title Print Power and Sample Size Results for Equivalence in
#' Two-Sample Risk Ratio
#' @description Prints the summary statistics from power calculation of
#' equivalence in two-sample risk ratio.
#'
#' @param x The designRiskRatioEquiv object to print.
#' @param ... Ensures that all arguments starting from "..." are named.
#'
#' @return A tabular printout of the design elements.
#'
#' @keywords internal
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @export
print.designRiskRatioEquiv <- function(x, ...) {
  a = x$overallResults
  s = x$byStageResults
  k = a$kMax

  if (k>1) {
    str1 = paste0("Group-sequential design with ", k, " stages")
  } else {
    str1 = "Fixed design"
  }

  str1 <- paste(str1, "for equivalence in two-sample risk ratio")

  str2 <- paste0("Lower limit for risk ratio: ",
                 round(a$riskRatioLower, 3), ", ",
                 "upper limit for risk ratio: ",
                 round(a$riskRatioUpper, 3))

  str3 <- paste0("Proportion on treatment: ", round(a$pi1, 3), ", ",
                 "proportion on control: ", round(a$pi2, 3), ", ",
                 "risk ratio: ", round(a$riskRatio, 3))

  str4 <- paste0("Overall power: ",
                 round(a$overallReject, 3), ", ",
                 "overall alpha: ",
                 round(a$alpha, 4), ", ",
                 "attained under H10: ",
                 round(a$attainedAlphaH10, 4), ", ",
                 "under H20: ",
                 round(a$attainedAlphaH20, 4))

  if (k>1) {
    str5 <- paste0("Max information: ",
                   round(a$information, 2), ", ",
                   "expected under H1: ",
                   round(a$expectedInformationH1, 2), ", ",
                   "under H10: ",
                   round(a$expectedInformationH10, 2), ", ",
                   "under H20: ",
                   round(a$expectedInformationH20, 2))

    str6 <- paste0("Max # subjects: ",
                   round(a$numberOfSubjects, 1), ", ",
                   "expected under H1: ",
                   round(a$expectedNumberOfSubjectsH1, 1), ", ",
                   "under H10: ",
                   round(a$expectedNumberOfSubjectsH10, 1), ", ",
                   "under H20: ",
                   round(a$expectedNumberOfSubjectsH20, 1))

  } else {
    str5 <- paste0("Information: ",
                   round(a$information, 2))

    str6 <- paste0("Number of subjects: ",
                   round(a$numberOfSubjects, 1))
  }

  str7 <- paste0("Allocation ratio: ",
                 round(x$settings$allocationRatioPlanned, 3), ", ",
                 "variance of standardized test statistic: ",
                 ifelse(x$settings$nullVariance, "under H0", "under H1"))

  if (k > 1) {
    asf = tolower(x$settings$typeAlphaSpending)
    asfpar = round(x$settings$parameterAlphaSpending, 3)
    asfuser = round(x$settings$userAlphaSpending, 4)

    if (asf == "of") {
      str8 = paste0("Alpha spending: O'Brien-Fleming")
    } else if (asf == "p") {
      str8 = paste0("Alpha spending: Pocock")
    } else if (asf == "wt") {
      str8 = paste0("Alpha spending: Wang-Tsiatis(Delta = ", asfpar, ")")
    } else if (asf == "sfof") {
      str8 = paste0("Alpha spending: Lan-DeMets O'Brien-Fleming")
    } else if (asf == "sfp") {
      str8 = paste0("Alpha spending: Lan-DeMets Pocock")
    } else if (asf == "sfkd") {
      str8 = paste0("Alpha spending: KD(rho = ", asfpar, ")")
    } else if (asf == "sfhsd") {
      str8 = paste0("Alpha spending: HSD(gamma = ", asfpar, ")")
    } else if (asf == "user") {
      str8 = paste0("Alpha spending: User defined(",
                    paste(asfuser, collapse = ","), ")")
    } else {
      str8 = "Alpha spending: None"
    }

    if (!any(is.na(x$settings$spendingTime)) &&
        !all(x$settings$spendingTime == s$informationRates)) {
      str9 = paste0("Spending time: ",
                    paste(round(x$settings$spendingTime, 3), collapse = ","))
      df1 = data.frame(x = rep("", 10))
      colnames(df1) = NULL
      rownames(df1) = c(str1, str2, str3, str4, str5, str6, str7,
                        str8, str9, "")
    } else {
      df1 = data.frame(x = rep("", 9))
      colnames(df1) = NULL
      rownames(df1) = c(str1, str2, str3, str4, str5, str6, str7, str8, "")
    }
  } else {
    df1 = data.frame(x = rep("", 8))
    colnames(df1) = NULL
    rownames(df1) = c(str1, str2, str3, str4, str5, str6, str7, "")
  }

  if (k>1) {
    b <- s[, c("informationRates", "efficacyBounds",
               "cumulativeRejection", "cumulativeAlphaSpent",
               "cumulativeAttainedAlphaH10",
               "cumulativeAttainedAlphaH20", "numberOfSubjects",
               "efficacyRiskRatioLower", "efficacyRiskRatioUpper",
               "efficacyP", "information")]

    # format number of digits after decimal for each column
    j1 <- 7
    j2 <- 11
    j3 <- c(1,2,3,8,9)
    j4 <- c(4,5,6,10)

    b[j1] <- lapply(b[j1], formatC, format = "f", digits = 1)
    b[j2] <- lapply(b[j2], formatC, format = "f", digits = 2)
    b[j3] <- lapply(b[j3], formatC, format = "f", digits = 3)
    b[j4] <- lapply(b[j4], formatC, format = "f", digits = 4)

    df = t(b)
    rownames(df) = c("Information rate",
                     "Boundary for each 1-sided test (Z)",
                     "Cumulative rejection",
                     "Cumulative alpha for each 1-sided test",
                     "Cumulative alpha attained under H10",
                     "Cumulative alpha attained under H20",
                     "Number of subjects",
                     "Boundary for lower limit (risk ratio)",
                     "Boundary for upper limit (risk ratio)",
                     "Boundary for each 1-sided test (p)",
                     "Information")

    colnames(df) <- paste("Stage", seq_len(ncol(df)), sep=" ")
  } else {
    b <- s[, c("efficacyBounds", "efficacyRiskRatioLower",
               "efficacyRiskRatioUpper", "efficacyP")]

    # format number of digits after decimal for each column
    j3 <- c(1,2,3)
    j4 <- 4

    b[j3] <- lapply(b[j3], formatC, format = "f", digits = 3)
    b[j4] <- lapply(b[j4], formatC, format = "f", digits = 4)

    df = t(b)

    rownames(df) = c("Boundary for each 1-sided test (Z)",
                     "Boundary for lower limit (risk ratio)",
                     "Boundary for upper limit (risk ratio)",
                     "Boundary for each 1-sided test (p)")
    colnames(df) <- NA
  }

  print(df1, ..., na.print = "" , quote = FALSE )
  print(df, ..., na.print = "" , quote = FALSE )
  invisible(x)
}


#' @title Print Power and Sample Size Results for Equivalence in
#' Two-Sample Odds Ratio
#' @description Prints the summary statistics from power calculation of
#' equivalence in two-sample odds ratio.
#'
#' @param x The designOddsRatioEquiv object to print.
#' @param ... Ensures that all arguments starting from "..." are named.
#'
#' @return A tabular printout of the design elements.
#'
#' @keywords internal
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @export
print.designOddsRatioEquiv <- function(x, ...) {
  a = x$overallResults
  s = x$byStageResults
  k = a$kMax

  if (k>1) {
    str1 = paste0("Group-sequential design with ", k, " stages")
  } else {
    str1 = "Fixed design"
  }

  str1 <- paste(str1, "for equivalence in two-sample odds ratio")

  str2 <- paste0("Lower limit for odds ratio: ",
                 round(a$oddsRatioLower, 3), ", ",
                 "upper limit for odds ratio: ",
                 round(a$oddsRatioUpper, 3))

  str3 <- paste0("Proportion on treatment: ", round(a$pi1, 3), ", ",
                 "proportion on control: ", round(a$pi2, 3), ", ",
                 "odds ratio: ", round(a$oddsRatio, 3))

  str4 <- paste0("Overall power: ",
                 round(a$overallReject, 3), ", ",
                 "overall alpha (1-sided): ",
                 round(a$alpha, 4), ", ",
                 "attained under H10: ",
                 round(a$attainedAlphaH10, 4), ", ",
                 "under H20: ",
                 round(a$attainedAlphaH20, 4))

  if (k>1) {
    str5 <- paste0("Max information: ",
                   round(a$information, 2), ", ",
                   "expected under H1: ",
                   round(a$expectedInformationH1, 2), ", ",
                   "under H10: ",
                   round(a$expectedInformationH10, 2), ", ",
                   "under H20: ",
                   round(a$expectedInformationH20, 2))

    str6 <- paste0("Max # subjects: ",
                   round(a$numberOfSubjects, 1), ", ",
                   "expected under H1: ",
                   round(a$expectedNumberOfSubjectsH1, 1), ", ",
                   "under H10: ",
                   round(a$expectedNumberOfSubjectsH10, 1), ", ",
                   "under H20: ",
                   round(a$expectedNumberOfSubjectsH20, 1))

  } else {
    str5 <- paste0("Information: ",
                   round(a$information, 2))

    str6 <- paste0("Number of subjects: ",
                   round(a$numberOfSubjects, 1))
  }

  str7 <- paste0("Allocation ratio: ",
                 round(x$settings$allocationRatioPlanned, 3), ", ",
                 "variance of standardized test statistic: ",
                 ifelse(x$settings$nullVariance, "under H0", "under H1"))

  if (k > 1) {
    asf = tolower(x$settings$typeAlphaSpending)
    asfpar = round(x$settings$parameterAlphaSpending, 3)
    asfuser = round(x$settings$userAlphaSpending, 4)

    if (asf == "of") {
      str8 = paste0("Alpha spending: O'Brien-Fleming")
    } else if (asf == "p") {
      str8 = paste0("Alpha spending: Pocock")
    } else if (asf == "wt") {
      str8 = paste0("Alpha spending: Wang-Tsiatis(Delta = ", asfpar, ")")
    } else if (asf == "sfof") {
      str8 = paste0("Alpha spending: Lan-DeMets O'Brien-Fleming")
    } else if (asf == "sfp") {
      str8 = paste0("Alpha spending: Lan-DeMets Pocock")
    } else if (asf == "sfkd") {
      str8 = paste0("Alpha spending: KD(rho = ", asfpar, ")")
    } else if (asf == "sfhsd") {
      str8 = paste0("Alpha spending: HSD(gamma = ", asfpar, ")")
    } else if (asf == "user") {
      str8 = paste0("Alpha spending: User defined(",
                    paste(asfuser, collapse = ","), ")")
    } else {
      str8 = "Alpha spending: None"
    }

    if (!any(is.na(x$settings$spendingTime)) &&
        !all(x$settings$spendingTime == s$informationRates)) {
      str9 = paste0("Spending time: ",
                    paste(round(x$settings$spendingTime, 3), collapse = ","))
      df1 = data.frame(x = rep("", 10))
      colnames(df1) = NULL
      rownames(df1) = c(str1, str2, str3, str4, str5, str6, str7,
                        str8, str9, "")
    } else {
      df1 = data.frame(x = rep("", 9))
      colnames(df1) = NULL
      rownames(df1) = c(str1, str2, str3, str4, str5, str6, str7, str8, "")
    }
  } else {
    df1 = data.frame(x = rep("", 8))
    colnames(df1) = NULL
    rownames(df1) = c(str1, str2, str3, str4, str5, str6, str7, "")
  }

  if (k>1) {
    b <- s[, c("informationRates", "efficacyBounds",
               "cumulativeRejection", "cumulativeAlphaSpent",
               "cumulativeAttainedAlphaH10",
               "cumulativeAttainedAlphaH20", "numberOfSubjects",
               "efficacyOddsRatioLower", "efficacyOddsRatioUpper",
               "efficacyP", "information")]

    # format number of digits after decimal for each column
    j1 <- 7
    j2 <- 11
    j3 <- c(1,2,3,8,9)
    j4 <- c(4,5,6,10)

    b[j1] <- lapply(b[j1], formatC, format = "f", digits = 1)
    b[j2] <- lapply(b[j2], formatC, format = "f", digits = 2)
    b[j3] <- lapply(b[j3], formatC, format = "f", digits = 3)
    b[j4] <- lapply(b[j4], formatC, format = "f", digits = 4)

    df = t(b)
    rownames(df) = c("Information rate",
                     "Boundary for each 1-sided test (Z)",
                     "Cumulative rejection",
                     "Cumulative alpha for each 1-sided test",
                     "Cumulative alpha attained under H10",
                     "Cumulative alpha attained under H20",
                     "Number of subjects",
                     "Boundary for lower limit (odds ratio)",
                     "Boundary for upper limit (odds ratio)",
                     "Boundary for each 1-sided test (p)",
                     "Information")

    colnames(df) <- paste("Stage", seq_len(ncol(df)), sep=" ")
  } else {
    b <- s[, c("efficacyBounds", "efficacyOddsRatioLower",
               "efficacyOddsRatioUpper", "efficacyP")]

    # format number of digits after decimal for each column
    j3 <- c(1,2,3)
    j4 <- 4

    b[j3] <- lapply(b[j3], formatC, format = "f", digits = 3)
    b[j4] <- lapply(b[j4], formatC, format = "f", digits = 4)

    df = t(b)

    rownames(df) = c("Boundary for each 1-sided test (Z)",
                     "Boundary for lower limit (odds ratio)",
                     "Boundary for upper limit (odds ratio)",
                     "Boundary for each 1-sided test (p)")
    colnames(df) <- NA
  }

  print(df1, ..., na.print = "" , quote = FALSE )
  print(df, ..., na.print = "" , quote = FALSE )
  invisible(x)
}


#' @title Print the Decision Table for the mTPI-2 Design
#' @description Prints the decision table of the modified toxicity
#' probability-2 (mTPI-2) design for MTD finding.
#'
#' @param x The mTPI2Table object to print.
#' @param ... Ensures that all arguments starting from "..." are named.
#'
#' @return A tabular printout of the design elements.
#'
#' @keywords internal
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @export
print.mTPI2Table <- function(x, ...) {

  str1 <- "Trial monitoring table for mTPI-2 design"
  df1 = data.frame(x = rep("", 2))
  colnames(df1) = NULL
  rownames(df1) = c(str1, "")

  df2 <- as.data.frame(x$decisionMatrix)

  str3 <- "Rows represent number of toxicities"
  str4 <- "Columns represent number of patients treated at current dose"
  str5 <- "E = Escalate to the next higher dose"
  str6 <- "S = Stay at the current dose"
  str7 <- "D = De-escalate to the next lower dose"
  str8 <- "DU = The current dose is unacceptably toxic"
  str9 <- paste0("Target toxicity: ", round(x$settings$pT, 3), ", ",
                 "epsilon1: ", round(x$settings$epsilon1, 3), ", ",
                 "epsilon2: ", round(x$settings$epsilon2, 3))

  df3 = data.frame(x = rep("", 7))
  colnames(df3) = NULL
  rownames(df3) = c(str3, str4, str5, str6, str7, str8, str9)

  print(df1, ..., na.print = "" , quote = FALSE )
  print(df2, ..., na.print = "" , quote = FALSE )
  print(df3, ..., na.print = "" , quote = FALSE )
  invisible(x)
}


#' @title Print the Decision Table for the BOIN Design
#' @description Prints the decision table of the Bayesian optimal interval
#' design for MTD finding.
#'
#' @param x The BOINTable object to print.
#' @param ... Ensures that all arguments starting from "..." are named.
#'
#' @return A tabular printout of the design elements.
#'
#' @keywords internal
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @export
print.BOINTable <- function(x, ...) {

  str1 <- "Trial monitoring table for BOIN design"
  df1 = data.frame(x = rep("", 2))
  colnames(df1) = NULL
  rownames(df1) = c(str1, "")

  df2 <- as.data.frame(x$decisionMatrix)

  str3 <- "Rows represent number of toxicities"
  str4 <- "Columns represent number of patients treated at current dose"
  str5 <- "E = Escalate to the next higher dose"
  str6 <- "S = Stay at the current dose"
  str7 <- "D = De-escalate to the next lower dose"
  str8 <- "DU = The current dose is unacceptably toxic"
  str9 <- paste0("Target toxicity: ", round(x$settings$pT, 3), ", ",
                 "phi1: ", round(x$settings$phi1, 3), ", ",
                 "phi2: ", round(x$settings$phi2, 3), ", ",
                 "lambda1: ", round(x$settings$lambda1, 3), ", ",
                 "lambda2: ", round(x$settings$lambda2, 3))

  df3 = data.frame(x = rep("", 7))
  colnames(df3) = NULL
  rownames(df3) = c(str3, str4, str5, str6, str7, str8, str9)

  print(df1, ..., na.print = "" , quote = FALSE )
  print(df2, ..., na.print = "" , quote = FALSE )
  print(df3, ..., na.print = "" , quote = FALSE )
  invisible(x)
}


#' @title Print the Point Estimate and Confidence Interval
#' @description Prints the point estimate and confidence interval.
#'
#' @param x The estimateCI object to print.
#' @param ... Ensures that all arguments starting from "..." are named.
#'
#' @return A tabular printout of the design elements.
#'
#' @keywords internal
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @export
print.estimateCI <- function(x, ...) {
  df1 <- x$estimates
  rownames(df1) = NULL
  print(df1, ..., na.print = "" , quote = FALSE )
  invisible(x)
}


#' @title Print Power and Sample Size Results for One-Sample Multinomial
#' Response
#' @description Prints the power and sample size for one-sample multinomial
#' response.
#'
#' @param x The designOneMultinom object to print.
#' @param ... Ensures that all arguments starting from "..." are named.
#'
#' @return A tabular printout of the design elements.
#'
#' @keywords internal
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @export
print.designOneMultinom <- function(x, ...) {
  df1 <- data.frame(alpha = x$alpha, power = x$power,
                    n = x$n, ncats = x$ncats,
                    effectsize = x$effectsize)
  rownames(df1) = NULL
  print(df1, ..., na.print = "" , quote = FALSE )
  invisible(x)
}


#' @title Print Power and Sample Size Results for Difference in Two-Sample
#' Multinomial Response
#' @description Prints the power and sample size for difference in
#' two-sample multinomial response.
#'
#' @param x The designTwoMultinom object to print.
#' @param ... Ensures that all arguments starting from "..." are named.
#'
#' @return A tabular printout of the design elements.
#'
#' @keywords internal
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @export
print.designTwoMultinom <- function(x, ...) {
  df1 <- data.frame(alpha = x$alpha, power = x$power,
                    n = x$n, ncats = x$ncats,
                    effectsize = x$effectsize)
  rownames(df1) = NULL
  print(df1, ..., na.print = "" , quote = FALSE )
  invisible(x)
}


#' @title Print Power and Sample Size Results for Wilcoxon Test for
#' Two-Sample Ordinal Response.
#' @description Prints the power and sample size for Wilcoxon test for
#' two-sample ordinal response.
#'
#' @param x The designTwoOrdinal object to print.
#' @param ... Ensures that all arguments starting from "..." are named.
#'
#' @return A tabular printout of the design elements.
#'
#' @keywords internal
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @export
print.designTwoOrdinal <- function(x, ...) {
  df1 <- data.frame(alpha = x$alpha, power = x$power,
                    n = x$n, ncats = x$ncats,
                    meanscore1 = x$meanscore1,
                    meanscore2 = x$meanscore2)
  rownames(df1) = NULL
  print(df1, ..., na.print = "" , quote = FALSE )
  invisible(x)
}


#' @title Print Power and Sample Size Results for Cochran-Armitage Trend Test
#' for Ordered Multi-Sample Binomial Response
#' @description Prints the power and sample size for Cochran-Armitage
#' trend test for ordered multi-sample binomial response.
#'
#' @param x The designOrderedBinom object to print.
#' @param ... Ensures that all arguments starting from "..." are named.
#'
#' @return A tabular printout of the design elements.
#'
#' @keywords internal
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @export
print.designOrderedBinom <- function(x, ...) {
  df1 <- data.frame(alpha = x$alpha, power = x$power,
                    n = x$n, ngroups = x$ngroups,
                    trendstat = x$trendstat)
  rownames(df1) = NULL
  print(df1, ..., na.print = "" , quote = FALSE )
  invisible(x)
}


#' @title Print Power and Sample Size Results for Unordered Multi-Sample
#' Binomial Response
#' @description Prints the power and sample size for unordered multi-sample
#' binomial response.
#'
#' @param x The designUnorderedBinom object to print.
#' @param ... Ensures that all arguments starting from "..." are named.
#'
#' @return A tabular printout of the design elements.
#'
#' @keywords internal
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @export
print.designUnorderedBinom <- function(x, ...) {
  df1 <- data.frame(alpha = x$alpha, power = x$power,
                    n = x$n, ngroups = x$ngroups,
                    effectsize = x$effectsize)
  rownames(df1) = NULL
  print(df1, ..., na.print = "" , quote = FALSE )
  invisible(x)
}


#' @title Print Power and Sample Size Results for Unordered Multi-Sample
#' Multinomial Response
#' @description Prints the power and sample size for unordered multi-sample
#' multinomial response.
#'
#' @param x The designUnorderedMultinom object to print.
#' @param ... Ensures that all arguments starting from "..." are named.
#'
#' @return A tabular printout of the design elements.
#'
#' @keywords internal
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @export
print.designUnorderedMultinom <- function(x, ...) {
  df1 <- data.frame(alpha = x$alpha, power = x$power,
                    n = x$n, ngroups = x$ngroups,
                    ncats = x$ncats, effectsize = x$effectsize)
  rownames(df1) = NULL
  print(df1, ..., na.print = "" , quote = FALSE )
  invisible(x)
}


#' @title Print Power and Sample Size Results for Logistic Regression
#' @description Prints the power and sample size for logistic regression.
#'
#' @param x The designLogistic object to print.
#' @param ... Ensures that all arguments starting from "..." are named.
#'
#' @return A tabular printout of the design elements.
#'
#' @keywords internal
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @export
print.designLogistic <- function(x, ...) {
  df1 <- data.frame(alpha = x$alpha, power = x$power,
                    n = x$n, ncovariates = x$ncovariates,
                    corr = x$corr,
                    responseprob = x$responseprob,
                    oddsratio = x$oddsratios[1],
                    effectsize = x$effectsize)
  rownames(df1) = NULL
  print(df1, ..., na.print = "" , quote = FALSE )
  invisible(x)
}


#' @title Print Power and Sample Size Results for Cohen's kappa.
#' @description Prints the power and sample size for Cohen's kappa.
#'
#' @param x The designAgreement object to print.
#' @param ... Ensures that all arguments starting from "..." are named.
#'
#' @return A tabular printout of the design elements.
#'
#' @keywords internal
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @export
print.designAgreement <- function(x, ...) {
  df1 <- data.frame(alpha = x$alpha, power = x$power,
                    n = x$n, kappaH0 = x$kappaH0,
                    kappa = x$kappa)
  rownames(df1) = NULL
  print(df1, ..., na.print = "" , quote = FALSE )
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
  a = x$overallResults
  s = x$byStageResults
  k = a$kMax

  if (k>1) {
    str1 = paste0("Group-sequential design with ", k, " stages")
  } else {
    str1 = "Fixed design"
  }

  str1 <- paste0(str1, " for difference in milestone survival")

  str2 <- paste0("Milestone: ", round(a$milestone, 3), ", ",
                 "survival difference under H0: ",
                 round(a$survDiffH0, 3))

  str3 <- paste0("Milestone survival on treatment: ",
                 round(a$surv1, 3), ", ",
                 "on control: ", round(a$surv2, 3))

  str4 <- paste0("Overall power: ",
                 round(a$overallReject, 3), ", ",
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
    asf = tolower(x$settings$typeAlphaSpending)
    asfpar = round(x$settings$parameterAlphaSpending, 3)
    asfuser = round(x$settings$userAlphaSpending, 4)

    bsf = tolower(x$settings$typeBetaSpending)
    bsfpar = round(x$settings$parameterBetaSpending, 3)

    if (asf == "of") {
      str11 = paste0("Alpha spending: O'Brien-Fleming")
    } else if (asf == "p") {
      str11 = paste0("Alpha spending: Pocock")
    } else if (asf == "wt") {
      str11 = paste0("Alpha spending: Wang-Tsiatis(Delta = ", asfpar, ")")
    } else if (asf == "sfof") {
      str11 = paste0("Alpha spending: Lan-DeMets O'Brien-Fleming")
    } else if (asf == "sfp") {
      str11 = paste0("Alpha spending: Lan-DeMets Pocock")
    } else if (asf == "sfkd") {
      str11 = paste0("Alpha spending: KD(rho = ", asfpar, ")")
    } else if (asf == "sfhsd") {
      str11 = paste0("Alpha spending: HSD(gamma = ", asfpar, ")")
    } else if (asf == "user") {
      str11 = paste0("Alpha spending: User defined(",
                     paste(asfuser, collapse = ","), ")")
    } else {
      str11 = "Alpha spending: None"
    }

    if (bsf == "of") {
      str12 = paste0("beta spending: O'Brien-Fleming")
    } else if (bsf == "p") {
      str12 = paste0("beta spending: Pocock")
    } else if (bsf == "wt") {
      str12 = paste0("beta spending: Wang-Tsiatis(Delta = ", bsfpar, ")")
    } else if (bsf == "sfof") {
      str12 = paste0("beta spending: Lan-DeMets O'Brien-Fleming")
    } else if (bsf == "sfp") {
      str12 = paste0("beta spending: Lan-DeMets Pocock")
    } else if (bsf == "sfkd") {
      str12 = paste0("beta spending: KD(rho = ", bsfpar, ")")
    } else if (bsf == "sfhsd") {
      str12 = paste0("beta spending: HSD(gamma = ", bsfpar, ")")
    } else if (bsf == "user") {
      str12 = paste0("beta spending: User defined")
    } else {
      str12 = "beta spending: None"
    }

    if (!any(is.na(x$settings$spendingTime)) &&
        !all(x$settings$spendingTime == s$informationRates)) {
      str13 = paste0("Spending time: ",
                     paste(round(x$settings$spendingTime, 3), collapse = ","))
      df1 = data.frame(x = rep("", 13))
      colnames(df1) = NULL
      rownames(df1) = c(str1, str2, str3, str4, str5, str6, str7,
                        str8, str9, str10, paste(str11, str12, sep = ", "),
                        str13, "")
    } else {
      df1 = data.frame(x = rep("", 12))
      colnames(df1) = NULL
      rownames(df1) = c(str1, str2, str3, str4, str5, str6, str7,
                        str8, str9, str10, paste(str11, str12, sep = ", "),
                        "")
    }
  } else {
    df1 = data.frame(x = rep("", 11))
    colnames(df1) = NULL
    rownames(df1) = c(str1, str2, str3, str4, str5, str6, str7,
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
        (k > 1 && any(x$byStageResults$futilityBounds[1:(k-1)] > -6))) {
      df = t(b)
      rownames(df) = c("Information rate",
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
      df = t(b[,c(1,2,4,6,7,8,9,10,11,12,14,16)])
      rownames(df) = c("Information rate",
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

    df = t(b)

    rownames(df) = c("Efficacy boundary (Z)",
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
  a = x$overallResults
  s = x$byStageResults
  k = a$kMax

  if (k>1) {
    str1 = paste0("Group-sequential design with ", k, " stages")
  } else {
    str1 = "Fixed design"
  }

  str1 <- paste0(str1, " for difference in restricted mean survival time")

  str2 <- paste0("Milestone: ", round(a$milestone, 3), ", ",
                 "restricted mean survival time difference under H0: ",
                 round(a$rmstDiffH0, 3))

  str3 <- paste0("Restricted mean survival time on treatment: ",
                 round(a$rmst1, 3), ", ",
                 "on control: ", round(a$rmst2, 3))

  str4 <- paste0("Overall power: ",
                 round(a$overallReject, 3), ", ",
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
    asf = tolower(x$settings$typeAlphaSpending)
    asfpar = round(x$settings$parameterAlphaSpending, 3)
    asfuser = round(x$settings$userAlphaSpending, 4)

    bsf = tolower(x$settings$typeBetaSpending)
    bsfpar = round(x$settings$parameterBetaSpending, 3)

    if (asf == "of") {
      str11 = paste0("Alpha spending: O'Brien-Fleming")
    } else if (asf == "p") {
      str11 = paste0("Alpha spending: Pocock")
    } else if (asf == "wt") {
      str11 = paste0("Alpha spending: Wang-Tsiatis(Delta = ", asfpar, ")")
    } else if (asf == "sfof") {
      str11 = paste0("Alpha spending: Lan-DeMets O'Brien-Fleming")
    } else if (asf == "sfp") {
      str11 = paste0("Alpha spending: Lan-DeMets Pocock")
    } else if (asf == "sfkd") {
      str11 = paste0("Alpha spending: KD(rho = ", asfpar, ")")
    } else if (asf == "sfhsd") {
      str11 = paste0("Alpha spending: HSD(gamma = ", asfpar, ")")
    } else if (asf == "user") {
      str11 = paste0("Alpha spending: User defined(",
                     paste(asfuser, collapse = ","), ")")
    } else {
      str11 = "Alpha spending: None"
    }

    if (bsf == "of") {
      str12 = paste0("beta spending: O'Brien-Fleming")
    } else if (bsf == "p") {
      str12 = paste0("beta spending: Pocock")
    } else if (bsf == "wt") {
      str12 = paste0("beta spending: Wang-Tsiatis(Delta = ", bsfpar, ")")
    } else if (bsf == "sfof") {
      str12 = paste0("beta spending: Lan-DeMets O'Brien-Fleming")
    } else if (bsf == "sfp") {
      str12 = paste0("beta spending: Lan-DeMets Pocock")
    } else if (bsf == "sfkd") {
      str12 = paste0("beta spending: KD(rho = ", bsfpar, ")")
    } else if (bsf == "sfhsd") {
      str12 = paste0("beta spending: HSD(gamma = ", bsfpar, ")")
    } else if (bsf == "user") {
      str12 = paste0("beta spending: User defined")
    } else {
      str12 = "beta spending: None"
    }

    if (!any(is.na(x$settings$spendingTime)) &&
        !all(x$settings$spendingTime == s$informationRates)) {
      str13 = paste0("Spending time: ",
                     paste(round(x$settings$spendingTime, 3), collapse = ","))
      df1 = data.frame(x = rep("", 13))
      colnames(df1) = NULL
      rownames(df1) = c(str1, str2, str3, str4, str5, str6, str7,
                        str8, str9, str10, paste(str11, str12, sep = ", "),
                        str13, "")
    } else {
      df1 = data.frame(x = rep("", 12))
      colnames(df1) = NULL
      rownames(df1) = c(str1, str2, str3, str4, str5, str6, str7,
                        str8, str9, str10, paste(str11, str12, sep = ", "),
                        "")
    }
  } else {
    df1 = data.frame(x = rep("", 11))
    colnames(df1) = NULL
    rownames(df1) = c(str1, str2, str3, str4, str5, str6, str7,
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
        (k > 1 && any(x$byStageResults$futilityBounds[1:(k-1)] > -6))) {
      df = t(b)
      rownames(df) = c("Information rate",
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
      df = t(b[,c(1,2,4,6,7,8,9,10,11,12,14,16)])
      rownames(df) = c("Information rate",
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

    df = t(b)

    rownames(df) = c("Efficacy boundary (Z)",
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
  a = x$overallResults
  s = x$byStageResults
  k = a$kMax

  if (k>1) {
    str1 = paste0("Group-sequential design with ", k, " stages")
  } else {
    str1 = "Fixed design"
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
                 round(a$overallReject, 3), ", ",
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
    asf = tolower(x$settings$typeAlphaSpending)
    asfpar = round(x$settings$parameterAlphaSpending, 3)
    asfuser = round(x$settings$userAlphaSpending, 4)

    if (asf == "of") {
      str11 = paste0("Alpha spending: O'Brien-Fleming")
    } else if (asf == "p") {
      str11 = paste0("Alpha spending: Pocock")
    } else if (asf == "wt") {
      str11 = paste0("Alpha spending: Wang-Tsiatis(Delta = ", asfpar, ")")
    } else if (asf == "sfof") {
      str11 = paste0("Alpha spending: Lan-DeMets O'Brien-Fleming")
    } else if (asf == "sfp") {
      str11 = paste0("Alpha spending: Lan-DeMets Pocock")
    } else if (asf == "sfkd") {
      str11 = paste0("Alpha spending: KD(rho = ", asfpar, ")")
    } else if (asf == "sfhsd") {
      str11 = paste0("Alpha spending: HSD(gamma = ", asfpar, ")")
    } else if (asf == "user") {
      str11 = paste0("Alpha spending: User defined(",
                     paste(asfuser, collapse = ","), ")")
    } else {
      str11 = "Alpha spending: None"
    }

    if (!any(is.na(x$settings$spendingTime)) &&
        !all(x$settings$spendingTime == s$informationRates)) {
      str12 = paste0("Spending time: ",
                     paste(round(x$settings$spendingTime, 3), collapse = ","))
      df1 = data.frame(x = rep("", 13))
      colnames(df1) = NULL
      rownames(df1) = c(str1, str2, str3, str4, str5, str6, str7,
                        str8, str9, str10, str11, str12, "")
    } else {
      df1 = data.frame(x = rep("", 12))
      colnames(df1) = NULL
      rownames(df1) = c(str1, str2, str3, str4, str5, str6, str7,
                        str8, str9, str10, str11, "")
    }
  } else {
    df1 = data.frame(x = rep("", 11))
    colnames(df1) = NULL
    rownames(df1) = c(str1, str2, str3, str4, str5, str6, str7,
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
    j3 <- c(1,2,3,12,13)
    j4 <- c(4,5,6,14)

    b[j1] <- lapply(b[j1], formatC, format = "f", digits = 1)
    b[j2] <- lapply(b[j2], formatC, format = "f", digits = 2)
    b[j3] <- lapply(b[j3], formatC, format = "f", digits = 3)
    b[j4] <- lapply(b[j4], formatC, format = "f", digits = 4)

    df = t(b)
    rownames(df) = c("Information rate",
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

    df = t(b)

    rownames(df) = c("Boundary for each 1-sided test (Z)",
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
  a = x$overallResults
  s = x$byStageResults
  k = a$kMax

  if (k>1) {
    str1 = paste0("Group-sequential design with ", k, " stages")
  } else {
    str1 = "Fixed design"
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
                 round(a$overallReject, 3), ", ",
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
    asf = tolower(x$settings$typeAlphaSpending)
    asfpar = round(x$settings$parameterAlphaSpending, 3)
    asfuser = round(x$settings$userAlphaSpending, 4)

    if (asf == "of") {
      str11 = paste0("Alpha spending: O'Brien-Fleming")
    } else if (asf == "p") {
      str11 = paste0("Alpha spending: Pocock")
    } else if (asf == "wt") {
      str11 = paste0("Alpha spending: Wang-Tsiatis(Delta = ", asfpar, ")")
    } else if (asf == "sfof") {
      str11 = paste0("Alpha spending: Lan-DeMets O'Brien-Fleming")
    } else if (asf == "sfp") {
      str11 = paste0("Alpha spending: Lan-DeMets Pocock")
    } else if (asf == "sfkd") {
      str11 = paste0("Alpha spending: KD(rho = ", asfpar, ")")
    } else if (asf == "sfhsd") {
      str11 = paste0("Alpha spending: HSD(gamma = ", asfpar, ")")
    } else if (asf == "user") {
      str11 = paste0("Alpha spending: User defined(",
                     paste(asfuser, collapse = ","), ")")
    } else {
      str11 = "Alpha spending: None"
    }

    if (!any(is.na(x$settings$spendingTime)) &&
        !all(x$settings$spendingTime == s$informationRates)) {
      str12 = paste0("Spending time: ",
                     paste(round(x$settings$spendingTime, 3), collapse = ","))
      df1 = data.frame(x = rep("", 13))
      colnames(df1) = NULL
      rownames(df1) = c(str1, str2, str3, str4, str5, str6, str7,
                        str8, str9, str10, str11, str12, "")
    } else {
      df1 = data.frame(x = rep("", 12))
      colnames(df1) = NULL
      rownames(df1) = c(str1, str2, str3, str4, str5, str6, str7,
                        str8, str9, str10, str11, "")
    }
  } else {
    df1 = data.frame(x = rep("", 11))
    colnames(df1) = NULL
    rownames(df1) = c(str1, str2, str3, str4, str5, str6, str7,
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
    j3 <- c(1,2,3,12,13)
    j4 <- c(4,5,6,14)

    b[j1] <- lapply(b[j1], formatC, format = "f", digits = 1)
    b[j2] <- lapply(b[j2], formatC, format = "f", digits = 2)
    b[j3] <- lapply(b[j3], formatC, format = "f", digits = 3)
    b[j4] <- lapply(b[j4], formatC, format = "f", digits = 4)

    df = t(b)
    rownames(df) = c("Information rate",
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

    df = t(b)

    rownames(df) = c("Boundary for each 1-sided test (Z)",
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
  a = x$overallResults
  s = x$byStageResults
  k = a$kMax

  if (k>1) {
    str1 = paste0("Group-sequential design with ", k, " stages")
  } else {
    str1 = "Fixed design"
  }

  str1 <- paste(str1, "for equivalence in hazard ratio")

  str2 <- paste0("Lower limit for hazard ratio: ",
                 round(a$hazardRatioLower, 3), ", ",
                 "upper limit for hazard ratio: ",
                 round(a$hazardRatioUpper, 3))

  str3 <- paste0("Overall power: ",
                 round(a$overallReject, 3), ", ",
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
    asf = tolower(x$settings$typeAlphaSpending)
    asfpar = round(x$settings$parameterAlphaSpending, 3)
    asfuser = round(x$settings$userAlphaSpending, 4)

    if (asf == "of") {
      str11 = paste0("Alpha spending: O'Brien-Fleming")
    } else if (asf == "p") {
      str11 = paste0("Alpha spending: Pocock")
    } else if (asf == "wt") {
      str11 = paste0("Alpha spending: Wang-Tsiatis(Delta = ", asfpar, ")")
    } else if (asf == "sfof") {
      str11 = paste0("Alpha spending: Lan-DeMets O'Brien-Fleming")
    } else if (asf == "sfp") {
      str11 = paste0("Alpha spending: Lan-DeMets Pocock")
    } else if (asf == "sfkd") {
      str11 = paste0("Alpha spending: KD(rho = ", asfpar, ")")
    } else if (asf == "sfhsd") {
      str11 = paste0("Alpha spending: HSD(gamma = ", asfpar, ")")
    } else if (asf == "user") {
      str11 = paste0("Alpha spending: User defined(",
                     paste(asfuser, collapse = ","), ")")
    } else {
      str11 = "Alpha spending: None"
    }

    if (!any(is.na(x$settings$spendingTime)) &&
        !all(x$settings$spendingTime == s$informationRates)) {
      str12 = paste0("Spending time: ",
                     paste(round(x$settings$spendingTime, 3), collapse = ","))
      df1 = data.frame(x = rep("", 13))
      colnames(df1) = NULL
      rownames(df1) = c(str1, str2, str3, str4, str5, str6, str7,
                        str8, str9, str10, str11, str12, "")
    } else {
      df1 = data.frame(x = rep("", 12))
      colnames(df1) = NULL
      rownames(df1) = c(str1, str2, str3, str4, str5, str6, str7,
                        str8, str9, str10, str11, "")
    }
  } else {
    df1 = data.frame(x = rep("", 11))
    colnames(df1) = NULL
    rownames(df1) = c(str1, str2, str3, str4, str5, str6, str7,
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
    j3 <- c(1,2,3,11,12)
    j4 <- c(4,5,6,13)

    b[j1] <- lapply(b[j1], formatC, format = "f", digits = 1)
    b[j2] <- lapply(b[j2], formatC, format = "f", digits = 2)
    b[j3] <- lapply(b[j3], formatC, format = "f", digits = 3)
    b[j4] <- lapply(b[j4], formatC, format = "f", digits = 4)

    df = t(b)
    rownames(df) = c("Information rate",
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

    df = t(b)

    rownames(df) = c("Boundary for each 1-sided test (Z)",
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
  a = x$overallResults
  s = x$byStageResults
  k = a$kMax

  if (k>1) {
    str1 = paste0("Group-sequential design with ", k, " stages")
  } else {
    str1 = "Fixed design"
  }

  str1 <- paste(str1, "for one-sample milestone survival probability")

  str2 <- paste0("Milestone: ", round(a$milestone, 3), ", ",
                 "survival probability under H0: ",
                 round(a$survH0, 3), ", ",
                 "under H1: ", round(a$surv, 3))

  str3 <- paste0("Overall power: ",
                 round(a$overallReject, 3), ", ",
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
    asf = tolower(x$settings$typeAlphaSpending)
    asfpar = round(x$settings$parameterAlphaSpending, 3)
    asfuser = round(x$settings$userAlphaSpending, 4)

    bsf = tolower(x$settings$typeBetaSpending)
    bsfpar = round(x$settings$parameterBetaSpending, 3)

    if (asf == "of") {
      str10 = paste0("Alpha spending: O'Brien-Fleming")
    } else if (asf == "p") {
      str10 = paste0("Alpha spending: Pocock")
    } else if (asf == "wt") {
      str10 = paste0("Alpha spending: Wang-Tsiatis(Delta = ", asfpar, ")")
    } else if (asf == "sfof") {
      str10 = paste0("Alpha spending: Lan-DeMets O'Brien-Fleming")
    } else if (asf == "sfp") {
      str10 = paste0("Alpha spending: Lan-DeMets Pocock")
    } else if (asf == "sfkd") {
      str10 = paste0("Alpha spending: KD(rho = ", asfpar, ")")
    } else if (asf == "sfhsd") {
      str10 = paste0("Alpha spending: HSD(gamma = ", asfpar, ")")
    } else if (asf == "user") {
      str10 = paste0("Alpha spending: User defined(",
                    paste(asfuser, collapse = ","), ")")
    } else {
      str10 = "Alpha spending: None"
    }

    if (bsf == "of") {
      str11 = paste0("beta spending: O'Brien-Fleming")
    } else if (bsf == "p") {
      str11 = paste0("beta spending: Pocock")
    } else if (bsf == "wt") {
      str11 = paste0("beta spending: Wang-Tsiatis(Delta = ", bsfpar, ")")
    } else if (bsf == "sfof") {
      str11 = paste0("beta spending: Lan-DeMets O'Brien-Fleming")
    } else if (bsf == "sfp") {
      str11 = paste0("beta spending: Lan-DeMets Pocock")
    } else if (bsf == "sfkd") {
      str11 = paste0("beta spending: KD(rho = ", bsfpar, ")")
    } else if (bsf == "sfhsd") {
      str11 = paste0("beta spending: HSD(gamma = ", bsfpar, ")")
    } else if (bsf == "user") {
      str11 = paste0("beta spending: User defined")
    } else {
      str11 = "beta spending: None"
    }

    if (!any(is.na(x$settings$spendingTime)) &&
        !all(x$settings$spendingTime == s$informationRates)) {
      str12 = paste0("Spending time: ",
                     paste(round(x$settings$spendingTime, 3), collapse = ","))
      df1 = data.frame(x = rep("", 12))
      colnames(df1) = NULL
      rownames(df1) = c(str1, str2, str3, str4, str5, str6, str7,
                        str8, str9, paste(str10, str11, sep = ", "),
                        str12, "")
    } else {
      df1 = data.frame(x = rep("", 11))
      colnames(df1) = NULL
      rownames(df1) = c(str1, str2, str3, str4, str5, str6, str7,
                        str8, str9, paste(str10, str11, sep = ", "), "")
    }
  } else {
    df1 = data.frame(x = rep("", 10))
    colnames(df1) = NULL
    rownames(df1) = c(str1, str2, str3, str4, str5, str6, str7,
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
        (k > 1 && any(x$byStageResults$futilityBounds[1:(k-1)] > -6))) {
      df = t(b)
      rownames(df) = c("Information rate",
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
      df = t(b[,c(1,2,4,6,7,8,9,10,11,12,14,16)])
      rownames(df) = c("Information rate",
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

    df = t(b)

    rownames(df) = c("Efficacy boundary (Z)",
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
  a = x$overallResults
  s = x$byStageResults
  k = a$kMax

  if (k>1) {
    str1 = paste0("Group-sequential design with ", k, " stages")
  } else {
    str1 = "Fixed design"
  }

  str1 <- paste(str1, "for one-sample restricted mean survival time")


  str2 <- paste0("Milestone: ", round(a$milestone, 3), ", ",
                 "restricted mean survival time under H0: ",
                 round(a$rmstH0, 3), ", ",
                 "under H1: ", round(a$rmst, 3))

  str3 <- paste0("Overall power: ",
                 round(a$overallReject, 3), ", ",
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
    asf = tolower(x$settings$typeAlphaSpending)
    asfpar = round(x$settings$parameterAlphaSpending, 3)
    asfuser = round(x$settings$userAlphaSpending, 4)

    bsf = tolower(x$settings$typeBetaSpending)
    bsfpar = round(x$settings$parameterBetaSpending, 3)

    if (asf == "of") {
      str10 = paste0("Alpha spending: O'Brien-Fleming")
    } else if (asf == "p") {
      str10 = paste0("Alpha spending: Pocock")
    } else if (asf == "wt") {
      str10 = paste0("Alpha spending: Wang-Tsiatis(Delta = ", asfpar, ")")
    } else if (asf == "sfof") {
      str10 = paste0("Alpha spending: Lan-DeMets O'Brien-Fleming")
    } else if (asf == "sfp") {
      str10 = paste0("Alpha spending: Lan-DeMets Pocock")
    } else if (asf == "sfkd") {
      str10 = paste0("Alpha spending: KD(rho = ", asfpar, ")")
    } else if (asf == "sfhsd") {
      str10 = paste0("Alpha spending: HSD(gamma = ", asfpar, ")")
    } else if (asf == "user") {
      str10 = paste0("Alpha spending: User defined(",
                     paste(asfuser, collapse = ","), ")")
    } else {
      str10 = "Alpha spending: None"
    }

    if (bsf == "of") {
      str11 = paste0("beta spending: O'Brien-Fleming")
    } else if (bsf == "p") {
      str11 = paste0("beta spending: Pocock")
    } else if (bsf == "wt") {
      str11 = paste0("beta spending: Wang-Tsiatis(Delta = ", bsfpar, ")")
    } else if (bsf == "sfof") {
      str11 = paste0("beta spending: Lan-DeMets O'Brien-Fleming")
    } else if (bsf == "sfp") {
      str11 = paste0("beta spending: Lan-DeMets Pocock")
    } else if (bsf == "sfkd") {
      str11 = paste0("beta spending: KD(rho = ", bsfpar, ")")
    } else if (bsf == "sfhsd") {
      str11 = paste0("beta spending: HSD(gamma = ", bsfpar, ")")
    } else if (bsf == "user") {
      str11 = paste0("beta spending: User defined")
    } else {
      str11 = "beta spending: None"
    }

    if (!any(is.na(x$settings$spendingTime)) &&
        !all(x$settings$spendingTime == s$informationRates)) {
      str12 = paste0("Spending time: ",
                     paste(round(x$settings$spendingTime, 3), collapse = ","))
      df1 = data.frame(x = rep("", 12))
      colnames(df1) = NULL
      rownames(df1) = c(str1, str2, str3, str4, str5, str6, str7,
                        str8, str9, paste(str10, str11, sep = ", "),
                        str12, "")
    } else {
      df1 = data.frame(x = rep("", 11))
      colnames(df1) = NULL
      rownames(df1) = c(str1, str2, str3, str4, str5, str6, str7,
                        str8, str9, paste(str10, str11, sep = ", "), "")
    }
  } else {
    df1 = data.frame(x = rep("", 10))
    colnames(df1) = NULL
    rownames(df1) = c(str1, str2, str3, str4, str5, str6, str7,
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
        (k > 1 && any(x$byStageResults$futilityBounds[1:(k-1)] > -6))) {
      df = t(b)
      rownames(df) = c("Information rate",
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
      df = t(b[,c(1,2,4,6,7,8,9,10,11,12,14,16)])
      rownames(df) = c("Information rate",
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

    df = t(b)

    rownames(df) = c("Efficacy boundary (Z)",
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
  lrchisq = -2*(x$sumstat$loglik0 - x$sumstat$loglik1)
  degrees = x$sumstat$nvar
  pvalue = sapply(1:nrow(x$sumstat), function(i) {
    ifelse(degrees[i] > 0,
           pchisq(lrchisq[i], degrees[i], 0, lower.tail = FALSE),
           NA)
  })
  df1 <- cbind(x$sumstat[, c("n", "nevents", "loglik0", "loglik1")],
               lrchisq = lrchisq, df = degrees, pvalue = pvalue,
               x$sumstat[, c("niter", "dist")])

  p = x$p
  if (p > 0) {
    nreps = nrow(x$parest)/p

    if (!x$robust) {
      if (x$plci) {
        df = data.frame(param = rep(x$param, nreps),
                        coef = x$parest$beta,
                        expcoef = x$parest$expbeta,
                        se = x$parest$sebeta,
                        z = x$parest$z,
                        lower = x$parest$lower,
                        upper = x$parest$upper,
                        p = x$parest$p,
                        method = x$parest$method)

        if (nreps > 1) {
          df = cbind(df, x$parest[, (p+10):ncol(x$parest)])
          colnames(df) <- c("param", "coef", "exp(coef)", "se(coef)", "z",
                            paste("lower", 1-x$alpha),
                            paste("upper", 1-x$alpha), "p", "method",
                            colnames(x$parest)[(p+10):ncol(x$parest)])
        } else {
          colnames(df) <- c("param", "coef", "exp(coef)", "se(coef)", "z",
                            paste("lower", 1-x$alpha),
                            paste("upper", 1-x$alpha), "p", "method")
        }
      } else {
        df = data.frame(param = rep(x$param, nreps),
                        coef = x$parest$beta,
                        expcoef = x$parest$expbeta,
                        se = x$parest$sebeta,
                        z = x$parest$z,
                        p = x$parest$p)

        if (nreps > 1) {
          df = cbind(df, x$parest[, (p+10):ncol(x$parest)])
          colnames(df) <- c("param", "coef", "exp(coef)", "se(coef)", "z",
                            "p", colnames(x$parest)[(p+10):ncol(x$parest)])
        } else {
          colnames(df) <- c("param", "coef", "exp(coef)", "se(coef)", "z",
                            "p")
        }
      }
    } else {
      if (x$plci) {
        df = data.frame(param = rep(x$param, nreps),
                        coef = x$parest$beta,
                        expcoef = x$parest$expbeta,
                        nse = x$parest$sebeta_naive,
                        se = x$parest$sebeta,
                        z = x$parest$z,
                        lower = x$parest$lower,
                        upper = x$parest$upper,
                        p = x$parest$p,
                        method = x$parest$method)

        if (nreps > 1) {
          df = cbind(df, x$parest[, (2*p+11):ncol(x$parest)])
          colnames(df) <- c("param", "coef", "exp(coef)", "se(coef)",
                            "robust se", "z", paste("lower", 1-x$alpha),
                            paste("upper", 1-x$alpha), "p", "method",
                            colnames(x$parest)[(2*p+11):ncol(x$parest)])
        } else {
          colnames(df) <- c("param", "coef", "exp(coef)", "se(coef)",
                            "robust se", "z", paste("lower", 1-x$alpha),
                            paste("upper", 1-x$alpha), "p", "method")
        }
      } else {
        df = data.frame(param = rep(x$param, nreps),
                        coef = x$parest$beta,
                        expcoef = x$parest$expbeta,
                        nse = x$parest$sebeta_naive,
                        se = x$parest$sebeta,
                        z = x$parest$z,
                        p = x$parest$p)

        if (nreps > 1) {
          df = cbind(df, x$parest[, (2*p+11):ncol(x$parest)])
          colnames(df) <- c("param", "coef", "exp(coef)", "se(coef)",
                            "robust se", "z", "p",
                            colnames(x$parest)[(2*p+11):ncol(x$parest)])
        } else {
          colnames(df) <- c("param", "coef", "exp(coef)", "se(coef)",
                            "robust se", "z", "p")
        }
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
  lrchisq = -2*(x$sumstat$loglik0 - x$sumstat$loglik1)
  degrees = x$sumstat$p
  pvalue = sapply(1:nrow(x$sumstat), function(i) {
    ifelse(degrees[i] > 0,
           pchisq(lrchisq[i], degrees[i], 0, lower.tail = FALSE),
           NA)
  })
  df1 <- cbind(x$sumstat[, c("n", "nevents", "loglik0", "loglik1")],
               lrchisq = lrchisq, df = degrees, pvalue = pvalue,
               x$sumstat[, c("scoretest", "niter", "ties")])
  print(df1, ..., na.print = "" , quote = FALSE )
  cat("\n")

  p = x$p
  if (p > 0) {
    nreps = nrow(x$parest)/p

    if (!x$robust) {
      if (x$plci) {
        df = data.frame(param = rep(x$param, nreps),
                        coef = x$parest$beta,
                        expcoef = x$parest$expbeta,
                        se = x$parest$sebeta,
                        z = x$parest$z,
                        lower = x$parest$lower,
                        upper = x$parest$upper,
                        p = x$parest$p,
                        method = x$parest$method)

        if (nreps > 1) {
          df = cbind(df, x$parest[, (p+10):ncol(x$parest)])
          colnames(df) <- c("param", "coef", "exp(coef)", "se(coef)", "z",
                            paste("lower", 1-x$alpha),
                            paste("upper", 1-x$alpha), "p", "method",
                            colnames(x$parest)[(p+10):ncol(x$parest)])
        } else {
          colnames(df) <- c("param", "coef", "exp(coef)", "se(coef)", "z",
                            paste("lower", 1-x$alpha),
                            paste("upper", 1-x$alpha), "p", "method")
        }
      } else {
        df = data.frame(param = rep(x$param, nreps),
                        coef = x$parest$beta,
                        expcoef = x$parest$expbeta,
                        se = x$parest$sebeta,
                        z = x$parest$z,
                        p = x$parest$p)

        if (nreps > 1) {
          df = cbind(df, x$parest[, (p+10):ncol(x$parest)])
          colnames(df) <- c("param", "coef", "exp(coef)", "se(coef)", "z",
                            "p", colnames(x$parest)[(p+10):ncol(x$parest)])
        } else {
          colnames(df) <- c("param", "coef", "exp(coef)", "se(coef)", "z",
                            "p")
        }
      }
    } else {
      if (x$plci) {
        df = data.frame(param = rep(x$param, nreps),
                        coef = x$parest$beta,
                        expcoef = x$parest$expbeta,
                        nse = x$parest$sebeta_naive,
                        se = x$parest$sebeta,
                        z = x$parest$z,
                        lower = x$parest$lower,
                        upper = x$parest$upper,
                        p = x$parest$p,
                        method = x$parest$method)

        if (nreps > 1) {
          df = cbind(df, x$parest[, (2*p+11):ncol(x$parest)])
          colnames(df) <- c("param", "coef", "exp(coef)", "se(coef)",
                            "robust se", "z", paste("lower", 1-x$alpha),
                            paste("upper", 1-x$alpha), "p", "method",
                            colnames(x$parest)[(2*p+11):ncol(x$parest)])
        } else {
          colnames(df) <- c("param", "coef", "exp(coef)", "se(coef)",
                            "robust se", "z", paste("lower", 1-x$alpha),
                            paste("upper", 1-x$alpha), "p", "method")
        }
      } else {
        df = data.frame(param = rep(x$param, nreps),
                        coef = x$parest$beta,
                        expcoef = x$parest$expbeta,
                        nse = x$parest$sebeta_naive,
                        se = x$parest$sebeta,
                        z = x$parest$z,
                        p = x$parest$p)

        if (nreps > 1) {
          df = cbind(df, x$parest[, (2*p+11):ncol(x$parest)])
          colnames(df) <- c("param", "coef", "exp(coef)", "se(coef)",
                            "robust se", "z", "p",
                            colnames(x$parest)[(2*p+11):ncol(x$parest)])
        } else {
          colnames(df) <- c("param", "coef", "exp(coef)", "se(coef)",
                            "robust se", "z", "p")
        }
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
  lrchisq = -2*(x$sumstat$loglik0 - x$sumstat$loglik1)
  degrees = x$sumstat$p - 1
  pvalue = sapply(1:nrow(x$sumstat), function(i) {
    ifelse(degrees[i] > 0,
           pchisq(lrchisq[i], degrees[i], 0, lower.tail = FALSE),
           NA)
  })
  df1 <- cbind(x$sumstat[, c("n", "nevents", "loglik0", "loglik1")],
               lrchisq = lrchisq, df = degrees, pvalue = pvalue,
               x$sumstat[, c("niter", "link", "firth", "flic")])

  p = x$p
  if (p > 0) {
    nreps = nrow(x$parest)/p

    if (!x$robust) {
      if (x$plci) {
        df = data.frame(param = rep(x$param, nreps),
                        coef = x$parest$beta,
                        expcoef = x$parest$expbeta,
                        se = x$parest$sebeta,
                        z = x$parest$z,
                        lower = x$parest$lower,
                        upper = x$parest$upper,
                        p = x$parest$p,
                        method = x$parest$method)

        if (nreps > 1) {
          df = cbind(df, x$parest[, (p+10):ncol(x$parest)])
          colnames(df) <- c("param", "coef", "exp(coef)", "se(coef)", "z",
                            paste("lower", 1-x$alpha),
                            paste("upper", 1-x$alpha), "p", "method",
                            colnames(x$parest)[(p+10):ncol(x$parest)])
        } else {
          colnames(df) <- c("param", "coef", "exp(coef)", "se(coef)", "z",
                            paste("lower", 1-x$alpha),
                            paste("upper", 1-x$alpha), "p", "method")
        }
      } else {
        df = data.frame(param = rep(x$param, nreps),
                        coef = x$parest$beta,
                        expcoef = x$parest$expbeta,
                        se = x$parest$sebeta,
                        z = x$parest$z,
                        p = x$parest$p)

        if (nreps > 1) {
          df = cbind(df, x$parest[, (p+10):ncol(x$parest)])
          colnames(df) <- c("param", "coef", "exp(coef)", "se(coef)", "z",
                            "p", colnames(x$parest)[(p+10):ncol(x$parest)])
        } else {
          colnames(df) <- c("param", "coef", "exp(coef)", "se(coef)", "z",
                            "p")
        }
      }
    } else {
      if (x$plci) {
        df = data.frame(param = rep(x$param, nreps),
                        coef = x$parest$beta,
                        expcoef = x$parest$expbeta,
                        nse = x$parest$sebeta_naive,
                        se = x$parest$sebeta,
                        z = x$parest$z,
                        lower = x$parest$lower,
                        upper = x$parest$upper,
                        p = x$parest$p,
                        method = x$parest$method)

        if (nreps > 1) {
          df = cbind(df, x$parest[, (2*p+11):ncol(x$parest)])
          colnames(df) <- c("param", "coef", "exp(coef)", "se(coef)",
                            "robust se", "z", paste("lower", 1-x$alpha),
                            paste("upper", 1-x$alpha), "p", "method",
                            colnames(x$parest)[(2*p+11):ncol(x$parest)])
        } else {
          colnames(df) <- c("param", "coef", "exp(coef)", "se(coef)",
                            "robust se", "z", paste("lower", 1-x$alpha),
                            paste("upper", 1-x$alpha), "p", "method")
        }
      } else {
        df = data.frame(param = rep(x$param, nreps),
                        coef = x$parest$beta,
                        expcoef = x$parest$expbeta,
                        nse = x$parest$sebeta_naive,
                        se = x$parest$sebeta,
                        z = x$parest$z,
                        p = x$parest$p)

        if (nreps > 1) {
          df = cbind(df, x$parest[, (2*p+11):ncol(x$parest)])
          colnames(df) <- c("param", "coef", "exp(coef)", "se(coef)",
                            "robust se", "z", "p",
                            colnames(x$parest)[(2*p+11):ncol(x$parest)])
        } else {
          colnames(df) <- c("param", "coef", "exp(coef)", "se(coef)",
                            "robust se", "z", "p")
        }
      }
    }
  }

  print(df1, ..., na.print = "" , quote = FALSE )
  cat("\n")
  print(df, ..., na.print = "" , quote = FALSE )
  invisible(x)
}


