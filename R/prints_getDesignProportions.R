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
  a <- x$overallResults
  s <- x$byStageResults
  k <- a$kMax

  if (k>1) {
    str1 <- paste0("Group-sequential design with ", k, " stages")
  } else {
    str1 <- "Fixed design"
  }

  str1 <- paste(str1, "for one-sample proportion")

  str2 <- paste0("Response probability under H0: ",
                 round(a$piH0, 3), ", ",
                 "response probability under H1: ",
                 round(a$pi, 3))

  str3 <- paste0("Overall power: ",
                 round(a$overallReject, 4), ", ",
                 "overall alpha (1-sided): ",
                 round(a$alpha, 4))

  if (x$settings$typeBetaSpending != 'none' ||
      (k > 1 && any(x$byStageResults$futilityBounds[1:(k-1)] > -8)) ||
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
    asf <- tolower(x$settings$typeAlphaSpending)
    asfpar <- round(x$settings$parameterAlphaSpending, 3)
    asfuser <- round(x$settings$userAlphaSpending, 4)

    bsf <- tolower(x$settings$typeBetaSpending)
    bsfpar <- round(x$settings$parameterBetaSpending, 3)
    bsfuser <- round(x$settings$userBetaSpending, 4)

    if (asf == "of") {
      str8 <- paste0("Alpha spending: O'Brien-Fleming")
    } else if (asf == "p") {
      str8 <- paste0("Alpha spending: Pocock")
    } else if (asf == "wt") {
      str8 <- paste0("Alpha spending: Wang-Tsiatis(Delta = ", asfpar, ")")
    } else if (asf == "sfof") {
      str8 <- paste0("Alpha spending: Lan-DeMets O'Brien-Fleming")
    } else if (asf == "sfp") {
      str8 <- paste0("Alpha spending: Lan-DeMets Pocock")
    } else if (asf == "sfkd") {
      str8 <- paste0("Alpha spending: KD(rho = ", asfpar, ")")
    } else if (asf == "sfhsd") {
      str8 <- paste0("Alpha spending: HSD(gamma = ", asfpar, ")")
    } else if (asf == "user") {
      str8 <- paste0("Alpha spending: User defined(",
                     paste(asfuser, collapse = ","), ")")
    } else {
      str8 <- "Alpha spending: None"
    }

    if (bsf == "of") {
      str9 <- paste0("beta spending: O'Brien-Fleming")
    } else if (bsf == "p") {
      str9 <- paste0("beta spending: Pocock")
    } else if (bsf == "wt") {
      str9 <- paste0("beta spending: Wang-Tsiatis(Delta = ", bsfpar, ")")
    } else if (bsf == "sfof") {
      str9 <- paste0("beta spending: Lan-DeMets O'Brien-Fleming")
    } else if (bsf == "sfp") {
      str9 <- paste0("beta spending: Lan-DeMets Pocock")
    } else if (bsf == "sfkd") {
      str9 <- paste0("beta spending: KD(rho = ", bsfpar, ")")
    } else if (bsf == "sfhsd") {
      str9 <- paste0("beta spending: HSD(gamma = ", bsfpar, ")")
    } else if (bsf == "user") {
      str9 <- paste0("beta spending: User defined(",
                     paste(bsfuser, collapse = ","), ")")
    } else {
      str9 <- "beta spending: None"
    }

    if (!any(is.na(x$settings$spendingTime)) &&
        !all(x$settings$spendingTime == s$informationRates)) {
      str10 <- paste0("Spending time: ",
                      paste(round(x$settings$spendingTime, 3), collapse = ","))
      df1 <- data.frame(x = rep("", 10))
      colnames(df1) <- NULL
      rownames(df1) <- c(str1, str2, str3, str4, str5, str6, str7,
                         paste(str8, str9, sep = ", "), str10, "")
    } else {
      df1 <- data.frame(x = rep("", 9))
      colnames(df1) <- NULL
      rownames(df1) <- c(str1, str2, str3, str4, str5, str6, str7,
                         paste(str8, str9, sep = ", "), "")
    }
  } else {
    df1 <- data.frame(x = rep("", 8))
    colnames(df1) <- NULL
    rownames(df1) <- c(str1, str2, str3, str4, str5, str6, str7, "")
  }

  if (k>1) {
    b <- s[, c("informationRates", "efficacyBounds", "futilityBounds",
               "cumulativeRejection", "cumulativeFutility",
               "cumulativeAlphaSpent", "numberOfSubjects",
               "efficacyResponses", "futilityResponses",
               "efficacyP", "futilityP", "information",
               "cumulativeRejectionH0", "cumulativeFutilityH0")]

    # format number of digits after decimal for each column
    j0 <- c(8,9)
    j1 <- 7
    j2 <- 12
    j3 <- c(1,2,3)
    j4 <- c(4,5,6,10,11,13,14)

    b[j0] <- lapply(b[j0], formatC, format = "f", digits = 0)
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
                        "Number of subjects",
                        "Efficacy boundary (# responses)",
                        "Futility boundary (# responses)",
                        "Efficacy boundary (p)",
                        "Futility boundary (p)",
                        "Information",
                        "Cumulative rejection under H0",
                        "Cumulative futility under H0")

    } else {
      df <- t(b[,c(1,2,4,6,7,8,10,12)])
      rownames(df) <- c("Information rate",
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
    j0 <- 2
    j3 <- 1
    j4 <- 3

    b[j0] <- lapply(b[j0], formatC, format = "f", digits = 0)
    b[j3] <- lapply(b[j3], formatC, format = "f", digits = 3)
    b[j4] <- lapply(b[j4], formatC, format = "f", digits = 4)

    df <- t(b)

    if (!x$settings$normalApproximation) {
      rownames(df) <- c("Efficacy boundary (Z)",
                        "Efficacy boundary (# responses)",
                        "Efficacy boundary (p)")
    } else {
      rownames(df) <- c("Efficacy boundary (Z)",
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
  a <- x$overallResults
  s <- x$byStageResults
  k <- a$kMax

  if (k>1) {
    str1 <- paste0("Group-sequential design with ", k, " stages")
  } else {
    str1 <- "Fixed design"
  }

  str1 <- paste(str1, "for McNemar's test")

  str2 <- paste0("Proportion of discordant pairs: ",
                 round(a$pDiscordant, 3), ", ",
                 "risk difference: ", round(a$riskDiff, 3))

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
    asf <- tolower(x$settings$typeAlphaSpending)
    asfpar <- round(x$settings$parameterAlphaSpending, 3)
    asfuser <- round(x$settings$userAlphaSpending, 4)

    bsf <- tolower(x$settings$typeBetaSpending)
    bsfpar <- round(x$settings$parameterBetaSpending, 3)
    bsfuser <- round(x$settings$userBetaSpending, 4)

    if (asf == "of") {
      str8 <- paste0("Alpha spending: O'Brien-Fleming")
    } else if (asf == "p") {
      str8 <- paste0("Alpha spending: Pocock")
    } else if (asf == "wt") {
      str8 <- paste0("Alpha spending: Wang-Tsiatis(Delta = ", asfpar, ")")
    } else if (asf == "sfof") {
      str8 <- paste0("Alpha spending: Lan-DeMets O'Brien-Fleming")
    } else if (asf == "sfp") {
      str8 <- paste0("Alpha spending: Lan-DeMets Pocock")
    } else if (asf == "sfkd") {
      str8 <- paste0("Alpha spending: KD(rho = ", asfpar, ")")
    } else if (asf == "sfhsd") {
      str8 <- paste0("Alpha spending: HSD(gamma = ", asfpar, ")")
    } else if (asf == "user") {
      str8 <- paste0("Alpha spending: User defined(",
                     paste(asfuser, collapse = ","), ")")
    } else {
      str8 <- "Alpha spending: None"
    }

    if (bsf == "of") {
      str9 <- paste0("beta spending: O'Brien-Fleming")
    } else if (bsf == "p") {
      str9 <- paste0("beta spending: Pocock")
    } else if (bsf == "wt") {
      str9 <- paste0("beta spending: Wang-Tsiatis(Delta = ", bsfpar, ")")
    } else if (bsf == "sfof") {
      str9 <- paste0("beta spending: Lan-DeMets O'Brien-Fleming")
    } else if (bsf == "sfp") {
      str9 <- paste0("beta spending: Lan-DeMets Pocock")
    } else if (bsf == "sfkd") {
      str9 <- paste0("beta spending: KD(rho = ", bsfpar, ")")
    } else if (bsf == "sfhsd") {
      str9 <- paste0("beta spending: HSD(gamma = ", bsfpar, ")")
    } else if (bsf == "user") {
      str9 <- paste0("beta spending: User defined(",
                     paste(bsfuser, collapse = ","), ")")
    } else {
      str9 <- "beta spending: None"
    }

    if (!any(is.na(x$settings$spendingTime)) &&
        !all(x$settings$spendingTime == s$informationRates)) {
      str10 <- paste0("Spending time: ",
                      paste(round(x$settings$spendingTime, 3), collapse = ","))
      df1 <- data.frame(x = rep("", 10))
      colnames(df1) <- NULL
      rownames(df1) <- c(str1, str2, str3, str4, str5, str6, str7,
                         paste(str8, str9, sep = ", "), str10, "")
    } else {
      df1 <- data.frame(x = rep("", 9))
      colnames(df1) <- NULL
      rownames(df1) <- c(str1, str2, str3, str4, str5, str6, str7,
                         paste(str8, str8, sep = ", "), "")
    }
  } else {
    df1 <- data.frame(x = rep("", 8))
    colnames(df1) <- NULL
    rownames(df1) <- c(str1, str2, str3, str4, str5, str6, str7, "")
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
    j3 <- c(1,2,3,8,9)
    j4 <- c(4,5,6,10,11,13,14)

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
                        "Number of subjects",
                        "Efficacy boundary (risk diff)",
                        "Futility boundary (risk diff)",
                        "Efficacy boundary (p)",
                        "Futility boundary (p)",
                        "Information",
                        "Cumulative rejection under H0",
                        "Cumulative futility under H0")

    } else {
      df <- t(b[,c(1,2,4,6,7,8,10,12)])
      rownames(df) <- c("Information rate",
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

    df <- t(b)

    rownames(df) <- c("Efficacy boundary (Z)",
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
  a <- x$overallResults
  s <- x$byStageResults
  k <- a$kMax

  if (k>1) {
    str1 <- paste0("Group-sequential design with ", k, " stages")
  } else {
    str1 <- "Fixed design"
  }

  str1 <- paste(str1,
                "for two-sample risk difference")

  str2 <- paste0("Risk difference under H0: ", round(a$riskDiffH0, 3), ", ",
                 "proportion on treatment: ", round(a$pi1, 3), ", ",
                 "proportion on control: ", round(a$pi2, 3))

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
    asf <- tolower(x$settings$typeAlphaSpending)
    asfpar <- round(x$settings$parameterAlphaSpending, 3)
    asfuser <- round(x$settings$userAlphaSpending, 4)

    bsf <- tolower(x$settings$typeBetaSpending)
    bsfpar <- round(x$settings$parameterBetaSpending, 3)
    bsfuser <- round(x$settings$userBetaSpending, 4)

    if (asf == "of") {
      str8 <- paste0("Alpha spending: O'Brien-Fleming")
    } else if (asf == "p") {
      str8 <- paste0("Alpha spending: Pocock")
    } else if (asf == "wt") {
      str8 <- paste0("Alpha spending: Wang-Tsiatis(Delta = ", asfpar, ")")
    } else if (asf == "sfof") {
      str8 <- paste0("Alpha spending: Lan-DeMets O'Brien-Fleming")
    } else if (asf == "sfp") {
      str8 <- paste0("Alpha spending: Lan-DeMets Pocock")
    } else if (asf == "sfkd") {
      str8 <- paste0("Alpha spending: KD(rho = ", asfpar, ")")
    } else if (asf == "sfhsd") {
      str8 <- paste0("Alpha spending: HSD(gamma = ", asfpar, ")")
    } else if (asf == "user") {
      str8 <- paste0("Alpha spending: User defined(",
                     paste(asfuser, collapse = ","), ")")
    } else {
      str8 <- "Alpha spending: None"
    }

    if (bsf == "of") {
      str9 <- paste0("beta spending: O'Brien-Fleming")
    } else if (bsf == "p") {
      str9 <- paste0("beta spending: Pocock")
    } else if (bsf == "wt") {
      str9 <- paste0("beta spending: Wang-Tsiatis(Delta = ", bsfpar, ")")
    } else if (bsf == "sfof") {
      str9 <- paste0("beta spending: Lan-DeMets O'Brien-Fleming")
    } else if (bsf == "sfp") {
      str9 <- paste0("beta spending: Lan-DeMets Pocock")
    } else if (bsf == "sfkd") {
      str9 <- paste0("beta spending: KD(rho = ", bsfpar, ")")
    } else if (bsf == "sfhsd") {
      str9 <- paste0("beta spending: HSD(gamma = ", bsfpar, ")")
    } else if (bsf == "user") {
      str9 <- paste0("beta spending: User defined(",
                     paste(bsfuser, collapse = ","), ")")
    } else {
      str9 <- "beta spending: None"
    }

    if (!any(is.na(x$settings$spendingTime)) &&
        !all(x$settings$spendingTime == s$informationRates)) {
      str10 <- paste0("Spending time: ",
                      paste(round(x$settings$spendingTime, 3), collapse = ","))
      df1 <- data.frame(x = rep("", 10))
      colnames(df1) <- NULL
      rownames(df1) <- c(str1, str2, str3, str4, str5, str6, str7,
                         paste(str8, str9, sep = ", "), str10, "")
    } else {
      df1 <- data.frame(x = rep("", 9))
      colnames(df1) <- NULL
      rownames(df1) <- c(str1, str2, str3, str4, str5, str6, str7,
                         paste(str8, str9, sep = ", "), "")
    }
  } else {
    df1 <- data.frame(x = rep("", 8))
    colnames(df1) <- NULL
    rownames(df1) <- c(str1, str2, str3, str4, str5, str6, str7, "")
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
    j3 <- c(1,2,3,8,9)
    j4 <- c(4,5,6,10,11,13,14)

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
                        "Number of subjects",
                        "Efficacy boundary (risk diff)",
                        "Futility boundary (risk diff)",
                        "Efficacy boundary (p)",
                        "Futility boundary (p)",
                        "Information",
                        "Cumulative rejection under H0",
                        "Cumulative futility under H0")

    } else {
      df <- t(b[,c(1,2,4,6,7,8,10,12)])
      rownames(df) <- c("Information rate",
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

    df <- t(b)

    rownames(df) <- c("Efficacy boundary (Z)",
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
  a <- x$overallResults
  s <- x$byStageResults
  k <- a$kMax

  if (k>1) {
    str1 <- paste0("Group-sequential design with ", k, " stages")
  } else {
    str1 <- "Fixed design"
  }

  str1 <- paste(str1, "for two-sample risk ratio")

  str2 <- paste0("Risk ratio under H0: ", round(a$riskRatioH0, 3), ", ",
                 "proportion on treatment: ", round(a$pi1, 3), ", ",
                 "proportion on control: ", round(a$pi2, 3))

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
    asf <- tolower(x$settings$typeAlphaSpending)
    asfpar <- round(x$settings$parameterAlphaSpending, 3)
    asfuser <- round(x$settings$userAlphaSpending, 4)

    bsf <- tolower(x$settings$typeBetaSpending)
    bsfpar <- round(x$settings$parameterBetaSpending, 3)
    bsfuser <- round(x$settings$userBetaSpending, 4)

    if (asf == "of") {
      str8 <- paste0("Alpha spending: O'Brien-Fleming")
    } else if (asf == "p") {
      str8 <- paste0("Alpha spending: Pocock")
    } else if (asf == "wt") {
      str8 <- paste0("Alpha spending: Wang-Tsiatis(Delta = ", asfpar, ")")
    } else if (asf == "sfof") {
      str8 <- paste0("Alpha spending: Lan-DeMets O'Brien-Fleming")
    } else if (asf == "sfp") {
      str8 <- paste0("Alpha spending: Lan-DeMets Pocock")
    } else if (asf == "sfkd") {
      str8 <- paste0("Alpha spending: KD(rho = ", asfpar, ")")
    } else if (asf == "sfhsd") {
      str8 <- paste0("Alpha spending: HSD(gamma = ", asfpar, ")")
    } else if (asf == "user") {
      str8 <- paste0("Alpha spending: User defined(",
                     paste(asfuser, collapse = ","), ")")
    } else {
      str8 <- "Alpha spending: None"
    }

    if (bsf == "of") {
      str9 <- paste0("beta spending: O'Brien-Fleming")
    } else if (bsf == "p") {
      str9 <- paste0("beta spending: Pocock")
    } else if (bsf == "wt") {
      str9 <- paste0("beta spending: Wang-Tsiatis(Delta = ", bsfpar, ")")
    } else if (bsf == "sfof") {
      str9 <- paste0("beta spending: Lan-DeMets O'Brien-Fleming")
    } else if (bsf == "sfp") {
      str9 <- paste0("beta spending: Lan-DeMets Pocock")
    } else if (bsf == "sfkd") {
      str9 <- paste0("beta spending: KD(rho = ", bsfpar, ")")
    } else if (bsf == "sfhsd") {
      str9 <- paste0("beta spending: HSD(gamma = ", bsfpar, ")")
    } else if (bsf == "user") {
      str9 <- paste0("beta spending: User defined(",
                     paste(bsfuser, collapse = ","), ")")
    } else {
      str9 <- "beta spending: None"
    }

    if (!any(is.na(x$settings$spendingTime)) &&
        !all(x$settings$spendingTime == s$informationRates)) {
      str10 <- paste0("Spending time: ",
                      paste(round(x$settings$spendingTime, 3), collapse = ","))
      df1 <- data.frame(x = rep("", 10))
      colnames(df1) <- NULL
      rownames(df1) <- c(str1, str2, str3, str4, str5, str6, str7,
                         paste(str8, str9, sep = ", "), str10, "")
    } else {
      df1 <- data.frame(x = rep("", 9))
      colnames(df1) <- NULL
      rownames(df1) <- c(str1, str2, str3, str4, str5, str6, str7,
                         paste(str8, str9, sep = ", "), "")
    }
  } else {
    df1 <- data.frame(x = rep("", 8))
    colnames(df1) <- NULL
    rownames(df1) <- c(str1, str2, str3, str4, str5, str6, str7, "")
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
    j3 <- c(1,2,3,8,9)
    j4 <- c(4,5,6,10,11,13,14)

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
                        "Number of subjects",
                        "Efficacy boundary (risk ratio)",
                        "Futility boundary (risk ratio)",
                        "Efficacy boundary (p)",
                        "Futility boundary (p)",
                        "Information",
                        "Cumulative rejection under H0",
                        "Cumulative futility under H0")

    } else {
      df <- t(b[,c(1,2,4,6,7,8,10,12)])
      rownames(df) <- c("Information rate",
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

    df <- t(b)

    rownames(df) <- c("Efficacy boundary (Z)",
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
  a <- x$overallResults
  s <- x$byStageResults
  k <- a$kMax

  if (k>1) {
    str1 <- paste0("Group-sequential design with ", k, " stages")
  } else {
    str1 <- "Fixed design"
  }

  str1 <- paste(str1, "for two-sample risk ratio based on the score test")

  str2 <- paste0("Risk ratio under H0: ", round(a$riskRatioH0, 3), ", ",
                 "proportion on treatment: ", round(a$pi1, 3), ", ",
                 "proportion on control: ", round(a$pi2, 3))

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
    asf <- tolower(x$settings$typeAlphaSpending)
    asfpar <- round(x$settings$parameterAlphaSpending, 3)
    asfuser <- round(x$settings$userAlphaSpending, 4)

    bsf <- tolower(x$settings$typeBetaSpending)
    bsfpar <- round(x$settings$parameterBetaSpending, 3)
    bsfuser <- round(x$settings$userBetaSpending, 4)

    if (asf == "of") {
      str8 <- paste0("Alpha spending: O'Brien-Fleming")
    } else if (asf == "p") {
      str8 <- paste0("Alpha spending: Pocock")
    } else if (asf == "wt") {
      str8 <- paste0("Alpha spending: Wang-Tsiatis(Delta = ", asfpar, ")")
    } else if (asf == "sfof") {
      str8 <- paste0("Alpha spending: Lan-DeMets O'Brien-Fleming")
    } else if (asf == "sfp") {
      str8 <- paste0("Alpha spending: Lan-DeMets Pocock")
    } else if (asf == "sfkd") {
      str8 <- paste0("Alpha spending: KD(rho = ", asfpar, ")")
    } else if (asf == "sfhsd") {
      str8 <- paste0("Alpha spending: HSD(gamma = ", asfpar, ")")
    } else if (asf == "user") {
      str8 <- paste0("Alpha spending: User defined(",
                     paste(asfuser, collapse = ","), ")")
    } else {
      str8 <- "Alpha spending: None"
    }

    if (bsf == "of") {
      str9 <- paste0("beta spending: O'Brien-Fleming")
    } else if (bsf == "p") {
      str9 <- paste0("beta spending: Pocock")
    } else if (bsf == "wt") {
      str9 <- paste0("beta spending: Wang-Tsiatis(Delta = ", bsfpar, ")")
    } else if (bsf == "sfof") {
      str9 <- paste0("beta spending: Lan-DeMets O'Brien-Fleming")
    } else if (bsf == "sfp") {
      str9 <- paste0("beta spending: Lan-DeMets Pocock")
    } else if (bsf == "sfkd") {
      str9 <- paste0("beta spending: KD(rho = ", bsfpar, ")")
    } else if (bsf == "sfhsd") {
      str9 <- paste0("beta spending: HSD(gamma = ", bsfpar, ")")
    } else if (bsf == "user") {
      str9 <- paste0("beta spending: User defined(",
                     paste(bsfuser, collapse = ","), ")")
    } else {
      str9 <- "beta spending: None"
    }

    if (!any(is.na(x$settings$spendingTime)) &&
        !all(x$settings$spendingTime == s$informationRates)) {
      str10 <- paste0("Spending time: ",
                      paste(round(x$settings$spendingTime, 3), collapse = ","))
      df1 <- data.frame(x = rep("", 10))
      colnames(df1) <- NULL
      rownames(df1) <- c(str1, str2, str3, str4, str5, str6, str7,
                         paste(str8, str9, sep = ", "), str10, "")
    } else {
      df1 <- data.frame(x = rep("", 9))
      colnames(df1) <- NULL
      rownames(df1) <- c(str1, str2, str3, str4, str5, str6, str7,
                         paste(str8, str9, sep = ", "), "")
    }
  } else {
    df1 <- data.frame(x = rep("", 8))
    colnames(df1) <- NULL
    rownames(df1) <- c(str1, str2, str3, str4, str5, str6, str7, "")
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
    j3 <- c(1,2,3,8,9)
    j4 <- c(4,5,6,10,11,13,14)

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
                        "Number of subjects",
                        "Efficacy boundary (risk ratio score)",
                        "Futility boundary (risk ratio score)",
                        "Efficacy boundary (p)",
                        "Futility boundary (p)",
                        "Information",
                        "Cumulative rejection under H0",
                        "Cumulative futility under H0")

    } else {
      df <- t(b[,c(1,2,4,6,7,8,10,12)])
      rownames(df) <- c("Information rate",
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

    df <- t(b)

    rownames(df) <- c("Efficacy boundary (Z)",
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
  a <- x$overallResults
  s <- x$byStageResults
  k <- a$kMax

  if (k>1) {
    str1 <- paste0("Group-sequential design with ", k, " stages")
  } else {
    str1 <- "Fixed design"
  }

  str1 <- paste(str1, "for two-sample odds ratio")

  str2 <- paste0("Odds ratio under H0: ", round(a$oddsRatioH0, 3), ", ",
                 "proportion on treatment: ", round(a$pi1, 3), ", ",
                 "proportion on control: ", round(a$pi2, 3))

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
    asf <- tolower(x$settings$typeAlphaSpending)
    asfpar <- round(x$settings$parameterAlphaSpending, 3)
    asfuser <- round(x$settings$userAlphaSpending, 4)

    bsf <- tolower(x$settings$typeBetaSpending)
    bsfpar <- round(x$settings$parameterBetaSpending, 3)
    bsfuser <- round(x$settings$userBetaSpending, 4)

    if (asf == "of") {
      str8 <- paste0("Alpha spending: O'Brien-Fleming")
    } else if (asf == "p") {
      str8 <- paste0("Alpha spending: Pocock")
    } else if (asf == "wt") {
      str8 <- paste0("Alpha spending: Wang-Tsiatis(Delta = ", asfpar, ")")
    } else if (asf == "sfof") {
      str8 <- paste0("Alpha spending: Lan-DeMets O'Brien-Fleming")
    } else if (asf == "sfp") {
      str8 <- paste0("Alpha spending: Lan-DeMets Pocock")
    } else if (asf == "sfkd") {
      str8 <- paste0("Alpha spending: KD(rho = ", asfpar, ")")
    } else if (asf == "sfhsd") {
      str8 <- paste0("Alpha spending: HSD(gamma = ", asfpar, ")")
    } else if (asf == "user") {
      str8 <- paste0("Alpha spending: User defined(",
                     paste(asfuser, collapse = ","), ")")
    } else {
      str8 <- "Alpha spending: None"
    }

    if (bsf == "of") {
      str9 <- paste0("beta spending: O'Brien-Fleming")
    } else if (bsf == "p") {
      str9 <- paste0("beta spending: Pocock")
    } else if (bsf == "wt") {
      str9 <- paste0("beta spending: Wang-Tsiatis(Delta = ", bsfpar, ")")
    } else if (bsf == "sfof") {
      str9 <- paste0("beta spending: Lan-DeMets O'Brien-Fleming")
    } else if (bsf == "sfp") {
      str9 <- paste0("beta spending: Lan-DeMets Pocock")
    } else if (bsf == "sfkd") {
      str9 <- paste0("beta spending: KD(rho = ", bsfpar, ")")
    } else if (bsf == "sfhsd") {
      str9 <- paste0("beta spending: HSD(gamma = ", bsfpar, ")")
    } else if (bsf == "user") {
      str9 <- paste0("beta spending: User defined(",
                     paste(bsfuser, collapse = ","), ")")
    } else {
      str9 <- "beta spending: None"
    }

    if (!any(is.na(x$settings$spendingTime)) &&
        !all(x$settings$spendingTime == s$informationRates)) {
      str10 <- paste0("Spending time: ",
                      paste(round(x$settings$spendingTime, 3), collapse = ","))
      df1 <- data.frame(x = rep("", 10))
      colnames(df1) <- NULL
      rownames(df1) <- c(str1, str2, str3, str4, str5, str6, str7,
                         paste(str8, str9, sep = ", "), str10, "")
    } else {
      df1 <- data.frame(x = rep("", 9))
      colnames(df1) <- NULL
      rownames(df1) <- c(str1, str2, str3, str4, str5, str6, str7,
                         paste(str8, str9, sep = ", "), "")
    }
  } else {
    df1 <- data.frame(x = rep("", 8))
    colnames(df1) <- NULL
    rownames(df1) <- c(str1, str2, str3, str4, str5, str6, str7, "")
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
    j3 <- c(1,2,3,8,9)
    j4 <- c(4,5,6,10,11,13,14)

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
                        "Number of subjects",
                        "Efficacy boundary (odds ratio)",
                        "Futility boundary (odds ratio)",
                        "Efficacy boundary (p)",
                        "Futility boundary (p)",
                        "Information",
                        "Cumulative rejection under H0",
                        "Cumulative futility under H0")

    } else {
      df <- t(b[,c(1,2,4,6,7,8,10,12)])
      rownames(df) <- c("Information rate",
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

    df <- t(b)

    rownames(df) <- c("Efficacy boundary (Z)",
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
  a <- x$overallResults
  s <- x$byStageResults
  k <- a$kMax

  if (k>1) {
    str1 <- paste0("Group-sequential design with ", k, " stages")
  } else {
    str1 <- "Fixed design"
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
                 round(a$overallReject, 4), ", ",
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
    asf <- tolower(x$settings$typeAlphaSpending)
    asfpar <- round(x$settings$parameterAlphaSpending, 3)
    asfuser <- round(x$settings$userAlphaSpending, 4)

    if (asf == "of") {
      str8 <- paste0("Alpha spending: O'Brien-Fleming")
    } else if (asf == "p") {
      str8 <- paste0("Alpha spending: Pocock")
    } else if (asf == "wt") {
      str8 <- paste0("Alpha spending: Wang-Tsiatis(Delta = ", asfpar, ")")
    } else if (asf == "sfof") {
      str8 <- paste0("Alpha spending: Lan-DeMets O'Brien-Fleming")
    } else if (asf == "sfp") {
      str8 <- paste0("Alpha spending: Lan-DeMets Pocock")
    } else if (asf == "sfkd") {
      str8 <- paste0("Alpha spending: KD(rho = ", asfpar, ")")
    } else if (asf == "sfhsd") {
      str8 <- paste0("Alpha spending: HSD(gamma = ", asfpar, ")")
    } else if (asf == "user") {
      str8 <- paste0("Alpha spending: User defined(",
                     paste(asfuser, collapse = ","), ")")
    } else {
      str8 <- "Alpha spending: None"
    }

    if (!any(is.na(x$settings$spendingTime)) &&
        !all(x$settings$spendingTime == s$informationRates)) {
      str9 <- paste0("Spending time: ",
                     paste(round(x$settings$spendingTime, 3), collapse = ","))
      df1 <- data.frame(x = rep("", 10))
      colnames(df1) <- NULL
      rownames(df1) <- c(str1, str2, str3, str4, str5, str6, str7,
                         str8, str9, "")
    } else {
      df1 <- data.frame(x = rep("", 9))
      colnames(df1) <- NULL
      rownames(df1) <- c(str1, str2, str3, str4, str5, str6, str7, str8, "")
    }
  } else {
    df1 <- data.frame(x = rep("", 8))
    colnames(df1) <- NULL
    rownames(df1) <- c(str1, str2, str3, str4, str5, str6, str7, "")
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
    j3 <- c(1,2,8,9)
    j4 <- c(3,4,5,6,10)

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

    df <- t(b)

    rownames(df) <- c("Boundary for each 1-sided test (Z)",
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
  a <- x$overallResults
  s <- x$byStageResults
  k <- a$kMax

  if (k>1) {
    str1 <- paste0("Group-sequential design with ", k, " stages")
  } else {
    str1 <- "Fixed design"
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
                 round(a$overallReject, 4), ", ",
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
    asf <- tolower(x$settings$typeAlphaSpending)
    asfpar <- round(x$settings$parameterAlphaSpending, 3)
    asfuser <- round(x$settings$userAlphaSpending, 4)

    if (asf == "of") {
      str8 <- paste0("Alpha spending: O'Brien-Fleming")
    } else if (asf == "p") {
      str8 <- paste0("Alpha spending: Pocock")
    } else if (asf == "wt") {
      str8 <- paste0("Alpha spending: Wang-Tsiatis(Delta = ", asfpar, ")")
    } else if (asf == "sfof") {
      str8 <- paste0("Alpha spending: Lan-DeMets O'Brien-Fleming")
    } else if (asf == "sfp") {
      str8 <- paste0("Alpha spending: Lan-DeMets Pocock")
    } else if (asf == "sfkd") {
      str8 <- paste0("Alpha spending: KD(rho = ", asfpar, ")")
    } else if (asf == "sfhsd") {
      str8 <- paste0("Alpha spending: HSD(gamma = ", asfpar, ")")
    } else if (asf == "user") {
      str8 <- paste0("Alpha spending: User defined(",
                     paste(asfuser, collapse = ","), ")")
    } else {
      str8 <- "Alpha spending: None"
    }

    if (!any(is.na(x$settings$spendingTime)) &&
        !all(x$settings$spendingTime == s$informationRates)) {
      str9 <- paste0("Spending time: ",
                     paste(round(x$settings$spendingTime, 3), collapse = ","))
      df1 <- data.frame(x = rep("", 10))
      colnames(df1) <- NULL
      rownames(df1) <- c(str1, str2, str3, str4, str5, str6, str7,
                         str8, str9, "")
    } else {
      df1 <- data.frame(x = rep("", 9))
      colnames(df1) <- NULL
      rownames(df1) <- c(str1, str2, str3, str4, str5, str6, str7, str8, "")
    }
  } else {
    df1 <- data.frame(x = rep("", 8))
    colnames(df1) <- NULL
    rownames(df1) <- c(str1, str2, str3, str4, str5, str6, str7, "")
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
    j3 <- c(1,2,8,9)
    j4 <- c(3,4,5,6,10)

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

    df <- t(b)

    rownames(df) <- c("Boundary for each 1-sided test (Z)",
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
  a <- x$overallResults
  s <- x$byStageResults
  k <- a$kMax

  if (k>1) {
    str1 <- paste0("Group-sequential design with ", k, " stages")
  } else {
    str1 <- "Fixed design"
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
                 round(a$overallReject, 4), ", ",
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
    asf <- tolower(x$settings$typeAlphaSpending)
    asfpar <- round(x$settings$parameterAlphaSpending, 3)
    asfuser <- round(x$settings$userAlphaSpending, 4)

    if (asf == "of") {
      str8 <- paste0("Alpha spending: O'Brien-Fleming")
    } else if (asf == "p") {
      str8 <- paste0("Alpha spending: Pocock")
    } else if (asf == "wt") {
      str8 <- paste0("Alpha spending: Wang-Tsiatis(Delta = ", asfpar, ")")
    } else if (asf == "sfof") {
      str8 <- paste0("Alpha spending: Lan-DeMets O'Brien-Fleming")
    } else if (asf == "sfp") {
      str8 <- paste0("Alpha spending: Lan-DeMets Pocock")
    } else if (asf == "sfkd") {
      str8 <- paste0("Alpha spending: KD(rho = ", asfpar, ")")
    } else if (asf == "sfhsd") {
      str8 <- paste0("Alpha spending: HSD(gamma = ", asfpar, ")")
    } else if (asf == "user") {
      str8 <- paste0("Alpha spending: User defined(",
                     paste(asfuser, collapse = ","), ")")
    } else {
      str8 <- "Alpha spending: None"
    }

    if (!any(is.na(x$settings$spendingTime)) &&
        !all(x$settings$spendingTime == s$informationRates)) {
      str9 <- paste0("Spending time: ",
                     paste(round(x$settings$spendingTime, 3), collapse = ","))
      df1 <- data.frame(x = rep("", 10))
      colnames(df1) <- NULL
      rownames(df1) <- c(str1, str2, str3, str4, str5, str6, str7,
                         str8, str9, "")
    } else {
      df1 <- data.frame(x = rep("", 9))
      colnames(df1) <- NULL
      rownames(df1) <- c(str1, str2, str3, str4, str5, str6, str7, str8, "")
    }
  } else {
    df1 <- data.frame(x = rep("", 8))
    colnames(df1) <- NULL
    rownames(df1) <- c(str1, str2, str3, str4, str5, str6, str7, "")
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
    j3 <- c(1,2,8,9)
    j4 <- c(3,4,5,6,10)

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

    df <- t(b)

    rownames(df) <- c("Boundary for each 1-sided test (Z)",
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
  df1 <- data.frame(x = rep("", 2))
  colnames(df1) <- NULL
  rownames(df1) <- c(str1, "")

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

  df3 <- data.frame(x = rep("", 7))
  colnames(df3) <- NULL
  rownames(df3) <- c(str3, str4, str5, str6, str7, str8, str9)

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
  df1 <- data.frame(x = rep("", 2))
  colnames(df1) <- NULL
  rownames(df1) <- c(str1, "")

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

  df3 <- data.frame(x = rep("", 7))
  colnames(df3) <- NULL
  rownames(df3) <- c(str3, str4, str5, str6, str7, str8, str9)

  print(df1, ..., na.print = "" , quote = FALSE )
  print(df2, ..., na.print = "" , quote = FALSE )
  print(df3, ..., na.print = "" , quote = FALSE )
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
  rownames(df1) <- NULL
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
  rownames(df1) <- NULL
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
  rownames(df1) <- NULL
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
  rownames(df1) <- NULL
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
  rownames(df1) <- NULL
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
  rownames(df1) <- NULL
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
  rownames(df1) <- NULL
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
  rownames(df1) <- NULL
  print(df1, ..., na.print = "" , quote = FALSE )
  invisible(x)
}
