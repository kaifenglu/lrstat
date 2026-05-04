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
  a <- x$overallResults
  s <- x$byStageResults
  k <- a$kMax

  if (k>1) {
    str1 <- paste0("Group-sequential design with ", k, " stages")
  } else {
    str1 <- "Fixed design"
  }

  str1 <- paste(str1, "for one-sample mean")

  str2 <- paste0("Mean under H0: ",
                 round(a$meanH0, 3), ", ",
                 "mean under H1: ",
                 round(a$mean, 3), ", ",
                 "standard deviation: ", round(a$stDev, 3))

  str3 <- paste0("Overall power: ",
                 round(a$overallReject, 4), ", ",
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
    asf <- tolower(x$settings$typeAlphaSpending)
    asfpar <- round(x$settings$parameterAlphaSpending, 3)
    asfuser <- round(x$settings$userAlphaSpending, 4)

    bsf <- tolower(x$settings$typeBetaSpending)
    bsfpar <- round(x$settings$parameterBetaSpending, 3)
    bsfuser <- round(x$settings$userBetaSpending, 4)

    if (asf == "of") {
      str7 <- paste0("Alpha spending: O'Brien-Fleming")
    } else if (asf == "p") {
      str7 <- paste0("Alpha spending: Pocock")
    } else if (asf == "wt") {
      str7 <- paste0("Alpha spending: Wang-Tsiatis(Delta = ", asfpar, ")")
    } else if (asf == "sfof") {
      str7 <- paste0("Alpha spending: Lan-DeMets O'Brien-Fleming")
    } else if (asf == "sfp") {
      str7 <- paste0("Alpha spending: Lan-DeMets Pocock")
    } else if (asf == "sfkd") {
      str7 <- paste0("Alpha spending: KD(rho = ", asfpar, ")")
    } else if (asf == "sfhsd") {
      str7 <- paste0("Alpha spending: HSD(gamma = ", asfpar, ")")
    } else if (asf == "user") {
      str7 <- paste0("Alpha spending: User defined(",
                     paste(asfuser, collapse = ","), ")")
    } else {
      str7 <- "Alpha spending: None"
    }

    if (bsf == "of") {
      str8 <- paste0("beta spending: O'Brien-Fleming")
    } else if (bsf == "p") {
      str8 <- paste0("beta spending: Pocock")
    } else if (bsf == "wt") {
      str8 <- paste0("beta spending: Wang-Tsiatis(Delta = ", bsfpar, ")")
    } else if (bsf == "sfof") {
      str8 <- paste0("beta spending: Lan-DeMets O'Brien-Fleming")
    } else if (bsf == "sfp") {
      str8 <- paste0("beta spending: Lan-DeMets Pocock")
    } else if (bsf == "sfkd") {
      str8 <- paste0("beta spending: KD(rho = ", bsfpar, ")")
    } else if (bsf == "sfhsd") {
      str8 <- paste0("beta spending: HSD(gamma = ", bsfpar, ")")
    } else if (bsf == "user") {
      str8 <- paste0("beta spending: User defined(",
                     paste(bsfuser, collapse = ","), ")")
    } else {
      str8 <- "beta spending: None"
    }

    if (!any(is.na(x$settings$spendingTime)) &&
        !all(x$settings$spendingTime == s$informationRates)) {
      str9 <- paste0("Spending time: ",
                     paste(round(x$settings$spendingTime, 3), collapse = ","))
      df1 <- data.frame(x = rep("", 9))
      colnames(df1) <- NULL
      rownames(df1) <- c(str1, str2, str3, str4, str5, str6,
                         paste(str7, str8, sep = ", "), str9, "")
    } else {
      df1 <- data.frame(x = rep("", 8))
      colnames(df1) <- NULL
      rownames(df1) <- c(str1, str2, str3, str4, str5, str6,
                         paste(str7, str8, sep = ", "), "")
    }
  } else {
    df1 <- data.frame(x = rep("", 7))
    colnames(df1) <- NULL
    rownames(df1) <- c(str1, str2, str3, str4, str5, str6, "")
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
    j3 <- c(1,2,3,8,9)
    j4 <- c(4,5,6,10,11,13,14)

    b[j1] <- lapply(b[j1], formatC, format = "f", digits = 1)
    b[j2] <- lapply(b[j2], formatC, format = "f", digits = 2)
    b[j3] <- lapply(b[j3], formatC, format = "f", digits = 3)
    b[j4] <- lapply(b[j4], formatC, format = "f", digits = 4)

    if (x$settings$typeBetaSpending != 'none' ||
        (k > 1 && any(x$byStageResults$futilityBounds[1:(k-1)] > -6))) {
      df <- t(b)
      rownames(df) <- c("Information rate",
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
      df <- t(b[,c(1,2,4,6,7,8,10,12)])
      rownames(df) <- c("Information rate",
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

    df <- t(b)

    if (!x$settings$normalApproximation) {
      rownames(df) <- c("Efficacy boundary (t)",
                        "Efficacy boundary (mean)",
                        "Efficacy boundary (p)")
    } else {
      rownames(df) <- c("Efficacy boundary (Z)",
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
  a <- x$overallResults
  s <- x$byStageResults
  k <- a$kMax

  if (k>1) {
    str1 <- paste0("Group-sequential design with ", k, " stages")
  } else {
    str1 <- "Fixed design"
  }

  str1 <- paste(str1, "for paired mean difference")

  str2 <- paste0("Paired difference under H0: ",
                 round(a$pairedDiffH0, 3), ", ",
                 "paired difference under H1: ",
                 round(a$pairedDiff, 3), ", ",
                 "standard deviation: ", round(a$stDev, 3))

  str3 <- paste0("Overall power: ",
                 round(a$overallReject, 4), ", ",
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
    asf <- tolower(x$settings$typeAlphaSpending)
    asfpar <- round(x$settings$parameterAlphaSpending, 3)
    asfuser <- round(x$settings$userAlphaSpending, 4)

    bsf <- tolower(x$settings$typeBetaSpending)
    bsfpar <- round(x$settings$parameterBetaSpending, 3)
    bsfuser <- round(x$settings$userBetaSpending, 4)

    if (asf == "of") {
      str7 <- paste0("Alpha spending: O'Brien-Fleming")
    } else if (asf == "p") {
      str7 <- paste0("Alpha spending: Pocock")
    } else if (asf == "wt") {
      str7 <- paste0("Alpha spending: Wang-Tsiatis(Delta = ", asfpar, ")")
    } else if (asf == "sfof") {
      str7 <- paste0("Alpha spending: Lan-DeMets O'Brien-Fleming")
    } else if (asf == "sfp") {
      str7 <- paste0("Alpha spending: Lan-DeMets Pocock")
    } else if (asf == "sfkd") {
      str7 <- paste0("Alpha spending: KD(rho = ", asfpar, ")")
    } else if (asf == "sfhsd") {
      str7 <- paste0("Alpha spending: HSD(gamma = ", asfpar, ")")
    } else if (asf == "user") {
      str7 <- paste0("Alpha spending: User defined(",
                     paste(asfuser, collapse = ","), ")")
    } else {
      str7 <- "Alpha spending: None"
    }

    if (bsf == "of") {
      str8 <- paste0("beta spending: O'Brien-Fleming")
    } else if (bsf == "p") {
      str8 <- paste0("beta spending: Pocock")
    } else if (bsf == "wt") {
      str8 <- paste0("beta spending: Wang-Tsiatis(Delta = ", bsfpar, ")")
    } else if (bsf == "sfof") {
      str8 <- paste0("beta spending: Lan-DeMets O'Brien-Fleming")
    } else if (bsf == "sfp") {
      str8 <- paste0("beta spending: Lan-DeMets Pocock")
    } else if (bsf == "sfkd") {
      str8 <- paste0("beta spending: KD(rho = ", bsfpar, ")")
    } else if (bsf == "sfhsd") {
      str8 <- paste0("beta spending: HSD(gamma = ", bsfpar, ")")
    } else if (bsf == "user") {
      str8 <- paste0("beta spending: User defined(",
                     paste(bsfuser, collapse = ","), ")")
    } else {
      str8 <- "beta spending: None"
    }

    if (!any(is.na(x$settings$spendingTime)) &&
        !all(x$settings$spendingTime == s$informationRates)) {
      str9 <- paste0("Spending time: ",
                     paste(round(x$settings$spendingTime, 3), collapse = ","))
      df1 <- data.frame(x = rep("", 9))
      colnames(df1) <- NULL
      rownames(df1) <- c(str1, str2, str3, str4, str5, str6,
                         paste(str7, str8, sep = ", "), str9, "")
    } else {
      df1 <- data.frame(x = rep("", 8))
      colnames(df1) <- NULL
      rownames(df1) <- c(str1, str2, str3, str4, str5, str6,
                         paste(str7, str8, sep = ", "), "")
    }
  } else {
    df1 <- data.frame(x = rep("", 7))
    colnames(df1) <- NULL
    rownames(df1) <- c(str1, str2, str3, str4, str5, str6, "")
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
    j3 <- c(1,2,3,8,9)
    j4 <- c(4,5,6,10,11,13,14)

    b[j1] <- lapply(b[j1], formatC, format = "f", digits = 1)
    b[j2] <- lapply(b[j2], formatC, format = "f", digits = 2)
    b[j3] <- lapply(b[j3], formatC, format = "f", digits = 3)
    b[j4] <- lapply(b[j4], formatC, format = "f", digits = 4)

    if (x$settings$typeBetaSpending != 'none' ||
        (k > 1 && any(x$byStageResults$futilityBounds[1:(k-1)] > -6))) {
      df <- t(b)
      rownames(df) <- c("Information rate",
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
      df <- t(b[,c(1,2,4,6,7,8,10,12)])
      rownames(df) <- c("Information rate",
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

    df <- t(b)

    if (!x$settings$normalApproximation) {
      rownames(df) <- c("Efficacy boundary (t)",
                        "Efficacy boundary (paired diff)",
                        "Efficacy boundary (p)")
    } else {
      rownames(df) <- c("Efficacy boundary (Z)",
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
  a <- x$overallResults
  s <- x$byStageResults
  k <- a$kMax

  if (k>1) {
    str1 <- paste0("Group-sequential design with ", k, " stages")
  } else {
    str1 <- "Fixed design"
  }

  str1 <- paste(str1, "for paired mean ratio")

  str2 <- paste0("Paired ratio under H0: ",
                 round(a$pairedRatioH0, 3), ", ",
                 "paired ratio under H1: ",
                 round(a$pairedRatio, 3), ", ",
                 "coefficient of variation: ", round(a$CV, 3))

  str3 <- paste0("Overall power: ",
                 round(a$overallReject, 4), ", ",
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
    asf <- tolower(x$settings$typeAlphaSpending)
    asfpar <- round(x$settings$parameterAlphaSpending, 3)
    asfuser <- round(x$settings$userAlphaSpending, 4)

    bsf <- tolower(x$settings$typeBetaSpending)
    bsfpar <- round(x$settings$parameterBetaSpending, 3)
    bsfuser <- round(x$settings$userBetaSpending, 4)

    if (asf == "of") {
      str7 <- paste0("Alpha spending: O'Brien-Fleming")
    } else if (asf == "p") {
      str7 <- paste0("Alpha spending: Pocock")
    } else if (asf == "wt") {
      str7 <- paste0("Alpha spending: Wang-Tsiatis(Delta = ", asfpar, ")")
    } else if (asf == "sfof") {
      str7 <- paste0("Alpha spending: Lan-DeMets O'Brien-Fleming")
    } else if (asf == "sfp") {
      str7 <- paste0("Alpha spending: Lan-DeMets Pocock")
    } else if (asf == "sfkd") {
      str7 <- paste0("Alpha spending: KD(rho = ", asfpar, ")")
    } else if (asf == "sfhsd") {
      str7 <- paste0("Alpha spending: HSD(gamma = ", asfpar, ")")
    } else if (asf == "user") {
      str7 <- paste0("Alpha spending: User defined(",
                     paste(asfuser, collapse = ","), ")")
    } else {
      str7 <- "Alpha spending: None"
    }

    if (bsf == "of") {
      str8 <- paste0("beta spending: O'Brien-Fleming")
    } else if (bsf == "p") {
      str8 <- paste0("beta spending: Pocock")
    } else if (bsf == "wt") {
      str8 <- paste0("beta spending: Wang-Tsiatis(Delta = ", bsfpar, ")")
    } else if (bsf == "sfof") {
      str8 <- paste0("beta spending: Lan-DeMets O'Brien-Fleming")
    } else if (bsf == "sfp") {
      str8 <- paste0("beta spending: Lan-DeMets Pocock")
    } else if (bsf == "sfkd") {
      str8 <- paste0("beta spending: KD(rho = ", bsfpar, ")")
    } else if (bsf == "sfhsd") {
      str8 <- paste0("beta spending: HSD(gamma = ", bsfpar, ")")
    } else if (bsf == "user") {
      str8 <- paste0("beta spending: User defined(",
                     paste(bsfuser, collapse = ","), ")")
    } else {
      str8 <- "beta spending: None"
    }

    if (!any(is.na(x$settings$spendingTime)) &&
        !all(x$settings$spendingTime == s$informationRates)) {
      str9 <- paste0("Spending time: ",
                     paste(round(x$settings$spendingTime, 3), collapse = ","))
      df1 <- data.frame(x = rep("", 9))
      colnames(df1) <- NULL
      rownames(df1) <- c(str1, str2, str3, str4, str5, str6,
                         paste(str7, str8, sep = ", "), str9, "")
    } else {
      df1 <- data.frame(x = rep("", 8))
      colnames(df1) <- NULL
      rownames(df1) <- c(str1, str2, str3, str4, str5, str6,
                         paste(str7, str8, sep = ", "), "")
    }
  } else {
    df1 <- data.frame(x = rep("", 7))
    colnames(df1) <- NULL
    rownames(df1) <- c(str1, str2, str3, str4, str5, str6, "")
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
    j3 <- c(1,2,3,8,9)
    j4 <- c(4,5,6,10,11,13,14)

    b[j1] <- lapply(b[j1], formatC, format = "f", digits = 1)
    b[j2] <- lapply(b[j2], formatC, format = "f", digits = 2)
    b[j3] <- lapply(b[j3], formatC, format = "f", digits = 3)
    b[j4] <- lapply(b[j4], formatC, format = "f", digits = 4)

    if (x$settings$typeBetaSpending != 'none' ||
        (k > 1 && any(x$byStageResults$futilityBounds[1:(k-1)] > -6))) {
      df <- t(b)
      rownames(df) <- c("Information rate",
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
      df <- t(b[,c(1,2,4,6,7,8,10,12)])
      rownames(df) <- c("Information rate",
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

    df <- t(b)

    if (!x$settings$normalApproximation) {
      rownames(df) <- c("Efficacy boundary (t)",
                        "Efficacy boundary (paired ratio)",
                        "Efficacy boundary (p)")
    } else {
      rownames(df) <- c("Efficacy boundary (Z)",
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
  a <- x$overallResults
  s <- x$byStageResults
  k <- a$kMax

  if (k>1) {
    str1 <- paste0("Group-sequential design with ", k, " stages")
  } else {
    str1 <- "Fixed design"
  }

  str1 <- paste(str1, "for two-sample mean difference")

  str2 <- paste0("Mean difference under H0: ",
                 round(a$meanDiffH0, 3), ", ",
                 "mean difference under H1: ",
                 round(a$meanDiff, 3), ", ",
                 "standard deviation: ", round(a$stDev, 3))

  str3 <- paste0("Overall power: ",
                 round(a$overallReject, 4), ", ",
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
      str9 = paste0("beta spending: Lan-DeMets Pocock")
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
               "efficacyMeanDiff", "futilityMeanDiff",
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
        (k > 1 && any(x$byStageResults$futilityBounds[1:(k-1)] > -6))) {
      df <- t(b)
      rownames(df) <- c("Information rate",
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
      df <- t(b[,c(1,2,4,6,7,8,10,12)])
      rownames(df) <- c("Information rate",
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

    df <- t(b)

    if (x$settings$normalApproximation) {
      rownames(df) <- c("Efficacy boundary (Z)",
                        "Efficacy boundary (mean diff)",
                        "Efficacy boundary (p)")
    } else {
      rownames(df) <- c("Efficacy boundary (t)",
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
  a <- x$overallResults
  s <- x$byStageResults
  k <- a$kMax

  if (k>1) {
    str1 <- paste0("Group-sequential design with ", k, " stages")
  } else {
    str1 <- "Fixed design"
  }

  str1 <- paste(str1, "for two-sample mean ratio")

  str2 <- paste0("Mean ratio under H0: ",
                 round(a$meanRatioH0, 3), ", ",
                 "mean ratio under H1: ",
                 round(a$meanRatio, 3), ", ",
                 "coefficient of variation: ", round(a$CV, 3))

  str3 <- paste0("Overall power: ",
                 round(a$overallReject, 4), ", ",
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
               "efficacyMeanRatio", "futilityMeanRatio",
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
        (k > 1 && any(x$byStageResults$futilityBounds[1:(k-1)] > -6))) {
      df <- t(b)
      rownames(df) <- c("Information rate",
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
      df <- t(b[,c(1,2,4,6,7,8,10,12)])
      rownames(df) <- c("Information rate",
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

    df <- t(b)

    if (x$settings$normalApproximation) {
      rownames(df) <- c("Efficacy boundary (Z)",
                        "Efficacy boundary (mean ratio)",
                        "Efficacy boundary (p)")
    } else {
      rownames(df) <- c("Efficacy boundary (t)",
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
  a <- x$overallResults
  s <- x$byStageResults
  k <- a$kMax

  if (k>1) {
    str1 <- paste0("Group-sequential design with ", k, " stages")
  } else {
    str1 <- "Fixed design"
  }

  str1 <- paste(str1, "for mean difference in 2x2 crossover")

  str2 <- paste0("Mean difference under H0: ",
                 round(a$meanDiffH0, 3), ", ",
                 "mean difference under H1: ",
                 round(a$meanDiff, 3), ", ",
                 "standard deviation: ",
                 round(a$stDev, 3))

  str3 <- paste0("Overall power: ",
                 round(a$overallReject, 4), ", ",
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
               "efficacyMeanDiff", "futilityMeanDiff",
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
        (k > 1 && any(x$byStageResults$futilityBounds[1:(k-1)] > -6))) {
      df <- t(b)
      rownames(df) <- c("Information rate",
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
      df <- t(b[,c(1,2,4,6,7,8,10,12)])
      rownames(df) <- c("Information rate",
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

    df <- t(b)

    if (x$settings$normalApproximation) {
      rownames(df) <- c("Efficacy boundary (Z)",
                        "Efficacy boundary (mean diff)",
                        "Efficacy boundary (p)")
    } else {
      rownames(df) <- c("Efficacy boundary (t)",
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
  a <- x$overallResults
  s <- x$byStageResults
  k <- a$kMax

  if (k>1) {
    str1 <- paste0("Group-sequential design with ", k, " stages")
  } else {
    str1 <- "Fixed design"
  }

  str1 <- paste(str1, "for mean ratio in 2x2 crossover")

  str2 <- paste0("Mean ratio under H0: ",
                 round(a$meanRatioH0, 3), ", ",
                 "mean ratio under H1: ",
                 round(a$meanRatio, 3), ", ",
                 "coefficient of variation: ", round(a$CV, 3))

  str3 <- paste0("Overall power: ",
                 round(a$overallReject, 4), ", ",
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
               "efficacyMeanRatio", "futilityMeanRatio",
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
        (k > 1 && any(x$byStageResults$futilityBounds[1:(k-1)] > -6))) {
      df <- t(b)
      rownames(df) <- c("Information rate",
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
      df <- t(b[,c(1,2,4,6,7,8,10,12)])
      rownames(df) <- c("Information rate",
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

    df <- t(b)

    if (!x$settings$normalApproximation) {
      rownames(df) <- c("Efficacy boundary (t)",
                        "Efficacy boundary (mean ratio)",
                        "Efficacy boundary (p)")
    } else {
      rownames(df) <- c("Efficacy boundary (Z)",
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
  a <- x$overallResults
  s <- x$byStageResults
  k <- a$kMax

  if (k>1) {
    str1 <- paste0("Group-sequential design with ", k, " stages")
  } else {
    str1 <- "Fixed design"
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
                 round(a$overallReject, 4), ", ",
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
    asf <- tolower(x$settings$typeAlphaSpending)
    asfpar <- round(x$settings$parameterAlphaSpending, 3)
    asfuser <- round(x$settings$userAlphaSpending, 4)

    if (asf == "of") {
      str7 <- paste0("Alpha spending: O'Brien-Fleming")
    } else if (asf == "p") {
      str7 <- paste0("Alpha spending: Pocock")
    } else if (asf == "wt") {
      str7 <- paste0("Alpha spending: Wang-Tsiatis(Delta = ", asfpar, ")")
    } else if (asf == "sfof") {
      str7 <- paste0("Alpha spending: Lan-DeMets O'Brien-Fleming")
    } else if (asf == "sfp") {
      str7 <- paste0("Alpha spending: Lan-DeMets Pocock")
    } else if (asf == "sfkd") {
      str7 <- paste0("Alpha spending: KD(rho = ", asfpar, ")")
    } else if (asf == "sfhsd") {
      str7 <- paste0("Alpha spending: HSD(gamma = ", asfpar, ")")
    } else if (asf == "user") {
      str7 <- paste0("Alpha spending: User defined(",
                     paste(asfuser, collapse = ","), ")")
    } else {
      str7 <- "Alpha spending: None"
    }

    if (!any(is.na(x$settings$spendingTime)) &&
        !all(x$settings$spendingTime == s$informationRates)) {
      str8 <- paste0("Spending time: ",
                     paste(round(x$settings$spendingTime, 3), collapse = ","))
      df1 <- data.frame(x = rep("", 9))
      colnames(df1) <- NULL
      rownames(df1) <- c(str1, str2, str3, str4, str5, str6, str7, str8, "")
    } else {
      df1 <- data.frame(x = rep("", 8))
      colnames(df1) <- NULL
      rownames(df1) <- c(str1, str2, str3, str4, str5, str6, str7, "")
    }
  } else {
    df1 <- data.frame(x = rep("", 7))
    colnames(df1) <- NULL
    rownames(df1) <- c(str1, str2, str3, str4, str5, str6, "")
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
    j3 <- c(1,2,7,8)
    j4 <- c(3,4,5,9)

    b[j1] <- lapply(b[j1], formatC, format = "f", digits = 1)
    b[j2] <- lapply(b[j2], formatC, format = "f", digits = 2)
    b[j3] <- lapply(b[j3], formatC, format = "f", digits = 3)
    b[j4] <- lapply(b[j4], formatC, format = "f", digits = 4)

    df <- t(b)
    rownames(df) <- c("Information rate",
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

    df <- t(b)

    if (x$settings$normalApproximation) {
      rownames(df) <- c("Boundary for each 1-sided test (Z)",
                        "Boundary for lower limit (paired diff)",
                        "Boundary for upper limit (paired diff)",
                        "Boundary for each 1-sided test (p)")
    } else {
      rownames(df) <- c("Boundary for each 1-sided test (t)",
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
  a <- x$overallResults
  s <- x$byStageResults
  k <- a$kMax

  if (k>1) {
    str1 <- paste0("Group-sequential design with ", k, " stages")
  } else {
    str1 <- "Fixed design"
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
                 round(a$overallReject, 4), ", ",
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
    asf <- tolower(x$settings$typeAlphaSpending)
    asfpar <- round(x$settings$parameterAlphaSpending, 3)
    asfuser <- round(x$settings$userAlphaSpending, 4)

    if (asf == "of") {
      str7 <- paste0("Alpha spending: O'Brien-Fleming")
    } else if (asf == "p") {
      str7 <- paste0("Alpha spending: Pocock")
    } else if (asf == "wt") {
      str7 <- paste0("Alpha spending: Wang-Tsiatis(Delta = ", asfpar, ")")
    } else if (asf == "sfof") {
      str7 <- paste0("Alpha spending: Lan-DeMets O'Brien-Fleming")
    } else if (asf == "sfp") {
      str7 <- paste0("Alpha spending: Lan-DeMets Pocock")
    } else if (asf == "sfkd") {
      str7 <- paste0("Alpha spending: KD(rho = ", asfpar, ")")
    } else if (asf == "sfhsd") {
      str7 <- paste0("Alpha spending: HSD(gamma = ", asfpar, ")")
    } else if (asf == "user") {
      str7 <- paste0("Alpha spending: User defined(",
                     paste(asfuser, collapse = ","), ")")
    } else {
      str7 <- "Alpha spending: None"
    }

    if (!any(is.na(x$settings$spendingTime)) &&
        !all(x$settings$spendingTime == s$informationRates)) {
      str8 <- paste0("Spending time: ",
                     paste(round(x$settings$spendingTime, 3), collapse = ","))
      df1 <- data.frame(x = rep("", 9))
      colnames(df1) <- NULL
      rownames(df1) <- c(str1, str2, str3, str4, str5, str6, str7, str8, "")
    } else {
      df1 <- data.frame(x = rep("", 8))
      colnames(df1) <- NULL
      rownames(df1) <- c(str1, str2, str3, str4, str5, str6, str7, "")
    }
  } else {
    df1 <- data.frame(x = rep("", 7))
    colnames(df1) <- NULL
    rownames(df1) <- c(str1, str2, str3, str4, str5, str6, "")
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
    j3 <- c(1,2,7,8)
    j4 <- c(3,4,5,9)

    b[j1] <- lapply(b[j1], formatC, format = "f", digits = 1)
    b[j2] <- lapply(b[j2], formatC, format = "f", digits = 2)
    b[j3] <- lapply(b[j3], formatC, format = "f", digits = 3)
    b[j4] <- lapply(b[j4], formatC, format = "f", digits = 4)

    df <- t(b)
    rownames(df) <- c("Information rate",
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

    df <- t(b)

    if (x$settings$normalApproximation) {
      rownames(df) <- c("Boundary for each 1-sided test (Z)",
                        "Boundary for lower limit (paired ratio)",
                        "Boundary for upper limit (paired ratio)",
                        "Boundary for each 1-sided test (p)")
    } else {
      rownames(df) <- c("Boundary for each 1-sided test (t)",
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
  a <- x$overallResults
  s <- x$byStageResults
  k <- a$kMax

  if (k>1) {
    str1 <- paste0("Group-sequential design with ", k, " stages")
  } else {
    str1 <- "Fixed design"
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
                 round(a$overallReject, 4), ", ",
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
               "cumulativeAttainedAlpha", "numberOfSubjects",
               "efficacyMeanDiffLower", "efficacyMeanDiffUpper",
               "efficacyP", "information")]

    # format number of digits after decimal for each column
    j1 <- 6
    j2 <- 10
    j3 <- c(1,2,7,8)
    j4 <- c(3,4,5,9)

    b[j1] <- lapply(b[j1], formatC, format = "f", digits = 1)
    b[j2] <- lapply(b[j2], formatC, format = "f", digits = 2)
    b[j3] <- lapply(b[j3], formatC, format = "f", digits = 3)
    b[j4] <- lapply(b[j4], formatC, format = "f", digits = 4)

    df <- t(b)
    rownames(df) <- c("Information rate",
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

    df <- t(b)

    if (x$settings$normalApproximation) {
      rownames(df) <- c("Boundary for each 1-sided test (Z)",
                        "Boundary for lower limit (mean diff)",
                        "Boundary for upper limit (mean diff)",
                        "Boundary for each 1-sided test (p)")
    } else {
      rownames(df) <- c("Boundary for each 1-sided test (t)",
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
  a <- x$overallResults
  s <- x$byStageResults
  k <- a$kMax

  if (k>1) {
    str1 <- paste0("Group-sequential design with ", k, " stages")
  } else {
    str1 <- "Fixed design"
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
                 round(a$overallReject, 4), ", ",
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
               "cumulativeAttainedAlpha", "numberOfSubjects",
               "efficacyMeanRatioLower", "efficacyMeanRatioUpper",
               "efficacyP", "information")]

    # format number of digits after decimal for each column
    j1 <- 6
    j2 <- 10
    j3 <- c(1,2,7,8)
    j4 <- c(3,4,5,9)

    b[j1] <- lapply(b[j1], formatC, format = "f", digits = 1)
    b[j2] <- lapply(b[j2], formatC, format = "f", digits = 2)
    b[j3] <- lapply(b[j3], formatC, format = "f", digits = 3)
    b[j4] <- lapply(b[j4], formatC, format = "f", digits = 4)

    df <- t(b)
    rownames(df) <- c("Information rate",
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

    df <- t(b)

    if (x$settings$normalApproximation) {
      rownames(df) <- c("Boundary for each 1-sided test (Z)",
                        "Boundary for lower limit (mean ratio)",
                        "Boundary for upper limit (mean ratio)",
                        "Boundary for each 1-sided test (p)")
    } else {
      rownames(df) <- c("Boundary for each 1-sided test (t)",
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
  a <- x$overallResults
  s <- x$byStageResults
  k <- a$kMax

  if (k>1) {
    str1 <- paste0("Group-sequential design with ", k, " stages")
  } else {
    str1 <- "Fixed design"
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
                 round(a$overallReject, 4), ", ",
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
      rownames(df1) <- c(str1, str2, str3, str4, str5, str6, str7,
                         str8, "")
    }
  } else {
    df1 <- data.frame(x = rep("", 8))
    colnames(df1) <- NULL
    rownames(df1) <- c(str1, str2, str3, str4, str5, str6, str7, "")
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
    j3 <- c(1,2,7,8)
    j4 <- c(3,4,5,9)

    b[j1] <- lapply(b[j1], formatC, format = "f", digits = 1)
    b[j2] <- lapply(b[j2], formatC, format = "f", digits = 2)
    b[j3] <- lapply(b[j3], formatC, format = "f", digits = 3)
    b[j4] <- lapply(b[j4], formatC, format = "f", digits = 4)

    df <- t(b)
    rownames(df) <- c("Information rate",
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

    df <- t(b)

    if (x$settings$normalApproximation) {
      rownames(df) <- c("Boundary for each 1-sided test (Z)",
                        "Boundary for lower limit (mean diff)",
                        "Boundary for upper limit (mean diff)",
                        "Boundary for each 1-sided test (p)")
    } else {
      rownames(df) <- c("Boundary for each 1-sided test (t)",
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
  a <- x$overallResults
  s <- x$byStageResults
  k <- a$kMax

  if (k>1) {
    str1 <- paste0("Group-sequential design with ", k, " stages")
  } else {
    str1 <- "Fixed design"
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
                 round(a$overallReject, 4), ", ",
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
               "cumulativeAttainedAlpha", "numberOfSubjects",
               "efficacyMeanRatioLower", "efficacyMeanRatioUpper",
               "efficacyP", "information")]

    # format number of digits after decimal for each column
    j1 <- 6
    j2 <- 10
    j3 <- c(1,2,7,8)
    j4 <- c(3,4,5,9)

    b[j1] <- lapply(b[j1], formatC, format = "f", digits = 1)
    b[j2] <- lapply(b[j2], formatC, format = "f", digits = 2)
    b[j3] <- lapply(b[j3], formatC, format = "f", digits = 3)
    b[j4] <- lapply(b[j4], formatC, format = "f", digits = 4)

    df <- t(b)
    rownames(df) <- c("Information rate",
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

    df <- t(b)

    if (x$settings$normalApproximation) {
      rownames(df) <- c("Boundary for each 1-sided test (Z)",
                        "Boundary for lower limit (mean ratio)",
                        "Boundary for upper limit (mean ratio)",
                        "Boundary for each 1-sided test (p)")
    } else {
      rownames(df) <- c("Boundary for each 1-sided test (t)",
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
  a <- x$overallResults
  s <- x$byStageResults
  k <- a$kMax

  if (k>1) {
    str1 <- paste0("Group-sequential design with ", k, " stages")
  } else {
    str1 <- "Fixed design"
  }

  str1 <- paste(str1, "for two-sample Wilcoxon test")

  str2 <- paste0("Probability of observations in treatment ",
                 "larger than those in control: ",
                 round(a$pLarger, 3))

  str3 <- paste0("Overall power: ",
                 round(a$overallReject, 4), ", ",
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
               "efficacyPLarger", "futilityPLarger",
               "efficacyP", "futilityP", "information",
               "cumulativeRejectionH0", "cumulativeFutilityH0")]

    # format number of digits after decimal for each column
    j1 <- 7
    j2 <- 12
    j3 <- c(1,2,3)
    j4 <- c(4,5,6,8,9,10,11,13,14)

    b[j1] <- lapply(b[j1], formatC, format = "f", digits = 1)
    b[j2] <- lapply(b[j2], formatC, format = "f", digits = 2)
    b[j3] <- lapply(b[j3], formatC, format = "f", digits = 3)
    b[j4] <- lapply(b[j4], formatC, format = "f", digits = 4)

    if (x$settings$typeBetaSpending != 'none' ||
        (k > 1 && any(x$byStageResults$futilityBounds[1:(k-1)] > -6))) {
      df <- t(b)
      rownames(df) <- c("Information rate",
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
      df <- t(b[,c(1,2,4,6,7,8,10,12)])
      rownames(df) <- c("Information rate",
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

    df <- t(b)

    rownames(df) <- c("Efficacy boundary (Z)",
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
  a <- x$overallResults
  s <- x$byStageResults
  kMax <- a$kMax

  if (kMax>1) {
    str1 <- paste0("Group-sequential design with ", kMax, " stages")
  } else {
    str1 <- "Fixed design"
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
                 round(a$overallReject, 4), ", ",
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
      df1 <- data.frame(x = rep("", 13))
      colnames(df1) <- NULL
      rownames(df1) <- c(str1, str2, str3, str4, str5, str6, str7, str8,
                         str9, str10, paste(str11, str12, sep = ", "),
                         str13, "")
    } else {
      df1 <- data.frame(x = rep("", 12))
      colnames(df1) <- NULL
      rownames(df1) <- c(str1, str2, str3, str4, str5, str6, str7, str8,
                         str9, str10, paste(str11, str12, sep = ", "),
                         "")
    }
  } else {
    df1 <- data.frame(x = rep("", 11))
    colnames(df1) <- NULL
    rownames(df1) <- c(str1, str2, str3, str4, str5, str6, str7, str8,
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
    j3 <- c(1,2,3,10,11)
    j4 <- c(4,5,6,12,13,15,16)

    b[j1] <- lapply(b[j1], formatC, format = "f", digits = 1)
    b[j2] <- lapply(b[j2], formatC, format = "f", digits = 2)
    b[j3] <- lapply(b[j3], formatC, format = "f", digits = 3)
    b[j4] <- lapply(b[j4], formatC, format = "f", digits = 4)

    if (x$settings$typeBetaSpending != 'none' ||
        (kMax > 1 && any(s$futilityBounds[1:(kMax-1)] > -6))) {
      df <- t(b)
      rownames(df) <- c("Information rate",
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
      df <- t(b[,c(1,2,4,6,7,8,9,10,12,14)])
      rownames(df) <- c("Information rate",
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

    df <- t(b)

    if (x$settings$normalApproximation) {
      rownames(df) <- c("Efficacy boundary (Z)",
                        "Efficacy boundary (mean diff)",
                        "Efficacy boundary (p)")
    } else {
      rownames(df) <- c("Efficacy boundary (t)",
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
  carryover <- x$carryover

  if (carryover) {
    str1 <- paste("Testing direct treatment effects accounting for",
                  "carryover effects")
  } else {
    str1 <- paste("Testing direct treatment effects without accounting for",
                  "carryover effects")
  }

  df1 <- data.frame(x = "")
  colnames(df1) <- NULL
  rownames(df1) <- str1

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

  str5 <- paste0("Power: ", round(x$power, 4), ", ",
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
    df3 <- data.frame(x = rep("", 11))
    colnames(df3) <- NULL
    rownames(df3) <- c(str2, str3, str4, str5, str6, str7, str8,
                       str9, str10, str11, str12)

  } else {
    df3 <- data.frame(x = rep("", 9))
    colnames(df3) <- NULL
    rownames(df3) <- c(str2, str3, str4, str5, str6, str7,
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
  carryover <- x$carryover

  if (carryover) {
    str1 <- paste("Testing equivalence in direct treatment effects",
                  "accounting for carryover effects")
  } else {
    str1 <- paste("Testing equivalence in direct treatment effects",
                  "without accounting for carryover effects")
  }

  df1 <- data.frame(x = "")
  colnames(df1) <- NULL
  rownames(df1) <- str1

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

  str6 <- paste0("Power: ", round(x$power, 4), ", ",
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
    df3 <- data.frame(x = rep("", 12))
    colnames(df3) <- NULL
    rownames(df3) <- c(str2, str3, str4, str5, str6, str7, str8,
                       str9, str10, str11, str12, str13)

  } else {
    df3 <- data.frame(x = rep("", 10))
    colnames(df3) <- NULL
    rownames(df3) <- c(str2, str3, str4, str5, str6, str7, str8,
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
  df1 <- data.frame(alpha = x$alpha, power = x$power,
                    n = x$n, ngroups = x$ngroups,
                    stDev = x$stDev, effectsize = x$effectsize)
  rownames(df1) <- NULL
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
  df1 <- x$powerdf
  df1$alpha <- x$alpha
  df1$nlevelsA <- x$nlevelsA
  df1$nlevelsB <- x$nlevelsB
  df1$stDev <- x$stDev
  df1$effectsizeA <- x$effectsizeA
  df1$effectsizeB <- x$effectsizeB
  df1$effectsizeAB <- x$effectsizeAB
  rownames(df1) <- NULL
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
  df1 <- data.frame(alpha = x$alpha, power = x$power,
                    n = x$n, ngroups = x$ngroups, stDev = x$stDev,
                    meanContrastH0 = x$meanContrastH0,
                    meanContrast = x$meanContrast,
                    effectsize = x$effectsize)
  rownames(df1) <- NULL
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
  df1 <- data.frame(alpha = x$alpha, power = x$power,
                    n = x$n, ngroups = x$ngroups,
                    stDev = x$stDev, corr = x$corr,
                    effectsize = x$effectsize)
  rownames(df1) <- NULL
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
  df1 <- data.frame(alpha = x$alpha, power = x$power,
                    n = x$n, ngroups = x$ngroups, stDev = x$stDev,
                    corr = x$corr,
                    meanContrastH0 = x$meanContrastH0,
                    meanContrast = x$meanContrast,
                    effectsize = x$effectsize)
  rownames(df1) <- NULL
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
  a <- x$overallResults
  s <- x$byStageResults
  k <- a$kMax

  if (k>1) {
    str1 <- paste0("Group-sequential design with ", k, " stages")
  } else {
    str1 <- "Fixed design"
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
                 round(a$overallReject, 4), ", ",
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
               "efficacySlope", "futilitySlope",
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
        (k > 1 && any(x$byStageResults$futilityBounds[1:(k-1)] > -6))) {
      df <- t(b)
      rownames(df) <- c("Information rate",
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
      df <- t(b[,c(1,2,4,6,7,8,10,12)])
      rownames(df) <- c("Information rate",
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

    df <- t(b)

    if (!x$settings$normalApproximation) {
      rownames(df) <- c("Efficacy boundary (t)",
                        "Efficacy boundary (slope)",
                        "Efficacy boundary (p)")
    } else {
      rownames(df) <- c("Efficacy boundary (Z)",
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
  a <- x$overallResults
  s <- x$byStageResults
  k <- a$kMax

  if (k>1) {
    str1 <- paste0("Group-sequential design with ", k, " stages")
  } else {
    str1 <- "Fixed design"
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
                 round(a$overallReject, 4), ", ",
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
    asf <- tolower(x$settings$typeAlphaSpending)
    asfpar <- round(x$settings$parameterAlphaSpending, 3)
    asfuser <- round(x$settings$userAlphaSpending, 4)

    bsf <- tolower(x$settings$typeBetaSpending)
    bsfpar <- round(x$settings$parameterBetaSpending, 3)
    bsfuser <- round(x$settings$userBetaSpending, 4)

    if (asf == "of") {
      str9 <- paste0("Alpha spending: O'Brien-Fleming")
    } else if (asf == "p") {
      str9 <- paste0("Alpha spending: Pocock")
    } else if (asf == "wt") {
      str9 <- paste0("Alpha spending: Wang-Tsiatis(Delta = ", asfpar, ")")
    } else if (asf == "sfof") {
      str9 <- paste0("Alpha spending: Lan-DeMets O'Brien-Fleming")
    } else if (asf == "sfp") {
      str9 <- paste0("Alpha spending: Lan-DeMets Pocock")
    } else if (asf == "sfkd") {
      str9 <- paste0("Alpha spending: KD(rho = ", asfpar, ")")
    } else if (asf == "sfhsd") {
      str9 <- paste0("Alpha spending: HSD(gamma = ", asfpar, ")")
    } else if (asf == "user") {
      str9 <- paste0("Alpha spending: User defined(",
                     paste(asfuser, collapse = ","), ")")
    } else {
      str9 <- "Alpha spending: None"
    }

    if (bsf == "of") {
      str10 <- paste0("beta spending: O'Brien-Fleming")
    } else if (bsf == "p") {
      str10 <- paste0("beta spending: Pocock")
    } else if (bsf == "wt") {
      str10 <- paste0("beta spending: Wang-Tsiatis(Delta = ", bsfpar, ")")
    } else if (bsf == "sfof") {
      str10 <- paste0("beta spending: Lan-DeMets O'Brien-Fleming")
    } else if (bsf == "sfp") {
      str10 <- paste0("beta spending: Lan-DeMets Pocock")
    } else if (bsf == "sfkd") {
      str10 <- paste0("beta spending: KD(rho = ", bsfpar, ")")
    } else if (bsf == "sfhsd") {
      str10 <- paste0("beta spending: HSD(gamma = ", bsfpar, ")")
    } else if (bsf == "user") {
      str10 <- paste0("beta spending: User defined(",
                      paste(bsfuser, collapse = ","), ")")
    } else {
      str10 <- "beta spending: None"
    }

    if (!any(is.na(x$settings$spendingTime)) &&
        !all(x$settings$spendingTime == s$informationRates)) {
      str11 <- paste0("Spending time: ",
                      paste(round(x$settings$spendingTime, 3), collapse = ","))
      df1 <- data.frame(x = rep("", 11))
      colnames(df1) <- NULL
      rownames(df1) <- c(str1, str2, str3, str4, str5, str6, str7, str8,
                         paste(str9, str10, sep = ", "), str11, "")
    } else {
      df1 <- data.frame(x = rep("", 10))
      colnames(df1) <- NULL
      rownames(df1) <- c(str1, str2, str3, str4, str5, str6, str7, str8,
                         paste(str9, str10, sep = ", "), "")
    }
  } else {
    df1 <- data.frame(x = rep("", 9))
    colnames(df1) <- NULL
    rownames(df1) <- c(str1, str2, str3, str4, str5, str6, str7, str8, "")
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
    j3 <- c(1,2,3,8,9)
    j4 <- c(4,5,6,10,11,13,14)

    b[j1] <- lapply(b[j1], formatC, format = "f", digits = 1)
    b[j2] <- lapply(b[j2], formatC, format = "f", digits = 2)
    b[j3] <- lapply(b[j3], formatC, format = "f", digits = 3)
    b[j4] <- lapply(b[j4], formatC, format = "f", digits = 4)

    if (x$settings$typeBetaSpending != 'none' ||
        (k > 1 && any(x$byStageResults$futilityBounds[1:(k-1)] > -6))) {
      df <- t(b)
      rownames(df) <- c("Information rate",
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
      df <- t(b[,c(1,2,4,6,7,8,10,12)])
      rownames(df) <- c("Information rate",
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

    df <- t(b)

    if (!x$settings$normalApproximation) {
      rownames(df) <- c("Efficacy boundary (t)",
                        "Efficacy boundary (slope diff)",
                        "Efficacy boundary (p)")
    } else {
      rownames(df) <- c("Efficacy boundary (Z)",
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
  a <- x$overallResults
  s <- x$byStageResults
  kMax <- a$kMax

  if (kMax>1) {
    str1 <- paste0("Group-sequential design with ", kMax, " stages")
  } else {
    str1 <- "Fixed design"
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
                 round(a$overallReject, 4), ", ",
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
    asf <- tolower(x$settings$typeAlphaSpending)
    asfpar <- round(x$settings$parameterAlphaSpending, 3)
    asfuser <- round(x$settings$userAlphaSpending, 4)

    bsf <- tolower(x$settings$typeBetaSpending)
    bsfpar <- round(x$settings$parameterBetaSpending, 3)
    bsfuser <- round(x$settings$userBetaSpending, 4)

    if (asf == "of") {
      str12 <- paste0("Alpha spending: O'Brien-Fleming")
    } else if (asf == "p") {
      str12 <- paste0("Alpha spending: Pocock")
    } else if (asf == "wt") {
      str12 <- paste0("Alpha spending: Wang-Tsiatis(Delta = ", asfpar, ")")
    } else if (asf == "sfof") {
      str12 = paste0("Alpha spending: Lan-DeMets O'Brien-Fleming")
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
      str13 <- paste0("beta spending: User defined(",
                      paste(bsfuser, collapse = ","), ")")
    } else {
      str13 <- "beta spending: None"
    }

    if (!any(is.na(x$settings$spendingTime)) &&
        !all(x$settings$spendingTime == s$informationRates)) {
      str14 <- paste0("Spending time: ",
                      paste(round(x$settings$spendingTime, 3), collapse = ","))
      df1 <- data.frame(x = rep("", 14))
      colnames(df1) <- NULL
      rownames(df1) <- c(str1, str2, str3, str4, str5, str6, str7, str8,
                         str9, str10, str11, paste(str12, str13, sep = ", "),
                         str14, "")
    } else {
      df1 <- data.frame(x = rep("", 13))
      colnames(df1) <- NULL
      rownames(df1) <- c(str1, str2, str3, str4, str5, str6, str7, str8,
                         str9, str10, str11, paste(str12, str13, sep = ", "),
                         "")
    }
  } else {
    df1 <- data.frame(x = rep("", 12))
    colnames(df1) <- NULL
    rownames(df1) <- c(str1, str2, str3, str4, str5, str6, str7, str8,
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
      j3 <- c(1,2,3,9,10)
      j4 <- c(4,5,6,11,12,14,15)

      b[j1] <- lapply(b[j1], formatC, format = "f", digits = 1)
      b[j2] <- lapply(b[j2], formatC, format = "f", digits = 2)
      b[j3] <- lapply(b[j3], formatC, format = "f", digits = 3)
      b[j4] <- lapply(b[j4], formatC, format = "f", digits = 4)

      if (x$settings$typeBetaSpending != 'none' ||
          (kMax > 1 && any(s$futilityBounds[1:(kMax-1)] > -6))) {
        df <- t(b)
        rownames(df) <- c("Information rate",
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
        df <- t(b[,c(1,2,4,6,7,8,9,11,13)])
        rownames(df) <- c("Information rate",
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
      j2 <- 14
      j3 <- c(1,2,3,10,11)
      j4 <- c(4,5,6,12,13,15,16)

      b[j1] <- lapply(b[j1], formatC, format = "f", digits = 1)
      b[j2] <- lapply(b[j2], formatC, format = "f", digits = 2)
      b[j3] <- lapply(b[j3], formatC, format = "f", digits = 3)
      b[j4] <- lapply(b[j4], formatC, format = "f", digits = 4)

      if (x$settings$typeBetaSpending != 'none' ||
          (kMax > 1 && any(s$futilityBounds[1:(kMax-1)] > -6))) {
        df <- t(b)
        rownames(df) <- c("Information rate",
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
        df <- t(b[,c(1,2,4,6,7,8,9,10,12,14)])
        rownames(df) <- c("Information rate",
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

    df <- t(b)

    if (x$settings$normalApproximation) {
      rownames(df) <- c("Efficacy boundary (Z)",
                        "Efficacy boundary (slope diff)",
                        "Efficacy boundary (p)")
    } else {
      rownames(df) <- c("Efficacy boundary (t)",
                        "Efficacy boundary (slope diff)",
                        "Efficacy boundary (p)")
    }

    colnames(df) <- NA
  }

  print(df1, ..., na.print = "" , quote = FALSE )
  print(df, ..., na.print = "" , quote = FALSE )
  invisible(x)
}
