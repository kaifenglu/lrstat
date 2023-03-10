#' @title Print group sequential design
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
#' @export
print.design <- function(x, ...) {
  a = x$overallResults;
  s = x$byStageResults;
  k = a$kMax;
  
  if (k>1) {
    str1 = paste0("Group-sequential trial with ", k, " stages")
  } else {
    str1 = "Fixed design"
  }
  
  
  str2 <- paste0("Overall power: ",
                 round(a$overallReject, 3), ", ",
                 "overall significance level (1-sided): ",
                 round(a$alpha, 4))
  
  str3 <- paste0("Drift parameter: ", round(a$drift, 3))
  str4 <- paste0("Inflation factor: ", round(a$inflationFactor, 3))
  
  
  df1 = data.frame(x = rep("", 5))
  colnames(df1) = NULL
  row.names(df1) = c(str1, str2, str3, str4, "")
  
  
  if (k>1) {
    b <- s[, c("informationRates", "efficacyBounds", "futilityBounds",
               "cumulativeRejection", "cumulativeFutility",
               "cumulativeAlphaSpent", 
               "efficacyP", "futilityP")]
    
    # format number of digits after decimal for each column
    j3 <- c(1,2,3,4,5)
    j4 <- c(6,7,8)
    
    b[j3] <- lapply(b[j3], formatC, format = "f", digits = 3)
    b[j4] <- lapply(b[j4], formatC, format = "f", digits = 4)
    
    if (x$settings$typeBetaSpending != 'none') {
      df = t(b)
      rownames(df) = c("Information rate",
                       "Efficacy boundary (Z-scale)",
                       "Futility boundary (Z-scale)",
                       "Cumulative rejection",
                       "Cumulative futility",
                       "Cumulative alpha spent",
                       "Efficacy boundary (p-scale)",
                       "Futility boundary (p-scale)")
      
    } else {
      df = t(b[,c(1,2,4,6,7)])
      rownames(df) = c("Information rate",
                       "Efficacy boundary (Z-scale)",
                       "Cumulative rejection",
                       "Cumulative alpha spent",
                       "Efficacy boundary (p-scale)")
    }
    
    
    colnames(df) <- paste("Stage", seq_len(ncol(df)), sep=" ")
  } else {
    b <- s[, c("efficacyBounds", "efficacyP")]
    
    # format number of digits after decimal for each column
    j3 <- 1
    j4 <- 2
    
    b[j3] <- lapply(b[j3], formatC, format = "f", digits = 3)
    b[j4] <- lapply(b[j4], formatC, format = "f", digits = 4)      
    
    df = t(b)
    
    rownames(df) = c("Efficacy boundary (Z-scale)",
                     "Efficacy boundary (p-scale)")
    
    colnames(df) <- NA
  }
  
  print(df1, ..., na.print = "" , quote = FALSE )
  print(df, ..., na.print = "" , quote = FALSE )
  invisible(x)
}


#' @title Print power and sample size results
#' @description Prints the summary statistics from power calculation.
#'
#' @param x The lrpower object to print.
#' @param ... Ensures that all arguments starting from "..." are named.
#'
#' @return A tabular printout of the summary statistics from power calculation.
#'
#' @keywords internal
#'
#' @export
print.lrpower <- function(x, ...) {
  a = x$overallResults;
  s = x$byStageResults;
  k = a$kMax;
  
  if (k>1) {
    str1 = paste0("Group-sequential trial with ", k, " stages")
  } else {
    str1 = "Fixed design"
  }
  
  if (a$rho1 != 0 || a$rho2 != 0) {
    str1 = paste0(str1, ", FH(", a$rho1, ", ", a$rho2, ")")
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
    
    str6 <- paste0("Total study duration: ",
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
    
    str6 <- paste0("Study duration: ",
                   round(a$studyDuration, 1))
  }
  
  str7 <- paste0("Accrual duration: ",
                 round(a$accrualDuration, 1), ", ",
                 "follow-up duration: ",
                 round(a$followupTime, 1), ", ",
                 "fixed follow-up: ", a$fixedFollowup)
  
  df1 = data.frame(x = rep("", 8))
  colnames(df1) = NULL
  row.names(df1) = c(str1, str2, str3, str4, str5, str6, str7, "")
  
  
  if (k>1) {
    if (a$estimateHazardRatio) {
      b <- s[, c("informationRates", "efficacyBounds", "futilityBounds",
                 "cumulativeRejection", "cumulativeFutility",
                 "cumulativeAlphaSpent", 
                 "numberOfEvents", "numberOfDropouts", "numberOfSubjects", 
                 "analysisTime", "efficacyHR", "futilityHR", 
                 "efficacyP", "futilityP", "information", "HR")]
      
      # format number of digits after decimal for each column
      j1 <- c(7,8,9,10)
      j2 <- 15
      j3 <- c(1,2,3,4,5,11,12,16)
      j4 <- c(6,13,14)
      
      b[j1] <- lapply(b[j1], formatC, format = "f", digits = 1)
      b[j2] <- lapply(b[j2], formatC, format = "f", digits = 2)
      b[j3] <- lapply(b[j3], formatC, format = "f", digits = 3)
      b[j4] <- lapply(b[j4], formatC, format = "f", digits = 4)
      
      if (x$settings$typeBetaSpending != 'none') {
        df = t(b)
        rownames(df) = c("Information rate",
                         "Efficacy boundary (Z-scale)",
                         "Futility boundary (Z-scale)",
                         "Cumulative rejection",
                         "Cumulative futility",
                         "Cumulative alpha spent",
                         "Number of events",
                         "Number of dropouts",
                         "Number of subjects",
                         "Analysis time",
                         "Efficacy boundary (HR-scale)",
                         "Futility boundary (HR-scale)",
                         "Efficacy boundary (p-scale)",
                         "Futility boundary (p-scale)",
                         "Information",
                         "HR")
        
      } else {
        df = t(b[,c(1,2,4,6,7,8,9,10,11,13,15,16)])
        rownames(df) = c("Information rate",
                         "Efficacy boundary (Z-scale)",
                         "Cumulative rejection",
                         "Cumulative alpha spent",
                         "Number of events",
                         "Number of dropouts",
                         "Number of subjects",
                         "Analysis time",
                         "Efficacy boundary (HR-scale)",
                         "Efficacy boundary (p-scale)",
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
      j3 <- c(1,2,3,4,5)
      j4 <- c(6,11,12)
      
      b[j1] <- lapply(b[j1], formatC, format = "f", digits = 1)
      b[j2] <- lapply(b[j2], formatC, format = "f", digits = 2)
      b[j3] <- lapply(b[j3], formatC, format = "f", digits = 3)
      b[j4] <- lapply(b[j4], formatC, format = "f", digits = 4)
      
      if (x$settings$typeBetaSpending != 'none') {
        df = t(b)
        rownames(df) = c("Information rate",
                         "Efficacy boundary (Z-scale)",
                         "Futility boundary (Z-scale)",
                         "Cumulative rejection",
                         "Cumulative futility",
                         "Cumulative alpha spent",
                         "Number of events",
                         "Number of dropouts",
                         "Number of subjects",
                         "Analysis time",
                         "Efficacy boundary (p-scale)",
                         "Futility boundary (p-scale)",
                         "Information")
        
      } else {
        df = t(b[,c(1,2,4,6,7,8,9,10,11,13)])
        rownames(df) = c("Information rate",
                         "Efficacy boundary (Z-scale)",
                         "Cumulative rejection",
                         "Cumulative alpha spent",
                         "Number of events",
                         "Number of dropouts",
                         "Number of subjects",
                         "Analysis time",
                         "Efficacy boundary (p-scale)",
                         "Information")
      }
      
    }
    
    colnames(df) <- paste("Stage", seq_len(ncol(df)), sep=" ")
  } else {
    if (a$estimateHazardRatio) {
      b <- s[, c("efficacyBounds", "efficacyHR", "efficacyP", 
                 "information", "HR")]
      
      # format number of digits after decimal for each column
      j2 <- 4
      j3 <- c(1,2,5)
      j4 <- 3
      
      b[j2] <- lapply(b[j2], formatC, format = "f", digits = 2)
      b[j3] <- lapply(b[j3], formatC, format = "f", digits = 3)
      b[j4] <- lapply(b[j4], formatC, format = "f", digits = 4)      
      
      df = t(b)
      
      rownames(df) = c("Efficacy boundary (Z-scale)",
                       "Efficacy boundary (HR-scale)",
                       "Efficacy boundary (p-scale)",
                       "Information",
                       "HR")
    } else {
      b <- s[, c("efficacyBounds", "efficacyP", "information")]
      
      # format number of digits after decimal for each column
      j2 <- 3
      j3 <- 1
      j4 <- 2
      
      b[j2] <- lapply(b[j2], formatC, format = "f", digits = 2)
      b[j3] <- lapply(b[j3], formatC, format = "f", digits = 3)
      b[j4] <- lapply(b[j4], formatC, format = "f", digits = 4)      
      
      df = t(b)
      
      rownames(df) = c("Efficacy boundary (Z-scale)",
                       "Efficacy boundary (p-scale)",
                       "Information")
    }
    
    colnames(df) <- NA
  }
  
  print(df1, ..., na.print = "" , quote = FALSE )
  print(df, ..., na.print = "" , quote = FALSE )
  invisible(x)
}


#' @title Print simulation results
#' @description Prints the summary statistics from simulation.
#'
#' @param x The lrsim object to print.
#' @param ... Ensures that all arguments starting from "..." are named.
#'
#' @return A tabular printout of the summary statistics from simulation runs.
#'
#' @keywords internal
#'
#' @export
print.lrsim <- function(x, ...) {
  a = x$overview;
  k = a$kMax;
  
  if (k>1) {
    str1 = paste0("Group-sequential trial with ", k, " stages")
  } else {
    str1 = "Fixed design"
  }
  
  if (a$rho1 != 0 || a$rho2 != 0) {
    str1 = paste0(str1, ", FH(", a$rho1, ", ", a$rho2, ")")
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
  row.names(df1) = c(str1, str2, str3, str4, str5, str6, str7, "")
  
  
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
