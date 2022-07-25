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
  s = x$byStageResults;
  t = x$overallResults;
  k = length(s$informationRates)
  if (k>1) {
    if (t$estimateHazardRatio) {
      df = t(data.frame(s$informationRates,
                        s$efficacyBounds,
                        s$futilityBounds,
                        s$cumulativeRejection,
                        s$cumulativeFutility,
                        s$cumulativeAlphaSpent,
                        s$numberOfEvents,
                        s$numberOfDropouts,
                        s$numberOfSubjects,
                        s$analysisTime,
                        s$efficacyHR,
                        s$futilityHR,
                        s$efficacyP,
                        s$futilityP,
                        s$information,
                        s$HR,
                        t$overallReject,
                        t$alpha,
                        t$numberOfEvents,
                        t$expectedNumberOfEvents,
                        t$numberOfDropouts,
                        t$expectedNumberOfDropouts,
                        t$numberOfSubjects,
                        t$expectedNumberOfSubjects,
                        t$studyDuration,
                        t$expectedStudyDuration,
                        t$accrualDuration,
                        t$followupTime,
                        t$fixedFollowup,
                        t$rho1,
                        t$rho2))
      df[seq(17,31), -1] <- NA # only show overall      
    } else {
      df = t(data.frame(s$informationRates,
                        s$efficacyBounds,
                        s$futilityBounds,
                        s$cumulativeRejection,
                        s$cumulativeFutility,
                        s$cumulativeAlphaSpent,
                        s$numberOfEvents,
                        s$numberOfDropouts,
                        s$numberOfSubjects,
                        s$analysisTime,
                        s$efficacyP,
                        s$futilityP,
                        s$information,
                        t$overallReject,
                        t$alpha,
                        t$numberOfEvents,
                        t$expectedNumberOfEvents,
                        t$numberOfDropouts,
                        t$expectedNumberOfDropouts,
                        t$numberOfSubjects,
                        t$expectedNumberOfSubjects,
                        t$studyDuration,
                        t$expectedStudyDuration,
                        t$accrualDuration,
                        t$followupTime,
                        t$fixedFollowup,
                        t$rho1,
                        t$rho2))
      df[seq(14,28), -1] <- NA # only show overall      
    }
    
    colnames(df) <- paste("stage", seq_len(ncol(df)), sep=" ")
  } else {
    if (t$estimateHazardRatio) {
      df = t(data.frame(t$overallReject,
                        t$alpha,
                        t$numberOfEvents,
                        t$numberOfDropouts,
                        t$numberOfSubjects,
                        t$studyDuration,
                        t$accrualDuration,
                        t$followupTime,
                        t$fixedFollowup,
                        t$rho1,
                        t$rho2,
                        s$efficacyBounds,
                        s$efficacyHR,
                        s$efficacyP,
                        s$information,
                        s$HR))
    } else {
      df = t(data.frame(t$overallReject,
                        t$alpha,
                        t$numberOfEvents,
                        t$numberOfDropouts,
                        t$numberOfSubjects,
                        t$studyDuration,
                        t$accrualDuration,
                        t$followupTime,
                        t$fixedFollowup,
                        t$rho1,
                        t$rho2,
                        s$efficacyBounds,
                        s$efficacyP,
                        s$information))
    }
    
    colnames(df) <- NA
  }
  rownames(df) <- sub("^[[:alpha:]][[:alnum:]]*.", "", rownames(df))
  print( round(df,3), ..., na.print = "" , quote = FALSE )
  invisible(x)
}
