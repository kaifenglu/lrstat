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
  s = x$overview
  k = length(s$numberOfEvents)
  if (k>1) {
    df = t(data.frame(s$cumulativeRejection,
                      s$cumulativeFutility,
                      s$numberOfEvents,
                      s$numberOfDropouts,
                      s$numberOfSubjects,
                      s$analysisTime,
                      s$overallReject,
                      s$expectedNumberOfEvents,
                      s$expectedNumberOfDropouts,
                      s$expectedNumberOfSubjects,
                      s$expectedStudyDuration))
    df[c(7,8,9,10,11), -1] <- NA # only show overall
    colnames(df) <- paste("stage", seq_len(ncol(df)), sep=" ")
  } else {
    df = t(data.frame(s$overallReject,
                      s$expectedNumberOfEvents,
                      s$expectedNumberOfDropouts,
                      s$expectedNumberOfSubjects,
                      s$expectedStudyDuration))
    colnames(df) <- NA
  }
  rownames(df) <- sub("^[[:alpha:]][[:alnum:]]*.", "", rownames(df))
  print( round(df,3), ..., na.print = "" , quote = FALSE )
  invisible(x)
}
