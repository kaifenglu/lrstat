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
  s = x$byStageResults;
  t = x$overallResults;
  k = length(s$informationRates)
  if (k>1) {
      df = t(data.frame(s$informationRates,
                        s$efficacyBounds,
                        s$futilityBounds,
                        s$cumulativeRejection,
                        s$cumulativeFutility,
                        s$cumulativeAlphaSpent,
                        s$efficacyP,
                        s$futilityP,
                        t$overallReject,
                        t$alpha,
                        t$drift,
                        t$inflationFactor))
      df[seq(9,12), -1] <- NA # only show overall      
    
    
    colnames(df) <- paste("stage", seq_len(ncol(df)), sep=" ")
  } else {

      df = t(data.frame(t$overallReject,
                        t$alpha,
                        s$efficacyBounds,
                        s$efficacyP,
                        s$drift))

    colnames(df) <- NA
  }
  rownames(df) <- sub("^[[:alpha:]][[:alnum:]]*.", "", rownames(df))
  print( round(df,3), ..., na.print = "" , quote = FALSE )
  invisible(x)
}
