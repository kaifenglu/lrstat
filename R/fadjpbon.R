#' @title Adjusted p-values for graphical approaches
#' @description Obtains the adjusted p-values for graphical approaches 
#' using weighted Bonferroni tests for fixed design.
#'
#' @param w The vector of initial weights for elementary hypotheses.
#' @param G The initial transition matrix.
#' @param p Raw p-values for elementary hypotheses.
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

