#' @title Run Shiny app
#' @description Runs the log-rank test power and sample size calculation 
#' Shiny app.
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @export
runShinyApp <- function() {
  shiny::shinyAppDir(system.file("shinyApp", package = "lrstat"))
}
