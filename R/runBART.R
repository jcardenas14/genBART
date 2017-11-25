#' Run BART Shiny App
#'
#' This function runs the BART shiny app.
#' @examples
#' if (interactive()) genBart::runBart()
#' @import ggplot2
#' @export
runBart <- function() {
  #appDir <- system.file("shiny-app", "BART", package = "genBart")
  #if (appDir == "") {
    #stop("Could not find BART directory. Try re-installing `genBart`.", 
         #call. = FALSE)
  #}
  #shiny::runApp(appDir, display.mode = "normal")
  shiny::runGitHub("BART", "jcardenas14")
}