#' Run BART Shiny App
#'
#' This function runs the BART shiny app 
#' @examples
#' ## Only run this example in interactive R sessions
#' if (interactive()) {
#'   runBart()
#' }
#' @import ggplot2
#' @export
runBart <- function() {
  shiny::runGitHub("BART", "jcardenas14")
}