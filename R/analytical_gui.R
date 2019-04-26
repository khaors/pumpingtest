#' @title
#' analytical_gui
#' @description
#' analytical GUI
#' @import shiny
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @export
analytical_gui <- function() {
  appDir <- system.file("Shiny", "analytical", package = "pumpingtest")
  if (appDir == "") {
    stop("Could not find example directory. Try re-installing `pumpingtest`.", call. = FALSE)
  }
  shiny::runApp(appDir, display.mode = "normal")
}
