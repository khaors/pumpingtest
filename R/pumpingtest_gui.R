#' @title
#' pumpingtest_gui
#' @description
#' pumpingtest GUI
#' @import shiny
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @export
pumpingtest_gui <- function() {
  appDir <- system.file("Shiny", "pumpingtest", package = "pumpingtest")
  if (appDir == "") {
    stop("Could not find example directory. Try re-installing `pumpingtest`.", call. = FALSE)
  }
  shiny::runApp(appDir, display.mode = "normal")
}
