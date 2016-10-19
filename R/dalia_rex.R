#' @title DALIA-REX: Dalia's app
#' 
#' @description Does something fancy
#' 
#' @import shiny gridExtra tools latex2exp
#' @export

dalia.rex <- function()
{
  
  appDir <- system.file("shiny-examples", "PIFapp", package = "pif")
  if (appDir == "") {
    stop("Could not find example directory. Try re-installing `pif`.", call. = FALSE)
  }
  
  shiny::runApp(appDir, display.mode = "normal")
}

