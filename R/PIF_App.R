#' @title PIF-APP: An app to calculate the potential impact fraction in a user friendly manner. 
#' 
#' @description This function loads the app to calculate the potential impact fraction in a friendly manner.
#'  One can calculate point estimates for PIF and PAF, confidence intervals for PAF.  
#'  Graphs are also outputs such as sensitivity analysis, a 3d graph for the counterfactual, and a graph of PAF and PIF for different theta values. 
#'  Most common relative risk functions and counterfactual functions are offered to the user. The user can load data in .csv files.
#' 
#' @import shiny
#' @import gridExtra 
#' @import tools 
#' @import latex2exp 
#' @import ggplot2
#' @export

pif.app <- function()
{
  appDir <- system.file("shiny-examples", "PIFapp", package = "pif")
  if (appDir == "") {
    stop("Could not find example directory. Try re-installing `pif`.", call. = FALSE)
  }
  
  shiny::runApp(appDir, display.mode = "normal")
}