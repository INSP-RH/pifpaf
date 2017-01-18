#' @title PIF-APP: An app to calculate the potential impact fraction in a user friendly manner. 
#' 
#' @description This function loads the app to calculate the Potential Impact Fraction \code{\link{pif}} in a friendly manner.
#'  One can calculate point estimates for \code{\link{pif}} and \code{\link{paf}}, confidence intervals for PAF.  
#'  Plots are also outputs such as sensitivity analysis, a heatmap for the counterfactual, a graphical display of the counterfactual,
#'  and a graph of PAF and PIF for different theta values. Most common relative risk functions and counterfactual functions are offered to the user. 
#'  The user can load data in .csv files.
#' 
#' @import shiny
#' @import gridExtra 
#' @import tools 
#' @import ggplot2
#' 
#' @author Rodrigo Zepeda Tello \email{rodrigo.zepeda@insp.mx}
#' @author Dalia Camacho García Formentí \email{daliaf172@gmail.com}
#' 
#' @export

pif.app <- function()
{
  # appDir <- system.file("shiny-examples", "PIFapp", package = "pif")
  # if (appDir == "") {
  #   stop("Could not find directory. Try re-installing `pif`.", call. = FALSE)
  # }
  # 
  # shiny::runApp(appDir, display.mode = "normal")
  stop("Shiny App comming soon")
}