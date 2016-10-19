#' @title Warnings plot
#' 
#' @description Please describe
#' 
#' @export

#Warning messages for Plots
Warnings2 <- function(Error){
  switch (Error,
          return(ggplot()+ggtitle("Please load file with exposure levels in .csv format.")),
          return(ggplot()+ggtitle("Warning: Incorrect file extension. File with exposure levels must be .csv format.")),
          return(ggplot()+ggtitle(TeX('Warning: Minimum value for $\\theta$ must be greater or equal to zero.'))),
          return(ggplot()+ggtitle("Please load Counterfactual function in an R file \n with the function specified as: \n Counterfactual <- function(X){function goes here}. \n Example: Counterfactual <- function(X){2*X} .")),
          return(ggplot()+ggtitle("Warning: Incorrect file extension. \n  Counterfactual function must be written in an R file.")),
          return(ggplot()+ggtitle("Warning: There is no function named Counterfactual. \n The counterfactual function loaded must be named Counterfactual \n and  it must be defined as \n Counterfactual <- function(X){function goes here} . \n  Example: Counterfactual <- function(X){2*X} .")),
          return(ggplot()+ggtitle("Warning: Counterfactual must be a function \n and must be defined as: \n  Counterfactual <- function(X){function goes here} .\n Example: Counterfactual <- function(X){2*X}.")),
          return(ggplot()+ggtitle("Warning: Counterfactual function could not be evaluated.")),
          return(ggplot()+ggtitle("Load weights in .csv format.")),
          return(ggplot()+ggtitle("Incorrect format for weights, please load in .csv format.")),
          return(ggplot()+ggtitle("Please load Relative Risk function in an R file \n with the function specified as: \n RRfunction <- function(X, theta){function goes here}. \n  Example: RRfunction <- function(X, theta){exp(theta*X)} .")),
          return(ggplot()+ggtitle("Warning: Incorrect file extension. \n Relative Risk function must be written in an R file.")),
          return(ggplot()+ggtitle("Warning: There is no function named RRfunction.\n The Relative Risk function loaded must be \n named RRfunction and it must be defined as: \n  RRfunction <- function(X, theta){function goes here}.\n Example: RRfunction <- function(X, theta){exp(theta*X)}.  ")),
          return(ggplot()+ggtitle("Warning: RRfunction must be a function and must be defined as RRfunction <- function(X, theta){function goes here}. Example: RRfunction <- function(X, theta){exp(theta*X)}. ")),
          return(ggplot()+ggtitle(TeX('Please load file with the values for $\\theta$.'))),
          return(ggplot()+ggtitle(TeX('Warning: Incorrect file extension.  $\\theta$ values must be in a .csv format.'))),
          return(ggplot()+ggtitle("Relative Risk Function could not be evaluated.")),
          return(ggplot()+ggtitle("Warning: Relative Risk function at minimum \n exposure level is not equal to one.")),
          return(ggplot()+ggtitle("Warning: Minimum exposure level must be greater or equal to zero.")),
          return(ggplot()+ggtitle(TeX('Minimum $\\theta$ value must be smaller  than the maximum $\\theta$ value.'))),
          return(ggplot()+ggtitle(TeX('Graph can only be plotted if  the Relative Risk function depends on a one dimensional $\\theta$.'))),
          return(ggplot()+ggtitle(TeX('Warning: Relative Risk function at minimum exposure level is not equal to one. At minimum or maximum $\\theta$.')))
  )
}