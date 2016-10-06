#Load libraries
library(shiny)
library(gridExtra)
library(pif)
library(tools)
library(matrixcalc)
library(ggplot2)
library(latex2exp)
library(MASS)

#Text Input
textInputRow<-function (inputId, label, value = "")
{
  div(style="display:inline-block",
      tags$label(label, `for` = inputId),
      tags$input(id = inputId, type = "text", value = value,class="input-small"))
}


#Warning messages for Numeric Results 
Warnings1 <- function(Error){
  switch (Error,
          return("Please load file with exposure levels in .csv format"),
          return("Warning: Incorrect file extension. <br/> File with exposure levels must be .csv format."),
          return("Warning: $$\\theta$$  values should be greater or equal to zero."),
          return("Please load Counterfactual function in an R file with the function specified as:  <br/> Counterfactual <- function(X){function goes here}.<br/> Example: Counterfactual <- function(X){2*X}."),
          return("Warning: Incorrect file extension. Counterfactual function must be written in an R file"),
          return("Warning: There is no function named Counterfactual. The counterfactual function loaded must be named Counterfactual and it must be defined as: <br/> Counterfactual <- function(X){function goes here}. <br/> Example: Counterfactual <- function(X){2*X}."),
          return("Warning: Counterfactual must be a function and must be defined as: <br/> Counterfactual <- function(X){function goes here}. <br/> Example: Counterfactual <- function(X){2*X}. "),
          return("Warning: Counterfactual function could not be evaluated."),
          return("Please Load weights in .csv format"),
          return("Warning: Incorrect format for weights, please load in .csv format"),
          return("Please load Relative Risk function in an R file with the function specified as: <br/> RRfunction <- function(X, theta){function goes here}. <br/> Example: RRfunction <- function(X, theta){exp(theta*X)}."),
          return("Warning: Incorrect file extension. <br/> Relative Risk function must be written in an R file."),
          return("Warning: There is no function named RRfunction. <br/> The Relative Risk function loaded must be named RRfunction and it must be defined as: <br/> RRfunction <- function(X, theta){function goes here}. <br/> Example: RRfunction <- function(X, theta){exp(theta*X)}."),
          return("Warning: RRfunction must be a function and must be defined as: <br/> RRfunction <- function(X, theta){function goes here}. <br/> Example: RRfunction <- function(X, theta){exp(theta*X)}."),
          return("Please load file with the values for $$ \\theta$$"),
          return("Warning: Incorrect file extension. <br/> $$\\theta$$ values must be in a .csv format."),
          return("Warning: Relative Risk Function could not be evaluated."),
          return("Warning: Relative Risk function at minimum exposure level is not equal to one."),
          return("Warning: Minimum exposure level must be greater or equal to zero."),
          return("Warning: The Relative Risk function at minimum exposure level must be equal to one."),
          return("Please load file with square root of variance and covariance matrix"),
          return("Warning: Incorrect file extension. <br/> File with square root of variance and covariance matrix must be .csv format"),
          return("Warning: Square root of variance and covariance matrix must be a matrix."),
          return("Warning: Square root of variance and covariance matrix must be positive definite.")
          
          
  )
}

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