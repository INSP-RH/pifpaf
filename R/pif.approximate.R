#' @title Point Estimate of the Potential Impact Fraction via the Approximate method.
#' 
#' @description Function that calculates the potential impact fraction via the approximate method.
#' 
#' @param Xmean      Mean value of exposure levels from a previous study.
#' 
#' @param Xvar       Variance of the exposure levels.
#' 
#' @param thetahat  Estimative of \code{theta} for the Relative Risk function.
#' 
#' @param rr        Function for relative risk.
#' 
#' 
#' **Optional**
#' 
#' @param cft       Function \code{cft(X)} for counterfactual. Leave empty for \code{PAF} where counterfactual is 0 exposure.
#' 
#' @author Rodrigo Zepeda Tello \email{rodrigo.zepeda@insp.mx}
#' @author Dalia Camacho García Formentí
#' 
#'@example 
#' Xmean     <- 2
#' thetahat  <- 0.12
#' Xvar      <- 0.2
#' rr        <- function(X,theta){exp(X*theta)}
#' cft       <- function(X){0.5*X}
#' pif.approximate(Xmean, Xvar, thetahat, rr, cft)
#' pif.approximate(Xmean, Xvar, thetahat, rr)
#' 
#' @import matrixStats
#' 
#' @import numDeriv
#'  
#' @export


pif.approximate <- function(Xmean, Xvar, thetahat, rr, 
                            cft = function(Xmean){matrix(0,ncol = ncol(as.matrix(Xmean)), nrow = nrow(as.matrix(Xmean)))}){
  
  .Xmean <- as.matrix(Xmean)
  
  #Check that rr is 1 when X = 0
  check.rr(.Xmean, thetahat, rr)
  
  #Rewrite functions as functions of X only
  rr.cft.fun <- function(X){
    Xcft.value  <- cft(X)
    rr(Xcft.value,thetahat)
  } 
  
  rr.fun.x <- function(X){
    rr(X,thetahat)
  }
  
  
  #Estimate weighted sums
  .mucft   <- rr.cft.fun(.Xmean) + 0.5*hessian(rr.cft.fun,.Xmean)*Xvar
  .mux     <- rr.fun.x(.Xmean) + 0.5*hessian(rr.fun.x,.Xmean)*Xvar

  
  #Check that integrals make sense
  check.integrals(.mux, .mucft)
  
  #Calculate PIF
  .pif   <- as.numeric( 1 - .mucft/.mux)
  
  return(.pif)
  
}
