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
#'@examples 
#'#Example 1
#' Xmean     <- 2
#' thetahat  <- 0.12
#' Xvar      <- 0.2
#' rr        <- function(X,theta){exp(X*theta)}
#' cft       <- function(X){0.5*X}
#' pif.approximate(Xmean, Xvar, thetahat, rr, cft)
#' pif.approximate(Xmean, Xvar, thetahat, rr)
#' 
#' #Example 2: Multivariate
#' #Multivariate example
#' 
#' X1        <- 2
#' X2        <- 1.1
#' Xmean     <- as.matrix(cbind(X1,X2))
#' Xvar      <- matrix(c(1,.4,.4,1),ncol = 2, byrow = TRUE)
#' Cft       <- function(X){.25*X}
#' thetahat  <- c(0.12, 0.03)
#' rr <- function(X, theta){exp(theta[1]*X[,1] + theta[2]*X[,2])}
#' pif.approximate(Xmean, Xvar, thetahat, rr) 
#' pif.approximate(Xmean, Xvar, thetahat, rr, Cft)
#' 
#' @import matrixStats numDeriv
#'  
#' @export


pif.approximate <- function(Xmean, Xvar, thetahat, rr, 
                            cft = function(Xmean){matrix(0,ncol = ncol(as.matrix(Xmean)), nrow = nrow(as.matrix(Xmean)))}){
  
  .Xmean  <- matrix(Xmean, ncol = length(Xmean))
  .Xvar   <- matrix(Xvar, ncol = sqrt(length(Xvar)))
  
  if(is.positive.semi.definite(.Xvar) == FALSE){
    stop("Covariance matrix must be positive semi-definite.")
  }
  
  #Check exposure values are greater than zero
  check.exposure(.Xmean)
  
  #Check that rr is 1 when X = 0
  check.rr(.Xmean, thetahat, rr)
  
  #Check counterfactual
  check.cft(cft, .Xmean)
  
  #Rewrite functions as functions of X only
  rr.cft.fun <- function(X){
    Xcft.value  <- cft(X)
    rr(Xcft.value,thetahat)
  } 
  
  rr.fun.x <- function(X){
    rr(X,thetahat)
  }
  
  
  #Estimate weighted sums
  .mucft   <- rr.cft.fun(.Xmean) + 0.5*EntryMult(hessian(rr.cft.fun,.Xmean), .Xvar)
  .mux     <- rr.fun.x(.Xmean) + 0.5*EntryMult(hessian(rr.fun.x,.Xmean), .Xvar)

  #Check that integrals make sense
  check.integrals(.mux, .mucft)
  
  #Calculate PIF
  .pif   <- as.numeric(1 - .mucft/.mux)
  
  return(.pif)
  
}
