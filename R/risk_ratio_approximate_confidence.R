#' @title Confidence intervals for the Risk Ratio Integral
#' 
#' @description Function that calculates confidence interval for the integral ∫RR(X;theta)f(x)dx
#' 
#' @param Xmean     Mean value of exposure levels from a previous study.
#' 
#' @param Xvar      Variance of exposure levels from a previous study.
#' 
#' @param thetahat  Estimative of \code{theta} for the Relative Risk function
#' 
#' @param thetavar   Estimator of variance of thetahat
#' 
#' @param rr        Function for relative risk
#' 
#' 
#' **Optional**
#'
#' 
#' @param nsim      Number of simulations
#' 
#' @param confidence Confidence level \% (default: 95)
#' 
#' @param check_thetas Checks that theta parameters are correctly inputed
#' 
#' @param force.min Boolean indicating whether to force the RR to have a 
#'                  minimum value of 1 instead of 0 (not recommended).
#' 
#' @author Rodrigo Zepeda Tello \email{rodrigo.zepeda@insp.mx}
#' @author Dalia Camacho García Formentí 
#' 
#' @examples 
#' set.seed(46987)
#' X       <- rnorm(100,3.2,1)
#' Xmean   <- 3.2
#' Xvar    <- 1
#' theta   <- 0.4
#' thetavar <- 0.001
#' rr      <- function(X,theta){exp(X*theta)}
#' risk.ratio.approximate.confidence(Xmean, Xvar, thetahat = theta, thetavar = thetavar, rr = rr)
#' @import MASS numDeriv
#' @export

risk.ratio.approximate.confidence <- function(Xmean, Xvar, thetahat, thetavar, rr,
                                  nsim = 1000, confidence = 95, check_thetas = TRUE,
                                  force.min = FALSE){
  
  #Get confidence
  check.confidence(confidence)
  .alpha <- max(0, 1 - confidence/100)
  .Z     <- qnorm(1-.alpha/2)
  
  #Check 
  .thetavar <- as.matrix(thetavar)
  if(check_thetas){ check.thetas(.thetavar, thetahat, NA, NA, "risk.ratio") }
  
  #Get number of simulations
  .nsim  <- max(10, ceiling(nsim))
  
  #To matrix
  .Xmean  <- matrix(Xmean, ncol = length(Xmean))
  .Xvar   <- matrix(Xvar, ncol = sqrt(length(Xvar)))
  
  #Calculate the conditional expected value as a function of theta
  .Risk  <- function(.theta){
    .paf  <- pif.approximate(Xmean =  .Xmean, Xvar = .Xvar, thetahat =.theta, rr=rr)
    return( 1/(1-.paf))
  }
  
  #Calculate the conditional variance as a function of theta
  .Variance <- function(.theta){
    rr.fun.x <- function(X){
      rr(X,.theta)
    }
    
    dR0   <- as.matrix(grad(rr.fun.x, .Xmean))
    .var  <- t(dR0)%*%.Xvar%*%(dR0)
    return(.var)
  }
  
  
  #Get expected value and variance of that
  .meanvec   <- rep(NA, .nsim)
  .varvec    <- rep(NA, .nsim)
  .thetasim  <- mvrnorm(.nsim, thetahat, thetavar, empirical = TRUE)
  for (i in 1:.nsim){
    .meanvec[i]  <- .Risk(.thetasim[i,])
    .varvec[i]   <- .Variance(.thetasim[i,])
  }
  
  #Get variance of that
  .inversevarpaf <- var(.meanvec) + mean(.varvec)
  
  #Create the confidence intervals
  .squareroot <- .Z*sqrt(.inversevarpaf)
  .ciup       <- .Risk(thetahat) + .squareroot
  .cilow      <- (.Risk(thetahat)^2)/.ciup #?
  
  #If minimum is forced to 1 correct CI
  if (force.min){
    .cilow      <- ((.Risk(thetahat) - 1)^2)/(.ciup-1) + 1
  }
  
  
  #Compute the Risk Ratio intervals
  .cirisk        <- c("Lower" = .cilow, "Point_Estimate" =  .Risk(thetahat) , "Upper" = .ciup )
  
  #Return variance
  return(.cirisk)
  
}
