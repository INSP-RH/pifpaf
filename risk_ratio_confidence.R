#' @title Confidence intervals for the Risk Ratio Integral
#' 
#' @description Function that calculates confidence interval for the integral ∫RR(X;theta)f(x)dx
#' 
#' @param X         Random sample (can be vector or matrix) which includes exposure and covariates.
#' 
#' @param thetahat  Estimative of \code{theta} for the Relative Risk function
#' 
#' @param thetasd   Estimator of standard deviation of thetahat (usually standard error)
#' 
#' @param rr        Function for relative risk
#' 
#' 
#' **Optional**
#' 
#' @param weights   Survey \code{weights} for the random sample \code{X}
#' 
#' @param nsim      Number of simulations
#' 
#' @param confidence Confidence level \% (default: 95)
#' 
#' @param force.min Boolean indicating whether to force the RR to have a 
#'                  minimum value of 1 instead of 0 (not recommended).
#' 
#' @author Rodrigo Zepeda Tello \email{rodrigo.zepeda@insp.mx}
#' @author Dalia Camacho García Formentí 
#' 
#' @examples 
#' 
#' #Example with risk given by HR
#' set.seed(18427)
#' X <- rnorm(100)
#' thetahat <- 0.4
#' thetavar <- 0.1
#' risk.ratio.confidence(X, thetahat, thetavar, function(X, theta){exp(theta*X)})
#' 
#' #With larger sample the variance reduces
#' set.seed(18427)
#' X <- rnorm(10000)
#' thetahat <- 0.12
#' thetavar <- 0.1
#' risk.ratio.confidence(X, thetahat, thetavar, function(X, theta){exp(theta*X)})
#' 
#' #We can force RR's CI to be >= 1
#' risk.ratio.confidence(X, thetahat, thetavar, function(X, theta){exp(theta*X)}, force.min = TRUE)
#' 
#' #Example with theta and X multivariate
#' set.seed(18427)
#' X1 <- rnorm(1000)
#' X2 <- rnorm(1000)
#' X  <- as.matrix(cbind(X1,X2))
#' thetahat <- c(0.12, 0.03)
#' thetasd <- matrix(c(0.1, 0, 0, 0.4), byrow = TRUE, nrow = 2)
#' rr <- function(X, theta){exp(theta[1]*X[,1] + theta[2]*X[,2])}
#' risk.ratio.confidence(X, thetahat, thetasd, rr) 
#' 
#' @import MASS
#' @export

risk.ratio.confidence <- function(X, thetahat, thetasd, rr, weights =  rep(1/nrow(as.matrix(X)),nrow(as.matrix(X))),
                                  nsim = 1000, confidence = 95, force.min = FALSE){
  
  #Get confidence
  .alpha <- max(0, 1 - confidence/100)
  .Z     <- qnorm(1-.alpha/2)
  
  #Get number of simulations
  .nsim  <- max(10, ceiling(nsim))
  
  #To matrix
  .X     <- as.matrix(X)
  
  #Calculate the conditional expected value as a function of theta
  .Risk  <- function(.theta){
    .RO  <- weighted.mean(rr(.X,.theta), weights)
    return( .RO )
  }
  
  #Calculate the conditional variance as a function of theta
  .Variance <- function(.theta){
    s   <- sum(weights)
    s2   <- sum(weights^2)
    .var <- ( s / (s^2 - s2) ) * weighted.mean((rr(.X,.theta) - .Risk(.theta))^2, weights)
    return(.var)
  }
  
  
  #Get expected value and variance of that
  .meanvec   <- rep(NA, .nsim)
  .varvec    <- rep(NA, .nsim)
  .thetasim  <- mvrnorm(.nsim, thetahat, thetasd)
  for (i in 1:.nsim){
    .meanvec[i]  <- .Risk(.thetasim[i,])
    .varvec[i]   <- .Variance(.thetasim[i,])
  }
  
  #Get variance of that
  .inversevarpaf <- var(.meanvec) + mean(.varvec)
  
  #Create the confidence intervals
  .squareroot <- .Z*sqrt(.inversevarpaf)
  .ciup       <- .Risk(thetahat) + .squareroot
  .cilow      <- (.Risk(thetahat)^2)/.ciup
  
  #If minimum is forced to 1 correct CI
  if (force.min){
    .cilow      <- ((.Risk(thetahat) - 1)^2)/(.ciup-1) + 1
  }
  
  
  #Compute the Risk Ratio intervals
  .cirisk        <- c("Lower" = .cilow, "Point_Estimate" =  .Risk(thetahat) , "Upper" = .ciup )
  
  #Return variance
  return(.cirisk)
  
}
