#' @title Approximate Confidence Intervals for the Potential Impact Fraction
#' 
#' @description Function that calculates confidence intervals to the potential impact fraction
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
#' @author Rodrigo Zepeda Tello \email{rodrigo.zepeda@insp.mx}
#' @author Dalia Camacho García Formentí 
#' 
#' @examples 
#' 
#' #Example with risk given by HR
#' set.seed(18427)
#' X <- rnorm(100)
#' thetahat <- 0.12
#' thetavar <- 0.1
#' paf.confidence.linear(X, thetahat, thetavar, function(X, theta){exp(theta*X)})
#' 
#' #With larger sample the variance reduces
#' set.seed(18427)
#' X <- rnorm(10000)
#' thetahat <- 0.12
#' thetavar <- 0.1
#' paf.confidence.linear(X, thetahat, thetavar, function(X, theta){exp(theta*X)})
#' 
#' @export


paf.variance.linear <- function(X, thetahat, thetasd, rr, weights =  rep(1/nrow(as.matrix(X)),nrow(as.matrix(X))), nsim = 1000){
  
  #Set X as matrix
  .X  <- as.matrix(X)
  
  #Set a minimum for nsim
  .nsim <- max(nsim,10)
  
  #Calculate the conditional expected value as a function of theta
  .pafexp <- function(.theta){
    return(1 - 1/weighted.mean(rr(.X,.theta), weights))
  }
  
  #Calculate the conditional variance as a function of theta
  s  <- sum(weights)
  s2 <- sum(weights^2)
  .pafvar <- function(.theta){
    .RO  <- weighted.mean(rr(.X,.theta), weights)
    .var <- ( s / (s^2 - s2) ) * weighted.mean((rr(.X,.theta) - .RO)^2, weights)
    return(.var)
  }
  
  #Get expected value and variance of that
  .meanvec   <- rep(NA, .nsim)
  .varvec    <- rep(NA, .nsim)
  .thetasim  <- rnorm(.nsim, thetahat, thetasd)
  for (i in 1:.nsim){
    .meanvec[i]  <- .pafexp(.thetasim[i])
    .varvec[i]   <- .pafvar(.thetasim[i])
  }
  
  #Get variance of that
  .varpaf <- var(.meanvec) + mean(.varvec)
  
  #Return variance
  return(.varpaf)
  
}


