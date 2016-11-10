#' @title Approximate Variance for the Population Attributable Fraction using the approximate method
#' 
#' @description Function that calculates approximate variance to the population attributable fraction.
#' 
#' @param Xmean     Mean value of exposure levels from a previous study.
#' 
#' @param Xvar      Variance of exposure levels from a previous study.
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
#' @param nsim      Number of simulations
#' 
#' @param check_thetas Checks that theta parameters are correctly inputed
#' 
#' @author Rodrigo Zepeda Tello \email{rodrigo.zepeda@insp.mx}
#' @author Dalia Camacho García Formentí 
#' 
#' @examples 
#' 
#' #Example 1
#' set.seed(46987)
#' rr      <- function(X,theta){exp(X*theta)}
#' Xmean   <- 3
#' Xvar    <- 1
#' theta   <- 1.2
#' thetasd <- 0.15
#' .Xmean  <- as.matrix(Xmean)
#' .Xvar   <- as.matrix(Xvar)
#' paf.variance.approximate(Xmean,Xvar,theta,thetasd,rr)
#' 
#' #Example 2: Compare paf.variance.approximate with paf.variance.linear
#' X1       <- rnorm(1000,3,.5)
#' X2       <- rnorm(1000,4,1)
#' X        <- as.matrix(cbind(X1,X2))
#' Xmean    <- colMeans(X)
#' Xvar     <- cov(X)
#' .Xmean   <- matrix(Xmean, ncol = length(Xmean))
#' .Xvar    <- matrix(Xvar, ncol = sqrt(length(Xvar)))
#' theta    <- c(0.12, 0.17)
#' thetasd  <- matrix(c(0.001, 0.00001, 0.00001, 0.004), byrow = TRUE, nrow = 2)
#' rr       <- function(X, theta){exp(theta[1]*X[,1] + theta[2]*X[,2])}
#' paf.variance.approximate(Xmean,Xvar,theta,thetasd,rr)
#' paf.variance.linear(X, theta, thetasd, rr)
#' 
#' 
#' @import MASS stats numDeriv
#' @export 


paf.variance.approximate <- function(Xmean, Xvar, thetahat, thetasd, rr, 
                                     check_thetas = TRUE, nsim = 1000){
  
  #Set X as matrix
  .Xmean  <- matrix(Xmean, ncol = length(Xmean))
  .Xvar   <- matrix(Xvar, ncol = sqrt(length(Xvar)))
  #Set a minimum for nsim
  .nsim        <- max(nsim,10)
  
  if(is.positive.semi.definite(.Xvar) == FALSE){
    stop("Variance matrix must be positive definite.")
  }
  
  #Function for checking that thetas are correctly inputed
  if(check_thetas){ check.thetas(thetasd, thetahat, NA, NA, "approximate") }
  
  #Check exposure values are greater than zero
  check.exposure(Xmean)
  
  #Check that rr is 1 when X = 0
  check.rr(.Xmean, thetahat, rr)
  
  
  
  #Rewrite functions as functions of X only
  rr.fun.x <- function(X){
    rr(X,thetahat)
  }
  
  #Get the expected paf
  .pafexp <- function(theta){
    pif.approximate(Xmean, Xvar, thetahat = theta, rr)
  }
  
  #Get the variance of PAF
  .pafvar <- function(theta){
    rr.fun.x <- function(X){
      rr(X,theta)
    }
    
    dR0 <- as.matrix(grad(rr.fun.x, .Xmean))
    R0  <- rr(.Xmean, theta)
    vr  <- t(dR0/(R0^2))%*%.Xvar%*%(dR0/(R0^2))
    return(vr)
  }
  
  #Get expected value and variance of that
  .meanvec   <- rep(NA, .nsim)
  .varvec    <- rep(NA, .nsim)
  .thetasim  <- mvrnorm(.nsim, thetahat, thetasd)
  for (i in 1:.nsim){
    .meanvec[i]  <- .pafexp(.thetasim[i,])
    .varvec[i]   <- .pafvar(.thetasim[i,])
  }
  
  #Get variance of that
  .varpaf <- var(.meanvec) + mean(.varvec)
  
  #Return variance
  return(.varpaf)
  
}


