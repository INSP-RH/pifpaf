#' @title Approximate Variance for the Population Attributable Fraction
#' 
#' @description Function that calculates approximate variance to the population attributable fraction.
#' 
#' @param X         Random sample (can be vector or matrix) which includes exposure and covariates.
#' 
#' @param thetahat  Estimative of \code{theta} for the Relative Risk function
#' 
#' @param thetavar   Estimator of standard deviation of thetahat (usually standard error)
#' 
#' @param rr        Function for relative risk
#' 
#' 
#' **Optional**
#' 
#' @param weights   Survey \code{weights} for the random sample \code{X}
#' 
#' @param check_thetas Checks that theta parameters are correctly inputed
#' 
#' @param nsim      Number of simulations
#' 
#' @author Rodrigo Zepeda Tello \email{rodrigo.zepeda@insp.mx}
#' @author Dalia Camacho García Formentí 
#' 
#' @examples 
#' 
#' #Example with risk given by HR
#' set.seed(18427)
#' X <- rnorm(100, 3, .6)
#' thetahat <- 0.12
#' thetavar <- 0.1
#' paf.variance.linear(X, thetahat, thetavar, function(X, theta){exp(theta*X)})
#' 
#' #With larger sample the variance reduces
#' set.seed(18427)
#' X <- rnorm(10000, 4, 1)
#' thetahat <- 0.12
#' thetavar <- 0.1
#' paf.variance.linear(X, thetahat, thetavar, function(X, theta){exp(theta*X)})
#' 
#' #Example with theta and X multivariate
#' set.seed(18427)
#' X1 <- rnorm(1000,4,1)
#' X2 <- rnorm(1000,4,1)
#' X  <- as.matrix(cbind(X1,X2))
#' thetahat <- c(0.12, 0.03)
#' thetavar <- matrix(c(0.1, 0, 0, 0.4), byrow = TRUE, nrow = 2)
#' paf.variance.linear(X, thetahat, thetavar, 
#' function(X, theta){exp(theta[1]*X[,1] + theta[2]*X[,2])}) 
#' 
#' @import MASS stats
#' @export


paf.variance.linear <- function(X, thetahat, thetavar, rr, 
                                weights =  rep(1/nrow(as.matrix(X)),nrow(as.matrix(X))), 
                                check_thetas = TRUE, nsim = 1000){
  
  #Set X as matrix
  .X  <- as.matrix(X)
  
  #Function for checking that thetas are correctly inputed
  if(check_thetas){ check.thetas(thetavar, thetahat, NA, NA, "linear") }
  
  #Set a minimum for nsim
  .nsim        <- max(nsim,10)
  
  #Get the expected paf
  .pafexp <- function(theta){
    R0 <- weighted.mean(rr(.X, theta), weights)
    return(1-1/R0)
  }
  
  #Get the variance of PAF
  .pafvar <- function(theta){
    s2  <- sum(weights^2)
    s   <- sum(weights)
    R0  <- weighted.mean(rr(.X, theta), weights)
    vr  <- (1/R0^4)*(s/(s^2-s2))*weighted.mean((rr(.X, theta) - R0)^2, weights)
    return(vr)
  }
  
  #Get expected value and variance of that
  .meanvec   <- rep(NA, .nsim)
  .varvec    <- rep(NA, .nsim)
  .thetasim  <- mvrnorm(.nsim, thetahat, thetavar, empirical = TRUE)
  for (i in 1:.nsim){
    .meanvec[i]  <- .pafexp(.thetasim[i,])
    .varvec[i]   <- .pafvar(.thetasim[i,])
  }
  
  #Get variance of that
  .varpaf <- var(.meanvec) + mean(.varvec)
  
  #Return variance
  return(.varpaf)
  
}


