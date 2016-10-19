#' @title Approximate Variance for the Population Attributable Fraction
#' 
#' @description Function that calculates approximate variance to the population attributable fraction
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
#' @param nsim      Number of simulations
#' 
#' @param confidence Confidence level \% (default 95)
#' 
#' @param check_thetas Check that theta parameters are correctly inputed
#' 
#' @author Rodrigo Zepeda Tello \email{rodrigo.zepeda@insp.mx}
#' @author Dalia Camacho García Formentí 
#' 
#' @examples 
#' 
#' #Example with risk given by HR
#' set.seed(18427)
#' X        <- rnorm(100)
#' thetahat <- 0.4
#' thetavar <- 0.1
#' paf.confidence.loglinear(X, thetahat, thetavar, function(X, theta){exp(theta*X)})
#' 
#' #With larger sample the variance reduces
#' set.seed(18427)
#' X        <- rnorm(10000)
#' thetahat <- 0.12
#' thetavar <- 0.1
#' paf.confidence.loglinear(X, thetahat, thetavar, function(X, theta){exp(theta*X)})
#' 
#' #Example with theta and X multivariate
#' set.seed(18427)
#' X1       <- rnorm(1000)
#' X2       <- rnorm(1000)
#' X        <- as.matrix(cbind(X1,X2))
#' thetahat <- c(0.12, 0.03)
#' thetavar  <- matrix(c(0.1, 0, 0, 0.4), byrow = TRUE, nrow = 2)
#' rr       <- function(X, theta){exp(theta[1]*X[,1] + theta[2]*X[,2])}
#' paf.confidence.loglinear(X, thetahat, thetavar, rr) 
#' 
#' @import MASS
#' @export

paf.confidence.loglinear <- function(X, thetahat, thetavar, rr, weights =  rep(1/nrow(as.matrix(X)),nrow(as.matrix(X))),
                                     nsim = 1000, confidence = 95, check_thetas = TRUE){
  
  #Get confidence
  .alpha <- max(0, 1 - confidence/100)
  .Z     <- qnorm(1-.alpha/2)
  
  #Get number of simulations
  .nsim  <- max(10, ceiling(nsim))
  
  #Check 
  .thetavar <- as.matrix(thetavar)
  if(check_thetas){ check.thetas(.thetavar, thetahat, NA, NA, "log") }
  
  #To matrix
  .X     <- as.matrix(X)
  
  #Calculate the conditional expected value as a function of theta
  .logpafexp <- function(.theta){
    .RO  <- weighted.mean(rr(.X,.theta), weights)
    return( -log(.RO) )
  }
  
  #Inverse
  .inverse <- 1/weighted.mean(rr(.X,thetahat), weights)
    
  #Calculate the conditional variance as a function of theta
  s  <- sum(weights)
  s2 <- sum(weights^2)
  .logpafvar <- function(.theta){
    .RO  <- weighted.mean(rr(.X,.theta), weights)
    .var <- (1/.RO^2)*( s / (s^2 - s2) ) * weighted.mean((rr(.X,.theta) - .RO)^2, weights)
    return(.var)
  }
  
  #Get expected value and variance of that
  .logmeanvec   <- rep(NA, .nsim)
  .logvarvec    <- rep(NA, .nsim)
  .thetasim     <- mvrnorm(.nsim, thetahat, thetavar)
  for (i in 1:.nsim){
    .logmeanvec[i]  <- .logpafexp(.thetasim[i,])
    .logvarvec[i]   <- .logpafvar(.thetasim[i,])
  }
  
  #Get variance of that
  .logvarpaf <- var(.logmeanvec) + mean(.logvarvec)
  
  #Create the confidence intervals
  .zqrt       <- .Z*sqrt(.logvarpaf)
  .paf        <- pif(.X, thetahat, rr, weights = weights)
  
  #Compute the PAF intervals
  .cipaf         <- 1-c("Lower" = .inverse*exp(.zqrt), "Point_Estimate" =  .inverse, "Upper" = .inverse*exp(-.zqrt) )
  
  #Return variance
  return(.cipaf)
  
}