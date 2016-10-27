#' @title Confidence Intervals for the Population Attributable Fraction when only mean and variance of exposure values is available using loglinear method
#' 
#' @description Confidence intervals for the Population Attributable Fraction for the approximate method where only mean and variance from a previous study is available.For relative risk inyective functions, the PAF is inyective, and intervals can be calculated for log(PAF), and then transformed to PAF CI.
#' 
#' @param Xmean     Mean value of exposure levels from a previous study.
#' 
#' @param Xvar      Variance of exposure levels from a previous study.
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
#' #Example 1
#' set.seed(46987)
#' rr      <- function(X,theta){exp(X*theta)}
#' Xmean   <- 3
#' Xvar    <- 1
#' theta   <- 0.4
#' thetasd <- 0.001
#' paf.confidence.approximate.loglinear(Xmean, Xvar, theta, thetasd, rr)
#'
#'#Example 2: Compare paf.variance.approximate with paf.variance.loglinear
#'X1       <- rnorm(1000,3,.5)
#'X2       <- rnorm(1000,4,1)
#'X        <- as.matrix(cbind(X1,X2))
#'Xmean    <- colMeans(X)
#'Xvar     <- cov(X)
#'theta    <- c(0.12, 0.17)
#'thetasd  <- matrix(c(0.001, 0.00001, 0.00001, 0.004), byrow = TRUE, nrow = 2)
#'rr       <- function(X, theta){exp(theta[1]*X[,1] + theta[2]*X[,2])}
#'
#'paf.confidence.approximate.loglinear(Xmean, Xvar, theta, thetasd, rr)
#'paf.confidence.loglinear(X, theta, thetasd, rr)
#' 
#' @import MASS numDeriv
#' @export

paf.confidence.approximate.loglinear <- function(Xmean, Xvar, thetahat, thetavar, rr,
                                     nsim = 1000, confidence = 95, check_thetas = TRUE){
  
  #Get confidence
  .alpha <- max(0, 1 - confidence/100)
  .Z     <- qnorm(1-.alpha/2)
  
  #Get number of simulations
  .nsim  <- max(10, ceiling(nsim))
  
  #Check 
  .thetavar <- as.matrix(thetavar)
  if(check_thetas){ check.thetas(.thetavar, thetahat, NA, NA, "log") }
  
  #Set X as matrix
  .Xmean  <- matrix(Xmean, ncol = length(Xmean))
  .Xvar   <- matrix(Xvar, ncol = sqrt(length(Xvar)))
  
  #Calculate the conditional expected value as a function of theta
  .logpafexp <- function(.theta){
    .paf  <- pif.approximate(.Xmean, .Xvar,.theta, rr)
    .R0   <- 1/(1-.paf)
    return( -log(.R0) )
  }
  
  #Inverse
  .paf      <- pif.approximate(.Xmean, .Xvar, thetahat, rr)
  .inverse  <- 1 - .paf
  
  #Calculate the conditional variance as a function of theta
 
  .logpafvar <- function(.theta){
    rr.fun.x <- function(X){
      rr(X,.theta)
    }
    
    dR0   <- as.matrix(grad(rr.fun.x, .Xmean))
    R0    <- rr(.Xmean, theta)
    .var  <- t(1/R0*dR0)%*%.Xvar%*%(1/R0*dR0)
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

  
  #Compute the PAF intervals
  .cipaf         <- 1-c("Lower" = .inverse*exp(.zqrt), "Point_Estimate" =  .inverse, "Upper" = .inverse*exp(-.zqrt) )
  
  #Return variance
  return(.cipaf)
  
}