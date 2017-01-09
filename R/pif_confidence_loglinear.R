#' @title Confidence intervals for the Potencial Impact Fraction, using the loglinear method
#' 
#' @description Confidence intervals for the potencial impact Fraction for relative risk inyective functions, the pif is inyective, and intervals can be calculated for log(pif), and then transformed to pif CI.
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
#' @param cft       Function \code{cft(X)} for counterfactual. Leave empty for 
#'                  the Population Attributable Fraction \code{PAF} where counterfactual is 0 exposure
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
#' #Example with square root counterfactual
#' set.seed(18427)
#' X        <- rnorm(100,5,1)
#' thetahat <- 0.4
#' thetavar <- 0.1
#' cft      <- function(X){sqrt(X)}
#' pif.confidence.loglinear(X, thetahat, thetavar, function(X, theta){exp(theta*X)}, cft)
#' 
#' #Example with linear counterfactual
#' a    <- 0.5
#' cft  <- function(X){a*X}
#' pif.confidence.loglinear(X, thetahat, thetavar, function(X, theta){exp(theta*X)})
#' 
#' #Example with theta and X multivariate
#' set.seed(18427)
#' X1        <- rnorm(100,4,1)
#' X2        <- rnorm(100,4,1)
#' X         <- as.matrix(cbind(X1,X2))
#' thetahat  <- c(0.12, 0.03)
#' thetavar  <- matrix(c(0.1, 0, 0, 0.4), byrow = TRUE, nrow = 2)
#' rr        <- function(X, theta){
#'            .X <- matrix(X, ncol = 2)
#'            exp(theta[1]*.X[,1] + theta[2]*.X[,2])
#'            }
#' pif.confidence.loglinear(X, thetahat, thetavar, rr) 
#' paf.confidence.loglinear(X, thetahat, thetavar, rr)
#' 
#' @import MASS
#' @export

pif.confidence.loglinear <- function(X, thetahat, thetavar, rr, 
                                     cft = function(Varx){matrix(0,ncol = ncol(as.matrix(Varx)), nrow = nrow(as.matrix(Varx)))},
                                     weights =  rep(1/nrow(as.matrix(X)),nrow(as.matrix(X))),
                                     nsim = 100, confidence = 95, check_thetas = TRUE){
  
  #Check exposure levels
  check.exposure(X)
  #Get confidence
  check.confidence(confidence)
  .alpha <- max(0, 1 - confidence/100)
  .Z     <- qnorm(1-.alpha/2)
  
  #Check counterfactual
  check.cft(cft, X)
  
  #Get number of simulations
  .nsim  <- max(10, ceiling(nsim))
  
  #Check 
  .thetavar <- as.matrix(thetavar)
  if(check_thetas){ check.thetas(.thetavar, thetahat, NA, NA, "log") }
  
  #To matrix
  .X     <- as.matrix(X)
  .cft.X <- as.matrix(cft(X))
  n      <- dim(.X)[1]
  #Calculate the conditional expected value as a function of theta
  .logpifexp <- function(.theta){
    .RO  <- weighted.mean(rr(.X,.theta), weights)
    .RC  <- weighted.mean(rr(.cft.X,.theta), weights)
    return(log(.RC) -log(.RO) )
  }
  
  #Inverse
  .inverse <- weighted.mean(rr(.cft.X, thetahat), weights)/weighted.mean(rr(.X,thetahat), weights)
  
  #Calculate the conditional variance as a function of theta
  

  .logpifvar <- function(.theta){
    s        <- sum(weights)
    s2       <- sum(weights^2)
    .RO      <- weighted.mean(rr(.X,.theta), weights)
    .RC      <- weighted.mean(rr(.cft.X,.theta), weights)
    
    .varRO   <- (1/.RO^2)*s2*( s / (s^2 - s2) ) * weighted.mean((rr(.X,.theta) - .RO)^2, weights)
    .varRC   <- (1/.RC^2)*s2*( s / (s^2 - s2) ) * weighted.mean((rr(.cft.X,.theta) - .RC)^2, weights)
    .covRORC <- (1/(.RO*.RC))*s2*s/(s^2 - s2)   * (weighted.mean((rr(.X, .theta))*(rr(.cft.X,theta)), weights)-.RO*.RC)
    
    .var     <-  .varRO + .varRC - 2*.covRORC
    return(.var)
  }
  
  #Get expected value and variance of that
  .logmeanvec   <- rep(NA, .nsim)
  .logvarvec    <- rep(NA, .nsim)
  .thetasim     <- mvrnorm(.nsim, thetahat, thetavar)
  for (i in 1:.nsim){
    .logmeanvec[i]  <- .logpifexp(.thetasim[i,])
    .logvarvec[i]   <- .logpifvar(.thetasim[i,])
  }
  
  #Get variance of that
  .logvarpif <- var(.logmeanvec) + mean(.logvarvec)
  
  #Create the confidence intervals
  .zqrt       <- .Z*sqrt(.logvarpif)

  #Compute the pif intervals
  .cipif         <- 1-c("Lower" = .inverse*exp(.zqrt), "Point_Estimate" =  .inverse, "Upper" = .inverse*exp(-.zqrt) )
  
  #Return variance
  return(.cipif)
  
}
