#' @title Confidence Intervals for the Potential Impact Fraction when only mean and variance of exposure values is available using loglinear method
#' 
#' @description Confidence intervals for the Population Attributable Fraction for the approximate method where only mean and variance from a previous study is available.For relative risk inyective functions, the pif is inyective, and intervals can be calculated for log(pif), and then transformed to pif CI.
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
#' @param cft       Function \code{cft(X)} for counterfactual. Leave empty for 
#'                  the Population Attributable Fraction \code{PAF} where counterfactual is 0 exposure
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
#' cft     <- function(X){0.4*X}
#' Xmean   <- 3
#' Xvar    <- 1
#' theta   <- 0.4
#' thetasd <- 0.001
#' pif.confidence.approximate.loglinear(Xmean, Xvar, theta, thetasd, rr)
#' paf.confidence.approximate.loglinear(Xmean, Xvar, theta, thetasd, rr)
#' pif.confidence.approximate.loglinear(Xmean, Xvar, theta, thetasd, rr, cft)
#' pif.confidence.approximate.loglinear(Xmean, Xvar, theta, thetasd, rr, cft = function(X){sqrt(X)})
#'
#'#Example 2: Compare pif.variance.approximate with paf.variance.loglinear
#'X1       <- rnorm(1000,3,.5)
#'X2       <- rnorm(1000,4,1)
#'X        <- as.matrix(cbind(X1,X2))
#'Xmean    <- colMeans(X)
#'Xvar     <- cov(X)
#'thetahat <- c(0.12, 0.17)
#'thetasd  <- matrix(c(0.001, 0.00001, 0.00001, 0.004), byrow = TRUE, nrow = 2)
#'rr       <- function(X, theta){exp(theta[1]*X[,1] + theta[2]*X[,2])}
#'
#' pif.confidence.approximate.loglinear(Xmean, Xvar, thetahat, thetasd, rr)
#' paf.confidence.loglinear(X, thetahat, thetasd, rr)
#' pif.confidence.approximate.loglinear(Xmean, Xvar, thetahat, thetasd, rr, cft = function(X){0.8*X})
#' 
#' @import MASS numDeriv
#' @export

pif.confidence.approximate.loglinear <- function(Xmean, Xvar, thetahat, thetavar, rr,
                                                 cft = function(Varx){matrix(0,ncol = ncol(as.matrix(Varx)), nrow = nrow(as.matrix(Varx)))}, 
                                                 nsim = 1000, confidence = 95, check_thetas = TRUE){
  
  #Get confidence
  check.confidence(confidence)
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
  .logpifexp <- function(.theta){
    rr.cft.fun <- function(X){
      Xcft.value  <- cft(X)
      rr(Xcft.value,.theta)
    } 
    rr.fun.x <- function(X){
      rr(X,.theta)
    }
    
    #Estimate weighted sums
    .R1   <- rr.cft.fun(.Xmean) + 0.5*EntryMult(hessian(rr.cft.fun,.Xmean), .Xvar)
    .R0   <- rr.fun.x(.Xmean) + 0.5*EntryMult(hessian(rr.fun.x,.Xmean), .Xvar)
    
    return( log(.R1)-log(.R0) )
  }
  
  #Inverse
  .pif      <- pif.approximate(.Xmean, .Xvar, thetahat, rr, cft)
  .inverse  <- 1 - .pif
  
  #Calculate the conditional variance as a function of theta
  
  .logpifvar <- function(.theta){
    rr.fun.x <- function(X){
      rr(X,.theta)
    }
    
    rr.cft.fun <- function(X){
      Xcft.value  <- cft(X)
      rr(Xcft.value,.theta)
    } 
    dR1    <- as.matrix(grad(rr.cft.fun, .Xmean))
    R1     <- rr(cft(.Xmean), .theta)
    
    dR0    <- as.matrix(grad(rr.fun.x, .Xmean))
    R0     <- rr(.Xmean, .theta)
    
    .var1  <- t(1/R1*dR1)%*%.Xvar%*%(1/R1*dR1)
    .var0  <- t(1/R0*dR0)%*%.Xvar%*%(1/R0*dR0)
    .var   <- .var0 +.var1 
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