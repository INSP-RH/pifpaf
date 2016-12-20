#' @title Approximate Variance for the Potential Impact Fraction using the approximate method
#' 
#' @description Function that calculates approximate variance to the potential impact fraction.
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
#' @param cft       Function \code{cft(X)} for counterfactual. Leave empty for 
#'                  the Population Attributable Fraction \code{PAF} where counterfactual is 0 exposure.
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
#' #Example 1
#' set.seed(46987)
#' rr      <- function(X,theta){exp(X*theta)}
#' cft     <- function(X){0.5*X}
#' Xmean   <- 3
#' Xvar    <- 1
#' theta   <- 1.2
#' thetavar <- 0.15
#' pif.variance.approximate.linear(Xmean,Xvar,theta,thetavar,rr)
#' pif.variance.approximate.linear(Xmean,Xvar,theta,thetavar,rr,cft) 
#' pif.variance.approximate.linear(Xmean,Xvar,theta,thetavar,rr,cft = function(X){sqrt(X)}) 
#' 
#' 
#' #Example 2: Compare pif.variance.approximate with paf.variance.approximate
#' X1       <- rnorm(1000,3,.5)
#' X2       <- rnorm(1000,4,1)
#' X        <- as.matrix(cbind(X1,X2))
#' Xmean    <- colMeans(X)
#' Xvar     <- cov(X)
#' theta    <- c(0.12, 0.17)
#' thetavar  <- matrix(c(0.001, 0.00001, 0.00001, 0.004), byrow = TRUE, nrow = 2)
#' rr       <- function(X, theta){exp(theta[1]*X[,1] + theta[2]*X[,2])}
#' pif.variance.approximate.linear(Xmean,Xvar,theta,thetavar,rr)
#' paf.variance.approximate (Xmean,Xvar,theta,thetavar,rr)
#' pif.variance.approximate.linear(Xmean,Xvar,theta,thetavar,rr, cft)
#' 
#' @import MASS stats numDeriv matrixcalc
#' @export 


pif.variance.approximate.linear <- function(Xmean, Xvar, thetahat, thetavar, rr,
                                     cft = function(Xmean){matrix(0,ncol = ncol(as.matrix(Xmean)), nrow = nrow(as.matrix(Xmean)))}, 
                                     check_thetas = TRUE, 
                                     nsim = 1000){
  
  #Function for checking that thetas are correctly inputed
  if(check_thetas){ check.thetas(thetavar, thetahat, NA, NA, "approximate") }
  
  #Set X as matrix
  .Xmean  <- matrix(Xmean, ncol = length(Xmean))
  .Xvar   <- matrix(Xvar, ncol = sqrt(length(Xvar)))
  #Set a minimum for nsim
  .nsim        <- max(nsim,10)
  
  if(is.positive.semi.definite(.Xvar) == FALSE){
    stop("Cobariance matrix must be positive semi-definite.")
  }
  
  #Check exposure values are greater than zero
  check.exposure(.Xmean)
  
  #Check that rr is 1 when X = 0
  check.rr(.Xmean, thetahat, rr)
  
  #Get the expected pif
  .pifexp <- function(theta){
    pif.approximate(Xmean, Xvar, thetahat = theta, rr = rr, cft = cft)
  }
  
  #Get the variance of pif
  .pifvar <- function(theta){
    #Rewrite functions as functions of X
    rr.fun.x <- function(X){
      rr(X,theta)
    }
    rr.fun.cft <- function(X){
      cftX  <- cft(X)
      rr(cftX, theta)
    }
    
    dR0 <- as.matrix(grad(rr.fun.x, .Xmean))
    dR1 <- as.matrix(grad(rr.fun.cft, .Xmean))
    R0  <- rr(.Xmean, theta)
    R1  <- rr(cft(.Xmean), theta)
    aux <- ((dR1*R0 - dR0*R1)/(R0^2))
    vr  <- t(aux)%*%.Xvar%*%(aux)
    return(vr)
  }
  
  #Get expected value and variance of that
  .meanvec   <- rep(NA, .nsim)
  .varvec    <- rep(NA, .nsim)
  .thetasim  <- mvrnorm(.nsim, thetahat, thetavar, empirical = TRUE)
  for (i in 1:.nsim){
    .meanvec[i]  <- .pifexp(.thetasim[i,])
    .varvec[i]   <- .pifvar(.thetasim[i,])
  }
  
  #Get variance of that
  .varpif <- var(.meanvec) + mean(.varvec)
  
  #Return variance
  return(.varpif)
  
}


