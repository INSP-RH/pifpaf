#' @title Confidence Intervals for the Population Attributable Fraction
#' 
#' @description Function that calculates confidence intervals of the Population Attributable Fraction
#' 
#' @param X          Random sample (can be vector or matrix) which includes exposure and covariates. Or sample mean if approximate method is selected.
#' 
#' @param thetahat   Estimative of \code{theta} for the Relative Risk function
#' 
#' @param thetavar   Estimator of standard deviation of thetahat (usually standard error)
#' 
#' @param rr         Function for relative risk
#' 
#' 
#' **Optional**
#' 
#' @param weights    Survey \code{weights} for the random sample \code{X}
#' 
#' @param nsim       Number of simulations
#' 
#' @param confidence Confidence level \% (default 95)
#' 
#' @param method     Method of estimation of confidence intervals "inverse", "log", "linear", "one2one"
#'                   (default: inverse)
#'                   
#' @param est.method Either \code{empirical} (default), \code{kernel} or \code{approximate}.
#' 
#' @param Xvar       Variance of exposure levels.
#'                   
#' @param thetamin   Minimum value of theta (needed if \code{method} = one2one)
#' 
#' @param thetamax   Maximum value of theta (needed if \code{method} = one2one)
#' 
#' @param force.min Boolean indicating whether to force the PAF to have a 
#'                  minimum value of 0 instead of allowing negative values (not recommended).
#'                  Only valid if \code{method} = "inverse"
#' 
#' @author Rodrigo Zepeda Tello \email{rodrigo.zepeda@insp.mx}
#' @author Dalia Camacho García Formentí 
#' 
#' @examples 
#' 
#' #Example with risk given by HR
#' set.seed(18427)
#' X <- rnorm(100,4,1)
#' thetahat <- 0.4
#' thetavar <- 0.01
#' paf.confidence(X, thetahat, thetavar, function(X, theta){exp(theta*X)})
#' 
#' #With larger sample the variance reduces
#' set.seed(18427)
#' X <- rnorm(10000,3,.7)
#' thetahat <- 0.12
#' thetavar <- 0.002
#' paf.confidence(X, thetahat, thetavar, function(X, theta){exp(theta*X)})
#' 
#' #We can force PAF's CI to be >= 0 
#' paf.confidence(X, thetahat, thetavar, function(X, theta){exp(theta*X)},
#'  method = "inverse", force.min = TRUE)
#' 
#' #Example with theta and X multivariate
#' set.seed(18427)
#' X1 <- rnorm(1000,3,.5)
#' X2 <- rnorm(1000,3,.6)
#' X  <- as.matrix(cbind(X1,X2))
#' thetahat <- c(0.12, 0.03)
#' thetasd <- matrix(c(0.001, 0, 0, 0.0024), byrow = TRUE, nrow = 2)
#' rr <- function(X, theta){exp(theta[1]*X[,1] + theta[2]*X[,2])}
#' paf.confidence(X, thetahat, thetasd, rr) 
#' 
#' #Examples with approximate method
#' set.seed(46987)
#' rr      <- function(X,theta){exp(X*theta)}
#' Xmean   <- 3
#' Xvar    <- 1
#' theta   <- 0.4
#' thetasd <- 0.001
#' .Xmean  <- as.matrix(Xmean)
#' .Xvar   <- as.matrix(Xvar)
#' paf.confidence(X = Xmean, Xvar = Xvar, thetahat = theta, thetavar = thetasd,rr = rr, est.method = "approximate")
#' 
#' #Compare approximate method with default method
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
#' paf.confidence(X = X, thetahat = theta, thetavar = thetasd, rr = rr)
#' paf.confidence(X = Xmean, thetahat = theta, thetavar = thetasd, rr = rr, Xvar = Xvar, est.method = "approximate")
#' 
#' @import MASS
#' @export

paf.confidence <- function(X, thetahat, thetavar = NA, rr,  Xvar = var(X), weights =  rep(1/nrow(as.matrix(X)),nrow(as.matrix(X))),
                           thetamin = NA, thetamax = NA, nsim = 1000, confidence = 95, 
                           method = c("inverse", "log", "linear", "one2one"),
                           est.method = c("empirical", "kernel", "approximate"),
                                   force.min = FALSE){
  
  #Get method from vector
  .method     <- as.vector(method)[1]
  .est.method <- as.vector(est.method)[1]
  
  #Change variance to matrix
  .thetavar <- as.matrix(thetavar)
  
  #Check that thetas are as they should
  check.thetas(.thetavar, thetahat, thetamin, thetamax, .method)
  
  if(est.method[1] == "approximate"){
    .cipaf <- paf.confidence.approximate(X, Xvar = Xvar, thetahat, .thetavar, rr, confidence, nsim,
                                         check_thetas = FALSE)
  }else{
    switch(.method,
           
           inverse = {
             .cipaf <- paf.confidence.inverse(X, thetahat, .thetavar, rr, weights, nsim, confidence, force.min,
                                              check_thetas = FALSE)
           }, 
           
           log = {
             .cipaf <- paf.confidence.loglinear(X, thetahat, .thetavar, rr, weights, nsim, confidence,
                                                check_thetas = FALSE)
           },
           
           linear = {
             .cipaf <- paf.confidence.linear(X, thetahat, .thetavar, rr, weights, nsim, confidence,
                                             check_thetas = FALSE)
           },
           
           one2one = {
             .cipaf <- paf.confidence.one2one(X, thetahat, thetamin, thetamax, rr, weights,confidence = confidence, nsim = nsim,
                                              check_thetas = FALSE)
           },
           
           {
             warning("Invalid method. Defaulting to method = inverse with no forcing")
             .cipaf <- paf.confidence.inverse(X, thetahat, thetavar, rr, weights, nsim, confidence, FALSE)
           }
    )
  }
 
  
  
  #Return variance
  return(.cipaf)
  
}