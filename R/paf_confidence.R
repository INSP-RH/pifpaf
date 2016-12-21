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
#' X <- rnorm(1000,3,.7)
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
#' X1 <- rnorm(100,3,.5)
#' X2 <- rnorm(100,3,.6)
#' X  <- as.matrix(cbind(X1,X2))
#' thetahat <- c(0.12, 0.03)
#' thetavar <- matrix(c(0.001, 0, 0, 0.0024), byrow = TRUE, nrow = 2)
#' rr <- function(X, theta){exp(theta[1]*X[,1] + theta[2]*X[,2])}
#' paf.confidence(X, thetahat, thetavar, rr) 
#' 
#' #Examples with approximate method
#' set.seed(46987)
#' rr       <- function(X,theta){exp(X*theta)}
#' Xmean    <- 3
#' Xvar     <- 1
#' theta    <- 0.4
#' thetavar <- 0.001
#' .Xmean   <- as.matrix(Xmean)
#' .Xvar    <- as.matrix(Xvar)
#' paf.confidence(X = Xmean, Xvar = Xvar, thetahat = theta,
#'  thetavar = thetavar, rr = rr, est.method = "approximate")
#' 
#' #Compare approximate method with default method
#' X1       <- rnorm(1000,3,.5)
#' X2       <- rnorm(1000,4,1)
#' X        <- as.matrix(cbind(X1,X2))
#' Xmean    <- colMeans(X)
#' Xvar     <- cov(X)
#' theta    <- c(0.12, 0.17)
#' thetavar  <- matrix(c(0.001, 0.00001, 0.00001, 0.004), byrow = TRUE, nrow = 2)
#' rr       <- function(X, theta){exp(theta[1]*X[,1] + theta[2]*X[,2])}
#' paf.confidence(X = X, thetahat = theta, thetavar = thetavar, rr = rr)
#' paf.confidence(X = Xmean, thetahat = theta, thetavar = thetavar,
#'  rr = rr, Xvar = Xvar, est.method = "approximate")
#' 
#' #Examples with approximate methods
#' rr      <- function(X,theta){exp(X*theta)}
#' X       <- rnorm(100,3.2,1)
#' Xmean   <- 3.2
#' Xvar    <- 1
#' theta   <- 0.4
#' thetavar <- 0.001
#' paf.confidence(Xmean, thetahat = theta, thetavar = thetavar,
#'  rr = rr, est.method = "approximate", Xvar = Xvar, method = "linear")
#'  
#' paf.confidence(Xmean, thetahat = theta, thetavar = thetavar,
#'  rr=rr, est.method = "approximate", Xvar = Xvar, method = "log")
#'  
#' paf.confidence(Xmean, thetahat = theta, thetavar = thetavar, rr=rr, 
#' est.method = "approximate", Xvar = Xvar, method = "inverse")
#' 
#' paf.confidence(Xmean, thetahat = .4, thetamin = .3, thetamax = .5, 
#' rr = rr, est.method = "approximate", method = "one2one",  Xvar = Xvar)
#' 
#' @import MASS
#' @export

paf.confidence <- function(X, thetahat, thetavar = NA, rr,  Xvar = NA, weights =  rep(1/nrow(as.matrix(X)),nrow(as.matrix(X))),
                           thetamin = NA, thetamax = NA, nsim = 1000, confidence = 95, 
                           method = c("inverse", "log", "linear", "one2one"),
                           est.method = c("empirical", "kernel", "approximate"),
                           force.min = FALSE, check_thetas = TRUE){
  
  #Check confidence
  check.confidence(confidence)
  
  #Get method from vector
  .method     <- as.vector(method)[1]
  .est.method <- as.vector(est.method)[1]
  
  #Change variance to matrix
  .thetavar <- as.matrix(thetavar)
  
  #Check that thetas are as they should
  check.thetas(.thetavar, thetahat, thetamin, thetamax, .method)
    
  switch (.method,
          inverse = {
            .cipaf <- paf.confidence.inverse(X = X, thetahat = thetahat, thetavar = .thetavar, 
                                             weights = weights, rr = rr, confidence = confidence, 
                                             nsim = nsim, Xvar = Xvar, method = est.method, 
                                             force.min = force.min, check_thetas = check_thetas)
          }, log     = {
            .cipaf <- paf.confidence.loglinear(X = X, Xvar = Xvar, thetahat = thetahat, weights = weights,
                                               thetavar = thetavar, rr=rr, nsim = nsim, 
                                               confidence = confidence, check_thetas = check_thetas,
                                               method = est.method)
          }, linear  = {
            .cipaf <- paf.confidence.linear(X, Xvar = Xvar, thetahat = thetahat, thetavar = .thetavar, 
                                            rr = rr, weights = weights, confidence = confidence, 
                                            nsim = nsim, check_thetas = check_thetas, method = est.method)
          }, one2one = {
            .cipaf <- paf.confidence.one2one(X, thetahat = thetahat, thetalow = thetamin, thetaup = thetamax, 
                                             rr= rr,confidence = confidence, weights = weights, 
                                             check_thetas = check_thetas, method = est.method, Xvar = Xvar)
          }, {
            warning(paste("Invalid method. Please choose one of the following: 'inverse',",
                          "'log', 'linear' or 'one2one'"))
            .cipaf <- NA
          }
    )
  
  #Return variance
  return(.cipaf)
  
}