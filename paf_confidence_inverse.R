#' @title Approximate Confidence Intervals for the Population Attributable Fraction
#' 
#' @description Function that calculates approximate confidence intervals of the Population Attributable Fraction
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
#' @param nsim      Number of simulations
#' 
#' @param confidence Confidence level \% (default 95)
#' 
#' @param force.min Boolean indicating whether to force the PAF to have a 
#'                  minimum value of 0 instead of allowing negative values (not recommended).
#' 
#' @author Rodrigo Zepeda Tello \email{rodrigo.zepeda@insp.mx}
#' @author Dalia Camacho García Formentí 
#' 
#' @examples 
#' 
#' #Example with risk given by HR
#' set.seed(18427)
#' X <- rnorm(100)
#' thetahat <- 0.4
#' thetavar <- 0.1
#' paf.confidence.inverse(X, thetahat, thetavar, function(X, theta){exp(theta*X)})
#' 
#' #With larger sample the variance reduces
#' set.seed(18427)
#' X <- rnorm(10000)
#' thetahat <- 0.12
#' thetavar <- 0.1
#' paf.confidence.inverse(X, thetahat, thetavar, function(X, theta){exp(theta*X)})
#' 
#' #We can force PAF's CI to be >= 0 
#' paf.confidence.inverse(X, thetahat, thetavar, function(X, theta){exp(theta*X)}, force.min = TRUE)
#' 
#' #Example with theta and X multivariate
#' set.seed(18427)
#' X1 <- rnorm(1000)
#' X2 <- rnorm(1000)
#' X  <- as.matrix(cbind(X1,X2))
#' thetahat <- c(0.12, 0.03)
#' thetasd <- matrix(c(0.1, 0, 0, 0.4), byrow = TRUE, nrow = 2)
#' rr <- function(X, theta){exp(theta[1]*X[,1] + theta[2]*X[,2])}
#' paf.confidence.inverse(X, thetahat, thetasd, rr) 
#' 
#' @import MASS
#' @export

paf.confidence.inverse <- function(X, thetahat, thetasd, rr, weights =  rep(1/nrow(as.matrix(X)),nrow(as.matrix(X))),
                                     nsim = 1000, confidence = 95, force.min = FALSE){
  
  #Compute the PAF intervals
  .cipaf         <- 1-1/risk.ratio.confidence(X = X, thetahat = thetahat, 
                                              thetasd = thetasd, rr = rr, 
                                              weights =  weights, nsim = nsim, 
                                              confidence = confidence, force.min = force.min)
  
  #Return variance
  return(.cipaf)
  
}