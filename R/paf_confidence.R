#' @title Confidence Intervals for the Population Attributable Fraction
#' 
#' @description Function that calculates confidence intervals of the Population Attributable Fraction
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
#' @param method     Method of estimation of confidence intervals "inverse", "log", "linear", "one2one"
#'                   (default: inverse)
#'                   
#' @param thetalow   Minimum value of theta (needed if \code{method} = one2one)
#' 
#' @param thetaup   Maximum value of theta (needed if \code{method} = one2one)
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
#' thetavar <- 0.1
#' paf.confidence(X, thetahat, thetavar, function(X, theta){exp(theta*X)})
#' 
#' #With larger sample the variance reduces
#' set.seed(18427)
#' X <- rnorm(10000,3,.7)
#' thetahat <- 0.12
#' thetavar <- 0.1
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
#' thetasd <- matrix(c(0.1, 0, 0, 0.4), byrow = TRUE, nrow = 2)
#' rr <- function(X, theta){exp(theta[1]*X[,1] + theta[2]*X[,2])}
#' paf.confidence(X, thetahat, thetasd, rr) 
#' 
#' @import MASS
#' @export

paf.confidence <- function(X, thetahat, thetavar = NA, rr, weights =  rep(1/nrow(as.matrix(X)),nrow(as.matrix(X))),
                           thetalow = NA, thetaup = NA, nsim = 1000, confidence = 95, 
                           method = c("inverse", "log", "linear", "one2one"), 
                                   force.min = FALSE){
  
  #Get method from vector
  .method   <- as.vector(method)[1]
  
  #Change variance to matrix
  .thetavar <- as.matrix(thetavar)
  
  #Check that thetas are as they should
  check.thetas(.thetavar, thetahat, thetalow, thetaup, .method)
  
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
           .cipaf <- paf.confidence.linear(X, thetahat, .thetavar, rr, weights, confidence, nsim,
                                           check_thetas = FALSE)
         },
         
         one2one = {
           .cipaf <- paf.confidence.one2one(X, thetahat, thetalow, thetaup, rr, weights, confidence,
                                            check_thetas = FALSE)
         },
         
         {
           warning("Invalid method. Defaulting to method = inverse with no forcing")
           .cipaf <- paf.confidence.inverse(X, thetahat, thetavar, rr, weights, nsim, confidence, FALSE)
         }
  )
  
  
  #Return variance
  return(.cipaf)
  
}