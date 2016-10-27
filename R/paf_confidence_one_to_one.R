#' @title Approximate Confidence Intervals for the Population Attributable Fraction with one to one RR function, only for unidimensional theta values
#' 
#' @description Function that calculates approximate confidence intervals of the Population Attributable Fraction
#' considering a one to one Relative Risk with unidimensional theta values.
#' 
#' @param X         Random sample (can be vector or matrix) which includes exposure and covariates. Or mean exposure from a previous study if no sample is available.
#' 
#' @param thetahat  Estimative of \code{theta} for the Relative Risk function
#' 
#' @param thetalow  Lower bound of the confidence interval (vector)
#' 
#' @param thetaup   Upper band of the confidence interval (vector)
#' 
#' @param rr        Function for relative risk
#' 
#' 
#' **Optional**
#' 
#' @param weights   Survey \code{weights} for the random sample \code{X}
#' 
#' @param method    Either \code{empirical} (default) or \code{approximate}. 
#' 
#' @param Xvar      Variance of exposure levels.
#' 
#' @param confidence Confidence level \% (default: 95)
#' 
#' @param check_thetas Check that thetas are correctly specified
#' 
#' @author Rodrigo Zepeda Tello \email{rodrigo.zepeda@insp.mx}
#' @author Dalia Camacho García Formentí 
#' 
#' @examples 
#' 
#' #Example with risk given by HR
#' set.seed(18427)
#' X <- rnorm(1000, 3,.7)
#' thetahat <- 0.4
#' thetalow <- 0.1
#' thetaup  <- 0.7
#' paf.confidence.one2one(X, thetahat, thetalow, thetaup, function(X, theta){exp(theta*X)})
#' 
#' 
#' #Example with theta and X multivariate
#' set.seed(18427)
#' X1 <- rnorm(1000,3,.7)
#' X2 <- rnorm(1000,3,.7)
#' X  <- as.matrix(cbind(X1,X2))
#' thetahat <- c(0.12, 0.03)
#' thetalow <- c(0.05, 0.01)
#' thetaup  <- c(0.25, 0.06)
#' rr <- function(X, theta){exp(theta[1]*X[,1] + theta[2]*X[,2])}
#' paf.confidence.one2one(X, thetahat, thetalow, thetaup, rr) 
#' 
#' #Example with approximate method
#' set.seed(46987)
#' rr      <- function(X,theta){exp(X*theta)}
#' X       <- rnorm(100,3.2,1)
#' Xmean   <- 3.2
#' Xvar    <- 1
#' theta   <- 0.4
#' thetasd <- 0.001
#' .Xmean  <- as.matrix(Xmean)
#' .Xvar   <- as.matrix(Xvar)
#' paf.confidence.one2one(Xmean, thetahat = .4, thetalow = .3, thetaup = .5, rr = rr, method = "approximate", Xvar = Xvar)

#' @import MASS
#' @export

paf.confidence.one2one <- function(X, thetahat, thetalow, thetaup, rr, 
                                   weights =  rep(1/nrow(as.matrix(X)),nrow(as.matrix(X))), 
                                   confidence = 95, check_thetas = TRUE,
                                   method = c("empirical","approximate"), Xvar = var(X)){
  
  #Get method from vector
  .method <- as.vector(method)[1]
  
  #Check that thetas apply
  if(check_thetas){ check.thetas(NA, thetahat, thetalow, thetaup, "one2one") }
  
  #Get thetavar as matrix
  thetavar <- matrix(0, ncol = length(thetahat), nrow = length(thetahat))
  
  #Calculate the PIF with confidence intervals
  switch (.method,
    empirical = {
    .upper <- paf.confidence.inverse(X, thetaup,  thetavar = thetavar, rr = rr, method = "empirical",
                                                  weights = weights, confidence = confidence, nsim = 0)
    .lower <- paf.confidence.inverse(X, thetalow, thetavar = thetavar, rr = rr, method = "empirical",
                                     weights = weights, confidence = confidence, nsim = 0)
    .point <- pif(X, thetahat, rr = rr, weights = weights)
    },
    approximate ={
      .upper <- paf.confidence.inverse(X, thetahat = thetaup, thetavar = thetavar, rr = rr, nsim = 0, 
                                       confidence = confidence, 
                                       method = "approximate", Xvar = Xvar)
      
      .lower <- paf.confidence.inverse(X, thetahat = thetalow, thetavar = thetavar, rr = rr, nsim = 0, 
                                       confidence = confidence, 
                                       method = "approximate", Xvar = Xvar)
      .point <- pif(X, thetahat, rr = rr, method = "approximate", Xvar = Xvar)
    }
  )
  
  
  #Return
  .confint <- c("Lower" = .lower["Lower"], "Point" = .point, "Upper" = .upper["Upper"])
  
  #Return
  return(.confint)
}
  
  