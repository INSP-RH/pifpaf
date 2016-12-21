#' @title Approximate Confidence Intervals for the Population Attributable Fraction
#' 
#' @description Function that calculates approximate confidence intervals to the population attributable fraction.
#' 
#' @param X         Random sample (can be vector or matrix) which includes exposure and covariates.
#' 
#' @param thetahat  Estimative of \code{theta} for the Relative Risk function
#' 
#' @param thetavar   Estimator of variance of thetahat 
#' 
#' @param rr        Function for relative risk
#' 
#' 
#' **Optional**
#' 
#' @param Xvar      Covariance matrix of X or variance value (for univariate X)
#' 
#' @param method    Method for estimation (either approximate or empirical)
#' 
#' @param nsim      Number of simulations for estimation of variance
#' 
#' @param weights    Survey \code{weights} for the random sample \code{X}
#' 
#' @param confidence Confidence level \% (default 95)
#' 
#' @param check_thetas Checks that theta parameters are correctly inputed
#' 
#' @author Rodrigo Zepeda Tello \email{rodrigo.zepeda@insp.mx}
#' @author Dalia Camacho García Formentí 
#' 
#' @examples 
#' 
#' #Example with risk given by HR
#' set.seed(18427)
#' X <- rnorm(100,3,.5)
#' thetahat <- 0.12
#' thetavar <- 0.1
#' paf.confidence.linear(X, thetahat, thetavar, function(X, theta){exp(theta*X)})
#' 
#' #Same example with approximate method 
#' paf.confidence.linear(3, thetahat, thetavar, function(X, theta){exp(theta*X)}, 
#' Xvar = 0.25, method = "approximate")
#' 
#' #Example with theta and X multivariate
#' set.seed(18427)
#' X1 <- rnorm(2000, 3,.5)
#' X2 <- rnorm(2000,3,.5)
#' X  <- as.matrix(cbind(X1,X2))
#' thetahat <- c(0.1, 0.03)
#' thetavar <- matrix(c(0.1, 0, 0, 0.05), byrow = TRUE, nrow = 2)
#' rr <- function(X, theta){exp(theta[1]*X[,1] + theta[2]*X[,2])}
#' paf.confidence.linear(X, thetahat, thetavar, rr) 
#' 
#' 
#' 
#' @export


paf.confidence.linear <- function(X, thetahat, thetavar, rr, Xvar = var(X),
                                  weights =  rep(1/nrow(as.matrix(X)),nrow(as.matrix(X))), 
                                  confidence = 95, nsim = 1000, check_thetas = TRUE,
                                  method = c("empirical", "kernel", "approximate")){
  
  #Make thetavar matrix
  .thetavar <- as.matrix(thetavar)
  
  #Get the method
  .method   <- method[1]
  
  #Get the point estimate and variance
  switch(.method,
         empirical = {
           
           #Check confidence
           check.confidence(confidence)
           
           #Set Z for the confidence interval
           Z <- qnorm(1 - ((100-confidence)/200))
           
           #Set the vector for the confidence intervals
           .ci <- c("Lower_CI" = NA, "Point_Estimate" = NA, "Upper_CI" = NA, "Estimated_Variance" = NA)
           
           #Create confidence intervals with asymptotic normality
           .ci["Point_Estimate"]      <- pif(X, thetahat, rr, weights = weights, method = .method)
           .ci["Estimated_Variance"]  <- paf.variance.linear(X, thetahat, .thetavar, rr, weights = weights, nsim = nsim, check_thetas)
           .ci["Lower_CI"]            <- .ci["Point_Estimate"] - Z*sqrt(.ci["Estimated_Variance"])
           .ci["Upper_CI"]            <- .ci["Point_Estimate"] + Z*sqrt(.ci["Estimated_Variance"])
           
         }, approximate = {
           .ci <- paf.confidence.approximate(Xmean = X, Xvar = Xvar, thetahat = thetahat, 
                                             thetavar = .thetavar, rr = rr, 
                                             confidence = confidence, nsim = nsim, 
                                             check_thetas = check_thetas)
         })
  
  return(.ci)
  
}


