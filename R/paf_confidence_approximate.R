#' @title Approximate Confidence Intervals for the Population Attributable Fraction
#' 
#' @description Function that calculates approximate confidence intervals to the population attributable fraction
#' 
#' @param Xmean     Mean value of exposure levels from a previous study.
#' 
#' @param Xvar      Variance of exposure levels from a previous study.
#' 
#' @param thetahat  Estimative of \code{theta} for the Relative Risk function
#' 
#' @param thetavar   Estimator of standard error of thetahat (usually standard error) 
#' 
#' @param rr        Function for relative risk
#' 
#' 
#' **Optional**
#' 
#' @param nsim      Number of simulations for estimation of variance
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
#' #Example 1
#' set.seed(46987)
#' rr      <- function(X,theta){exp(X*theta)}
#' Xmean   <- 3
#' Xvar    <- 1
#' theta   <- 0.4
#' thetasd <- 0.001
#' .Xmean  <- as.matrix(Xmean)
#' .Xvar   <- as.matrix(Xvar)
#' paf.confidence.approximate(Xmean,Xvar,theta,thetasd,rr)
#' 
#' #Example 2: Compare paf.variance.approximate with paf.variance.linear
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
#' paf.confidence.approximate(Xmean,Xvar,theta,thetasd,rr)
#' paf.confidence.linear(X, theta, thetasd, rr)
#' 
#' @export


paf.confidence.approximate <- function(Xmean, Xvar, thetahat, thetavar, rr,
                                  confidence = 95, nsim = 1000, check_thetas = TRUE){
  
  #Make thetavar matrix
  .thetavar <- as.matrix(thetavar)
  
  #Function checking confidence is correctly specified
  check.confidence(confidence)
  
  #Function for checking that thetas are correctly inputed
  if(check_thetas){ check.thetas(.thetavar, thetahat, NA, NA, "approximate") }
  
  #Set the vector for the confidence intervals
  .ci <- c("Lower_CI" = NA, "Point_Estimate" = NA, "Upper_CI" = NA, "Estimated_Variance" = NA)
  
  #Set Z for the confidence interval
  Z <- qnorm(1 - ((100-confidence)/200))
  
  #Get the point estimate and variance
  .ci["Point_Estimate"]     <- pif.approximate(Xmean = Xmean, Xvar = Xvar, thetahat = thetahat, rr = rr)
  .ci["Estimated_Variance"] <- paf.variance.approximate(Xmean = Xmean, Xvar = Xvar, thetahat = thetahat, thetasd = .thetavar, rr = rr, nsim = nsim)
  .ci["Lower_CI"]           <- .ci["Point_Estimate"] - Z*sqrt(.ci["Estimated_Variance"])
  .ci["Upper_CI"]           <- .ci["Point_Estimate"] + Z*sqrt(.ci["Estimated_Variance"])
  
  return(.ci)
  
}


