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
#' @param cft       Function \code{cft(X)} for counterfactual. Leave empty for 
#'                  the Population Attributable Fraction \code{PAF} where counterfactual is 0 exposure
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
#' cft     <- function(X){0.5*X}
#' Xmean   <- 3
#' Xvar    <- 1
#' theta   <- 0.4
#' thetasd <- 0.001
#' .Xmean  <- as.matrix(Xmean)
#' .Xvar   <- as.matrix(Xvar)
#' pif.confidence.approximate(Xmean,Xvar,theta,thetasd,rr)
#' paf.confidence.approximate(Xmean,Xvar,theta,thetasd,rr)
#' pif.confidence.approximate(Xmean,Xvar,theta,thetasd,rr, cft)
#'
#' #Example 2: Compare pif.variance.approximate with pif.variance.linear
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
#' pif.confidence.approximate(Xmean, Xvar, theta, thetasd, rr)
#' paf.confidence.approximate(Xmean, Xvar, theta, thetasd, rr)
#' pif.confidence.approximate(Xmean, Xvar, theta, thetasd, rr, cft = function(X){sqrt(X)})
#' pif.confidence.approximate(Xmean, Xvar, theta, thetasd, rr, cft = function(X){0.5*X})
#' 
#' @export


pif.confidence.approximate <- function(Xmean, Xvar, thetahat, thetavar, rr,
                                       cft = function(Xmean){matrix(0,ncol = ncol(as.matrix(Xmean)), nrow = nrow(as.matrix(Xmean)))}, 
                                       confidence = 95, nsim = 1000, check_thetas = TRUE){
  check.confidence(confidence)
  
  #Make thetavar matrix
  .thetavar <- as.matrix(thetavar)
  
  #Function for checking that thetas are correctly inputed
  if(check_thetas){ check.thetas(.thetavar, thetahat, NA, NA, "approximate") }
  
  #Set the vector for the confidence intervals
  .ci <- c("Lower_CI" = NA, "Point_Estimate" = NA, "Upper_CI" = NA, "Estimated_Variance" = NA)
  
  #Set Z for the confidence interval
  Z <- qnorm(1 - ((100-confidence)/200))
  
  #Get the point estimate and variance
  .ci["Point_Estimate"]     <- pif.approximate(Xmean = Xmean, Xvar = Xvar, thetahat = thetahat, rr = rr, cft = cft)
  .ci["Estimated_Variance"] <- pif.variance.approximate.linear(Xmean = Xmean, Xvar = Xvar, thetahat = thetahat, thetasd = .thetavar, rr = rr, nsim = nsim)
  .ci["Lower_CI"]           <- .ci["Point_Estimate"] - Z*sqrt(.ci["Estimated_Variance"])
  .ci["Upper_CI"]           <- .ci["Point_Estimate"] + Z*sqrt(.ci["Estimated_Variance"])
  
  return(.ci)
  
}


