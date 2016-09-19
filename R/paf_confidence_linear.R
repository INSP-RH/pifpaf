#' @title Approximate Confidence Intervals for the Potential Impact Fraction
#' 
#' @description Function that calculates confidence intervals to the potential impact fraction
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
#' 
#' @param weights   Survey \code{weights} for the random sample \code{X}
#' 
#' @author Rodrigo Zepeda Tello \email{rodrigo.zepeda@insp.mx}
#' @author Dalia Camacho García Formentí 
#' 
#' @examples 
#' 
#' #Example with risk given by HR
#' set.seed(18427)
#' X <- rnorm(100)
#' thetahat <- 0.12
#' thetavar <- 0.1
#' paf.confidence.linear(X, thetahat, thetavar, function(X, theta){exp(theta*X)})
#' 
#' #With larger sample the confidence interval reduces
#' set.seed(18427)
#' X <- rnorm(10000)
#' thetahat <- 0.12
#' thetavar <- 0.1
#' paf.confidence.linear(X, thetahat, thetavar, function(X, theta){exp(theta*X)})
#' 
#' 
#' @import matrixStats
#' 
#' @export


paf.confidence.linear <- function(X, thetahat, thetasd, rr, weights =  rep(1/nrow(as.matrix(X)),nrow(as.matrix(X))), confidence = 95){
  
  #Set the vector for the confidence intervals
  .ci <- c("Lower_CI" = NA, "Point_Estimate" = NA, "Upper_CI" = NA, "Estimated_Variance" = NA)
  
  #Set Z for the confidence interval
  Z <- qnorm(1 - ((100-confidence)/200))
  
  #Get the point estimate and variance
  .ci["Point_Estimate"]     <- pif(X, thetahat, rr, weights = weights)
  .ci["Estimated_Variance"] <- paf.variance.linear(X, thetahat, thetasd, rr, weights = weights)
  .ci["Lower_CI"]           <- max(0,.ci["Point_Estimate"] - Z*sqrt(.ci["Estimated_Variance"]))
  .ci["Upper_CI"]           <- min(1,.ci["Point_Estimate"] + Z*sqrt(.ci["Estimated_Variance"]))
  
  return(.ci)
  
}


