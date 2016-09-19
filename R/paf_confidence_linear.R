#' @title Approximate Confidence Intervals for the Potential Impact Fraction
#' 
#' @description Function that calculates approximate confidence intervals to the potential impact fraction
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
#' #Example with theta and X multivariate
#' set.seed(18427)
#' X1 <- rnorm(2000)
#' X2 <- rnorm(2000)
#' X  <- as.matrix(cbind(X1,X2))
#' thetahat <- c(0.1, 0.03)
#' thetasd <- matrix(c(0.1, 0, 0, 0.05), byrow = T, nrow = 2)
#' paf.confidence.linear(X, thetahat, thetasd, function(X, theta){exp(theta[1]*X[,1] + theta[2]*X[,2])}) 
#' 
#' @export


paf.confidence.linear <- function(X, thetahat, thetasd, rr, weights =  rep(1/nrow(as.matrix(X)),nrow(as.matrix(X))), confidence = 95, nsim = 1000){
  
  #Set the vector for the confidence intervals
  .ci <- c("Lower_CI" = NA, "Point_Estimate" = NA, "Upper_CI" = NA, "Estimated_Variance" = NA)
  
  #Set Z for the confidence interval
  Z <- qnorm(1 - ((100-confidence)/200))
  
  #Get the point estimate and variance
  .ci["Point_Estimate"]     <- pif(X, thetahat, rr, weights = weights, eval.cvx = FALSE)
  .ci["Estimated_Variance"] <- paf.variance.linear(X, thetahat, thetasd, rr, weights = weights, nsim = nsim)
  .ci["Lower_CI"]           <- max(.ci["Point_Estimate"] - Z*sqrt(.ci["Estimated_Variance"]), 0)
  .ci["Upper_CI"]           <- min(.ci["Point_Estimate"] + Z*sqrt(.ci["Estimated_Variance"]), 1)
  
  return(.ci)
  
}


