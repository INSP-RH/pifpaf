#' @title Confidence Intervals for the Potential Impact Fraction for the approximate method
#' 
#' @description Function that calculates confidence intervals of the potential impact fraction for the approximate method
#' 
#' @param Xmean     Mean value of exposure levels.
#' 
#' @param Xvar      Variance of exposure levels.
#' 
#' @param thetahat  Estimator (vector or matrix) of \code{theta} for the 
#'   Relative Risk function \code{rr}
#' 
#' @param thetavar   Estimator of variance of \code{thetahat}
#' 
#' @param rr        Function for Relative Risk which uses parameter 
#'   \code{theta}. The order of the parameters shound be \code{rr(X, theta)}.
#' 
#' 
#' **Optional**
#' 
#' @param cft       Function \code{cft(X)} for counterfactual. Leave empty for 
#'   the Population Attributable Fraction \code{\link{paf}} where counterfactual
#'   is 0 exposure.
#' 
#' @param confidence Concidence level (0 to 100) default = \code{95} \%
#' 
#' @param nsim      Number of simulations for estimation of variance
#' 
#' @param check_thetas Checks that theta parameters are correctly inputed
#' 
#' @param check_cft  Check if counterfactual function \code{cft} reduces exposure.
#' 
#' @param check_xvar Check if it is covariance matrix.
#' 
#'@param check_exposure  Check that exposure \code{X} is positive and numeric
#'  
#'@param check_rr        Check that Relative Risk function \code{rr} equals 
#'  \code{1} when evaluated at \code{0}
#'  
#'@param deriv.method.args \code{method.args} for
#'  \code{\link[numDeriv]{hessian}}.
#'  
#'@param deriv.method      \code{method} for \code{\link[numDeriv]{hessian}}.
#'  Don't change this unless you know what you are doing.
#'  
#'@param check_integrals Check that counterfactual and relative risk's expected 
#'  values are well defined for this scenario.
#'  
 #' @author Rodrigo Zepeda Tello \email{rodrigo.zepeda@insp.mx}
#' @author Dalia Camacho García Formentí \email{daliaf172@gmail.com}
#' 
#' @examples 
#' 
#' #Example 1: Exponential Relative risk
#' #--------------------------------------------
#' set.seed(46987)
#' rr      <- function(X,theta){exp(X*theta)}
#' cft     <- function(X){0.5*X}
#' X       <- runif(1000)
#' Xmean   <- mean(X)
#' Xvar    <- var(X)
#' theta   <-  0.2
#' thetavar <- 0.015
#' pif.confidence.approximate(Xmean, Xvar, theta, thetavar, rr)
#' pif.confidence.approximate(Xmean, Xvar, theta, thetavar, rr, cft) 
#'
#' #Example 2: Multivariate example
#' #--------------------------------------------
#' X1       <- rnorm(1000,3,.5)
#' X2       <- rnorm(1000,4,1)
#' X        <- as.matrix(cbind(X1,X2))
#' Xmean    <- colMeans(X)
#' Xvar     <- cov(X)
#' theta    <- c(0.12, 0.17)
#' thetavar  <- matrix(c(0.001, 0.00001, 0.00001, 0.004), byrow = TRUE, nrow = 2)
#' rr       <- function(X, theta){exp(theta[1]*X[,1] + theta[2]*X[,2])}
#' pif.confidence.approximate(Xmean, Xvar, theta, thetavar, rr, 
#' cft = function(X){cbind(0.5*X[,1],0.4*X[,2])})
#' 
#' @importFrom stats qnorm
#' 
#' @export


pif.confidence.approximate <- function(Xmean, Xvar, thetahat, thetavar, rr, 
                                       cft = function(Xmean){matrix(0,ncol = ncol(as.matrix(Xmean)), nrow = nrow(as.matrix(Xmean)))}, 
                                       check_thetas = TRUE, check_cft = TRUE, check_xvar = TRUE, check_rr = TRUE, 
                                       check_integrals = TRUE, check_exposure = TRUE, deriv.method.args = list(), 
                                       deriv.method = c("Richardson", "complex"), nsim = 1000, confidence = 95){
  
  #Check confidence interval
  check.confidence(confidence)
  
  #Set the vector for the confidence intervals
  .ci <- c("Lower_CI" = NA, "Point_Estimate" = NA, "Upper_CI" = NA, "Estimated_Variance" = NA)
  
  #Set Z for the confidence interval
  Z <- qnorm(1 - ((100-confidence)/200))
  
  #Get the point estimate and variance
  .ci["Point_Estimate"]     <- pif.approximate(X = Xmean, Xvar = Xvar, thetahat = thetahat, rr = rr, cft = cft,
                                               deriv.method.args = deriv.method.args, deriv.method = deriv.method,
                                               check_exposure = check_exposure, check_rr = check_rr, 
                                               check_integrals = check_integrals)
  
  .ci["Estimated_Variance"] <- pif.variance.approximate.linear(Xmean = Xmean, Xvar = Xvar, thetahat = thetahat, 
                                                               thetavar = thetavar, rr = rr, cft = cft, 
                                                               check_thetas = check_thetas, check_cft = check_cft, 
                                                               check_xvar = check_xvar, deriv.method.args = deriv.method.args, 
                                                               check_rr = FALSE, check_integrals = FALSE, check_exposure = FALSE,  #False as pif.approximate already checked them
                                                               deriv.method = deriv.method, nsim = nsim)
  
  #Get lower and upper ci (asymptotically normal)
  .ci["Lower_CI"]           <- .ci["Point_Estimate"] - Z*sqrt(.ci["Estimated_Variance"])
  .ci["Upper_CI"]           <- .ci["Point_Estimate"] + Z*sqrt(.ci["Estimated_Variance"])
  
  #Check ci < 1
  if (.ci["Upper_CI"] >= 1){
    
    #Transform the problem to 0 <= 1 - pif to apply bounded CIs
    .transf_ci             <- 1 - .ci
    names(.transf_ci)      <- c("Upper_CI", "Point_Estimate", "Lower_CI", "Estimated_Variance")
    .transf_ci["Lower_CI"] <- (.transf_ci["Point_Estimate"]^2)/.transf_ci["Upper_CI"]           #Bound 1 - .ci below 
    .ci["Upper_CI"]        <- 1 - .transf_ci["Lower_CI"]                                        #Transform back
    
  }
  
  return(.ci)
  
}


