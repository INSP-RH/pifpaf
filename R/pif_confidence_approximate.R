#' @title Approximate Confidence Intervals for the Population Attributable Fraction
#' 
#' @description Function that calculates approximate confidence intervals to the population attributable fraction
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
#' @param is_paf Force evaluation of paf  
#'  
#' @author Rodrigo Zepeda-Tello \email{rzepeda17@gmail.com}
#' @author Dalia Camacho-García-Formentí \email{daliaf172@gmail.com}
#' 
#' @examples 
#' \dontrun{
#' #Example 1: Exponential Relative risk
#' #--------------------------------------------
#' set.seed(46987)
#' rr      <- function(X,theta){exp(X*theta)}
#' cft     <- function(X){0.5*X}
#' X       <- runif(1000)
#' Xmean   <- data.frame(mean(X))
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
#' Xmean    <- data.frame(t(colMeans(X)))
#' Xvar     <- cov(X)
#' theta    <- c(0.12, 0.17)
#' thetavar  <- matrix(c(0.001, 0.00001, 0.00001, 0.004), byrow = TRUE, nrow = 2)
#' rr       <- function(X, theta){exp(theta[1]*X[,1] + theta[2]*X[,2])}
#' pif.confidence.approximate(Xmean, Xvar, theta, thetavar, rr, 
#' cft = function(X){cbind(0.5*X[,1],0.4*X[,2])})
#' }
#' @importFrom stats qnorm
#' @keywords internal
#' @export


pif.confidence.approximate <- function(Xmean, Xvar, thetahat, thetavar, rr, 
                                       cft = function(Xmean){matrix(0,ncol = ncol(as.matrix(Xmean)), nrow = nrow(as.matrix(Xmean)))}, 
                                       check_thetas = TRUE, check_cft = TRUE, check_xvar = TRUE, check_rr = TRUE, 
                                       check_integrals = TRUE, check_exposure = TRUE, deriv.method.args = list(), 
                                       deriv.method = c("Richardson", "complex"), nsim = 1000, confidence = 95,
                                       is_paf = FALSE){
  
  #Check confidence interval
  check.confidence(confidence)
  
  #Set the vector for the confidence intervals
  .ci <- c("Lower_CI" = NA, "Point_Estimate" = NA, "Upper_CI" = NA, "Estimated_Variance" = NA)
  
  #Set Z for the confidence interval
  Z <- qnorm(1 - ((100-confidence)/200))
  
  #Get the point estimate and variance
  .ci["Point_Estimate"]     <- pif.approximate(X = Xmean, Xvar = Xvar, thetahat = thetahat, 
                                               rr = rr, cft = cft, 
                                               deriv.method.args = deriv.method.args, 
                                               deriv.method = deriv.method,
                                               check_exposure = check_exposure, 
                                               check_rr = check_rr, 
                                               check_integrals = check_integrals, 
                                               is_paf = is_paf)
  
  .ci["Estimated_Variance"] <- 
    pif.variance.approximate.linear(X = Xmean, thetahat = thetahat, rr = rr,  
                                    thetavar = thetavar, Xvar = Xvar, cft = cft, 
                                    check_thetas = check_thetas, check_cft = check_cft, 
                                    check_xvar = check_xvar, check_rr = FALSE, 
                                    check_integrals = FALSE, check_exposure = FALSE,
                                    deriv.method.args = deriv.method.args,   
                                    deriv.method = deriv.method, nsim = nsim, is_paf = is_paf)
  
  #Get lower and upper ci (asymptotically normal)
  .ci["Lower_CI"]           <- .ci["Point_Estimate"] - Z*sqrt(.ci["Estimated_Variance"])
  .ci["Upper_CI"]           <- .ci["Point_Estimate"] + Z*sqrt(.ci["Estimated_Variance"])
  
  .ci["Upper_CI"]  <- min(.ci["Upper_CI"] , 1)
  
  return(.ci)
  
}


