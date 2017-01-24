#' @title Approximate Confidence Intervals for the Population Attributable
#'   Fraction with one to one RR function, only for unidimensional theta values
#'   
#' @description Function that calculates approximate confidence intervals of the
#'   Population Attributable Fraction considering a one to one Relative Risk
#'   with unidimensional theta values.
#'   
#' @param X         Random sample (can be vector or matrix) which includes
#'   exposure and covariates.
#'   
#' @param thetahat  Estimative of \code{theta} for the Relative Risk function
#'   
#' @param thetalow  Lower bound of the confidence interval (vector)
#'   
#' @param thetaup   Upper band of the confidence interval (vector)
#'   
#' @param rr        Function for Relative Risk which uses parameter 
#'   \code{theta}. The order of the parameters shound be \code{rr(X, theta)}.
#'   
#'   
#'   **Optional**
#'   
#' @param weights    Survey \code{weights} for the random sample \code{X}.
#'   
#' @param method    Either \code{empirical} (default) or \code{approximate}.
#'   
#' @param Xvar      Variance of exposure levels.
#'   
#' @param confidence Confidence level \% (default: 95)
#'   
#' @param check_thetas Check that thetas are correctly specified
#'   
#' @param check_integrals Check that counterfactual and relative risk's expected
#'   values are well defined for this scenario
#'   
#' @param check_exposure  Check that exposure \code{X} is positive and numeric
#'   
#' @param check_rr        Check that Relative Risk function \code{rr} equals 
#'   \code{1} when evaluated at \code{0}
#'   
#' @param deriv.method.args \code{method.args} for
#'  \code{\link[numDeriv]{hessian}}.
#'  
#' @param deriv.method      \code{method} for \code{\link[numDeriv]{hessian}}.
#'  Don't change this unless you know what you are doing.
#'
#' @author Rodrigo Zepeda Tello \email{rodrigo.zepeda@insp.mx}
#' @author Dalia Camacho García Formentí \email{daliaf172@gmail.com}
#'   
#' @examples 
#' 
#' #Example 1: Exponential Relative Risk
#' #--------------------------------------------
#' set.seed(18427)
#' X <- rnorm(1000, 3,.7)
#' thetahat <- 0.4
#' thetalow <- 0.1
#' thetaup  <- 0.7
#' paf.confidence.one2one(X, thetahat, thetalow, thetaup, function(X, theta){exp(theta*X)})
#' 
#' #Approximate method:
#' Xmean <- mean(X)
#' Xvar  <- var(X)
#' paf.confidence.one2one(Xmean, thetahat, thetalow, thetaup, function(X, theta){exp(theta*X)}, 
#' Xvar = Xvar, method = "approximate")
#' 
#' #Example 2: Multivariate example
#' #--------------------------------------------
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
#' #Approximate method:
#' Xmean <- colMeans(X)
#' Xvar  <- var(X)
#' paf.confidence.one2one(Xmean, thetahat, thetalow, thetaup, 
#' function(X, theta){exp(theta[1]*X[,1] + theta[2]*X[,2])}, 
#' Xvar = Xvar, method = "approximate")
#' 
#' @export

paf.confidence.one2one <- function(X, thetahat, thetalow, thetaup, rr, 
                                   weights =  rep(1/nrow(as.matrix(X)),nrow(as.matrix(X))), 
                                   confidence = 95, check_thetas = TRUE,
                                   deriv.method.args = list(), deriv.method = c("Richardson", "complex"),
                                   method = c("empirical","approximate"), Xvar = var(X),
                                   check_exposure = TRUE, check_rr = TRUE, check_integrals = TRUE){
  
  #Get method from vector
  .method <- as.vector(method)[1]
  
  #Get thetavar as matrix
  thetavar <- matrix(0, ncol = length(thetahat), nrow = length(thetahat))
  
  #Calculate the PIF with confidence intervals
  .upper <- paf.confidence.inverse(X, thetahat = thetaup,  thetavar = thetavar, rr = rr,
                                   method = .method, weights = weights, confidence = confidence,
                                   nsim = 1, deriv.method.args = deriv.method.args, deriv.method = deriv.method,
                                   force.min = FALSE, check_thetas = check_thetas, Xvar = Xvar)
  .lower <- paf.confidence.inverse(X, thetahat = thetalow,  thetavar = thetavar, rr = rr,
                                   method = .method, weights = weights, confidence = confidence,
                                   nsim = 1, deriv.method.args = deriv.method.args, deriv.method = deriv.method,
                                   force.min = FALSE, check_thetas = FALSE, Xvar = Xvar)
  .point <- pif(X, thetahat, rr = rr, weights = weights, method = .method, check_exposure = check_exposure, 
                check_rr = check_rr, check_integrals = check_integrals, Xvar = Xvar, 
                deriv.method.args = deriv.method.args, deriv.method = deriv.method)
  
  
  #Return
  .confint        <- c(.lower["Lower_CI"], .point, .upper["Upper_CI"])
  names(.confint) <- c("Lower_CI","Point_Estimate","Upper_CI")
  
  #Return
  return(.confint)
}
  
  