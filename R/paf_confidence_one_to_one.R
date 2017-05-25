#' @title Approximate Confidence Intervals for the Population Attributable 
#'   Fraction with one to one expected value of RR function, and unidimensional
#'   theta values
#'   
#' @description Function that calculates approximate confidence intervals of the
#'   Population Attributable Fraction \code{\link{paf}} considering a one to one
#'   Relative Risk \code{rr} with unidimensional \code{theta} parameter values
#'   
#' @param X         Random sample (\code{data.frame}) which includes exposure
#'   and covariates.
#'   
#' @param thetahat  Estimative of \code{theta} for the Relative Risk function
#'   
#' @param thetalow  Lower bound of the confidence interval.
#'   
#' @param thetaup   Upper bound of the confidence interval.
#'   
#' @param rr        Function for Relative Risk which uses parameter 
#'   \code{theta}. The order of the parameters shound be \code{rr(X, theta)}.
#'   
#'   
#'   **Optional**
#'   
#' @param weights    Survey \code{weights} for the random sample \code{X}.
#'   
#' @param method    Either \code{"empirical"} (default) or \code{"approximate"}.
#'   
#' @param Xvar      Variance of exposure levels (for \code{"approximate"}
#'   method).
#'   
#' @param confidence Confidence level \% (default: \code{95})
#' 
#' @param confidence_theta Confidence level \% of \code{theta} corresponding to
#' the interval [\code{thetalow}, \code{thetaup}] (default: \code{99}\%).
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
#'   \code{\link[numDeriv]{hessian}}.
#'   
#' @param deriv.method      \code{method} for \code{\link[numDeriv]{hessian}}. 
#'   Don't change this unless you know what you are doing.
#'   
#' @author Rodrigo Zepeda-Tello \email{rzepeda17@gmail.com}
#' @author Dalia Camacho-García-Formentí \email{daliaf172@gmail.com}
#'
#'   
#' @examples 
#' 
#' #Example 1: Exponential Relative Risk
#' #--------------------------------------------
#' set.seed(18427)
#' X <- data.frame(rnorm(1000, 3,.7))
#' thetahat <- 0.4
#' thetalow <- 0.1
#' thetaup  <- 0.7
#' paf.confidence.one2one(X, thetahat, function(X, theta){exp(theta*X)}, 
#' thetalow, thetaup)
#' 
#' #Approximate method:
#' Xmean <- data.frame(mean(X[,1]))
#' Xvar  <- var(X[,1])
#' paf.confidence.one2one(Xmean, thetahat, function(X, theta){exp(theta*X)}, 
#' thetalow, thetaup, Xvar = Xvar, method = "approximate")
#' 
#' #Example 2: Multivariate example
#' #--------------------------------------------
#' set.seed(18427)
#' X1 <- rnorm(1000,3,.7)
#' X2 <- rnorm(1000,3,.7)
#' X  <- data.frame(X1,X2)
#' thetahat <- c(0.12, 0.03)
#' thetalow <- c(0.05, 0.01)
#' thetaup  <- c(0.25, 0.06)
#' rr <- function(X, theta){exp(theta[1]*X[,1] + theta[2]*X[,2])}
#' paf.confidence.one2one(X, thetahat, rr, thetalow, thetaup) 
#' 
#' #Approximate method:
#' Xmean <- data.frame(t(colMeans(X)))
#' Xvar  <- var(X)
#' paf.confidence.one2one(Xmean, thetahat, 
#' function(X, theta){exp(theta[1]*X[,1] + theta[2]*X[,2])}, 
#' thetalow, thetaup, 
#' Xvar = Xvar, method = "approximate")
#' 
#' @keywords internal
#' @export

paf.confidence.one2one <- function(X, thetahat, rr, thetalow, thetaup, 
                                   weights =  rep(1/nrow(as.matrix(X)),nrow(as.matrix(X))), 
                                   confidence = 95, 
                                   confidence_theta = 99,
                                   check_thetas = TRUE,
                                   deriv.method.args = list(), 
                                   deriv.method = c("Richardson", "complex"),
                                   method = c("empirical","approximate"), Xvar = var(X),
                                   check_exposure = TRUE, check_rr = TRUE, 
                                   check_integrals = TRUE){
  # Check theta values are correctly defined
  
  if(check_thetas){
  check.thetas(thetavar = NA, thetahat = thetahat, thetalow = thetalow, thetaup = thetaup,
               method = "one2one")
  }
  #Get method from vector
  .method <- as.vector(method)[1]
  
  #Get thetavar as matrix
  .thetavar <- matrix(0, ncol = length(thetahat), nrow = length(thetahat))
  
  #Get new confidence assuming worst-case scenario
  .confidence <- confidence/confidence_theta
    
  #Calculate the PIF with confidence intervals
  .upper <- paf.confidence.inverse(X, thetahat = thetaup, rr = rr, thetavar = .thetavar, 
                                    weights = weights, method = .method,  nsim = 1,
                                   confidence = .confidence,
                                   deriv.method.args = deriv.method.args, deriv.method = deriv.method,
                                   force.min = FALSE, check_thetas = check_thetas, Xvar = Xvar)
  .lower <- paf.confidence.inverse(X, thetahat = thetalow, rr = rr, thetavar = .thetavar, 
                                   weights = weights,  method = .method, nsim = 1,
                                   confidence = .confidence,
                                    deriv.method.args = deriv.method.args, deriv.method = deriv.method,
                                   force.min = FALSE, check_thetas = FALSE, Xvar = Xvar)
  .point <- paf(X, thetahat, rr = rr,  method = .method,  weights = weights,
                Xvar = Xvar,
                deriv.method.args = deriv.method.args, deriv.method = deriv.method,
                check_exposure = check_exposure, check_rr = check_rr, 
                check_integrals = check_integrals)
  
  
  #Rename
  .confint        <- c(min(.lower["Lower_CI"], .upper["Upper_CI"]), .point, 
                       max(.lower["Lower_CI"], .upper["Upper_CI"]))
  names(.confint) <- c("Lower_CI","Point_Estimate","Upper_CI")
  
  
  #Return
  return(.confint)
}
  
  