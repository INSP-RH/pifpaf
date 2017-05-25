#' @title Confidence intervals for the Population Attributable Fraction, using the inverse method
#' 
#' @description Confidence intervals for the Population Attributable Fraction for relative risk inyective functions, the PAF is inyective, and intervals can be calculated for the relative risk, and then transformed to PAF CI. 
#' 
#' @param X         Random sample (\code{data.frame}) which includes exposure
#'   and covariates.
#' 
#' @param thetahat  Estimative of \code{theta} for the Relative Risk function.
#' 
#' @param thetavar   Estimator of variance of \code{thetahat}.
#' 
#' @param rr        Function for Relative Risk which uses parameter 
#'   \code{theta}. The order of the parameters shound be \code{rr(X, theta)}.
#' 
#' 
#' **Optional**
#' 
#' @param weights    Survey \code{weights} for the random sample \code{X}.
#' 
#' @param method    Either \code{empirical} (default) or \code{approximate}. 
#' 
#' @param Xvar      Variance of exposure levels.
#' 
#' @param confidence Confidence level \% (default: \code{95})
#' 
#' @param nsim      Number of simulations (default: \code{1000})
#' 
#' @param force.min Boolean indicating whether to force the \code{rr} to have a 
#'                  minimum value of 1 instead of 0 (not recommended).
#'                  
#' @param check_thetas Checks that theta parameters are correctly inputed    
#'             
#' @param deriv.method.args \code{method.args} for
#'  \code{\link[numDeriv]{hessian}}. Only if \code{"approximate"} method is chosen.
#'  
#' @param deriv.method      \code{method} for \code{\link[numDeriv]{hessian}}.
#'  Don't change this unless you know what you are doing. Only if \code{"approximate"} method is chosen.
#'  
#' @note The \code{force.min} option forces the relative risk \code{rr} to have a minimum of \code{1} and thus
#' an \code{rr < 1} is NOT possible. This is only for when absolute certainty is had that \code{rr > 1} and should
#' be used under careful consideration. The confidence interval to acheive such an \code{rr} is based on the paper
#' by Do Le Minh and Y. .s. Sherif                  
#' 
#' @author Rodrigo Zepeda-Tello \email{rzepeda17@gmail.com}
#' @author Dalia Camacho-García-Formentí \email{daliaf172@gmail.com}
#' 
#' @examples 
#'  \dontrun{
#' #Example 1: Exponential Relative Risk
#' #--------------------------------------------
#' set.seed(18427)
#' X <- rnorm(100,0.3,.05)
#' thetahat <- 0.4
#' thetavar <- 0.1
#' paf.confidence.inverse(X, thetahat, function(X, theta){exp(theta*X)}, thetavar)
#' 
#'
#' #With approximate method
#' Xmean <- mean(X)
#' Xvar  <- var(X)
#' paf.confidence.inverse(Xmean, thetahat, 
#' function(X, theta){exp(theta*X)}, thetavar, Xvar = Xvar, method = "approximate")
#' 
#' #We can force PAF's CI to be >= 0 (only if it is certain)
#' paf.confidence.inverse(X, thetahat, 
#' function(X, theta){exp(theta*X)}, thetavar, force.min = TRUE)
#' 
#' #Example 2: Multivariate Relative Risk
#' #--------------------------------------------
#' set.seed(18427)
#' X1 <- rnorm(1000,0.3,.05)
#' X2 <- rnorm(1000,0.3,.05)
#' X  <- as.matrix(cbind(X1,X2))
#' thetahat <- c(0.12, 0.03)
#' thetavar <- matrix(c(0.1, 0, 0, 0.4), byrow = TRUE, nrow = 2)
#' rr <- function(X, theta){exp(theta[1]*X[,1] + theta[2]*X[,2])}
#' paf.confidence.inverse(X, thetahat, rr, thetavar) 
#' 
#' #Same example with approximate method
#' Xmean    <- matrix(colMeans(X), ncol = 2)
#' Xvar     <- cov(X)
#'paf.confidence.inverse(Xmean, thetahat, rr=rr, thetavar = thetavar, 
#'method = "approximate", Xvar = Xvar)
#'}
#' 
#' @keywords internal
#' @export

paf.confidence.inverse <- function(X, thetahat, rr, thetavar,
                                   weights =  rep(1/nrow(as.matrix(X)),nrow(as.matrix(X))),
                                   method = c("empirical", "approximate"),
                                   nsim = 1000, confidence = 95, 
                                   deriv.method.args = list(), 
                                   deriv.method = c("Richardson", "complex"), 
                                   force.min = FALSE, 
                                   check_thetas = TRUE, Xvar = var(X)){
  
  #Get method from vector
  .method   <- as.vector(method)[1]
  
  #Change variance to matrix
  .thetavar <- as.matrix(thetavar)
  
  rr.CI <- switch (.method,
          empirical = {
            risk.ratio.confidence(X = X, thetahat = thetahat, 
                                  rr = rr, thetavar = .thetavar,
                                  weights =  weights, nsim = nsim, 
                                  confidence = confidence, check_thetas = check_thetas, 
                                  force.min = force.min)
            },
          approximate = {
            .Xvar <- check.xvar(Xvar)
            risk.ratio.approximate.confidence(X = X, Xvar = .Xvar, thetahat = thetahat, 
                                              rr = rr,  thetavar = .thetavar, nsim = nsim, 
                                              confidence = confidence,
                                              deriv.method.args = deriv.method.args,
                                              deriv.method = deriv.method,
                                              check_thetas = check_thetas, 
                                              force.min = force.min)
          },
          stop("Incorrect method. Please specify 'approximate' or 'empirical'.")
  )
  #Compute the PAF intervals
  .cipaf         <- 1-1/rr.CI
  names(.cipaf)  <- c("Lower_CI","Point_Estimate","Upper_CI") 
  
  #Return ci
  return(.cipaf)
  
}