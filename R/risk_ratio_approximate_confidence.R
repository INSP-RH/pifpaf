#'@title Approximate Confidence intervals for the Risk Ratio Integral
#'  
#'@description Function that approximates the confidence interval for the 
#'  integral \deqn{ \int RR(x;\theta)f(x)dx } where \eqn{f(x)} is the density 
#'  function of the exposure X, \eqn{ RR(x;\theta)} the relative risk of the 
#'  exposure with associated parameter \eqn{\theta}. In particular this is an 
#'  approximation when only mean and variance of exposure known
#'  
#'@param X      Mean value of exposure levels from a cross-sectional random 
#'  sample.
#'  
#'@param Xvar      Variance of exposure levels.
#'  
#'@param thetahat  Estimator (vector or matrix) of \code{theta} for the Relative
#'  Risk function.
#'  
#'@param rr        Function for Relative Risk which uses parameter \code{theta}.
#'  The order of the parameters shound be \code{rr(X, theta)}.
#'  
#'@param thetavar   Estimator of variance of \code{thetahat}
#'  
#'  **Optional**
#'  
#'@param nsim      Number of simulations
#'  
#'@param deriv.method.args \code{method.args} for 
#'  \code{\link[numDeriv]{hessian}}.
#'  
#'@param deriv.method      \code{method} for \code{\link[numDeriv]{hessian}}. 
#'  Don't change this unless you know what you are doing.
#'  
#'@param confidence Confidence level \% (default: \code{95})
#'  
#'@param check_thetas Checks that \code{theta} parameters are correctly inputed.
#'  
#'@param force.min Boolean indicating whether to force the \code{rr} to have a 
#'  minimum value of \code{1} instead of \code{0} (not recommended).
#'  
#'@note When a sample of the exposure \code{X} is available the method 
#'  \link{risk.ratio.confidence} should be prefered.
#'  
#'@note The \code{force.min} option forces the relative risk \code{rr} to have a
#'  minimum of \code{1} and thus an \code{rr < 1} is NOT possible. This is only 
#'  for when absolute certainty is had that \code{rr > 1} and should be used 
#'  under careful consideration. The confidence interval to acheive such an 
#'  \code{rr} is based on the paper by Do Le Minh and Y. .s. Sherif
#'  
#'@seealso \link{risk.ratio.confidence} for a method when there is a sample of 
#'  the exposure.
#'  
#'@references Sherif, Y. .s. (1989). The lower confidence limit for the mean of 
#'  positive random variables. Microelectronics Reliability, 29(2), 151-152.
#'  
#' @author Rodrigo Zepeda-Tello \email{rzepeda17@gmail.com}
#' @author Dalia Camacho-García-Formentí \email{daliaf172@gmail.com}
#'  
#' @examples 
#' \dontrun{
#' #' #Example 1: Exponential Relative Risk
#' #--------------------------------------------
#' set.seed(18427)
#' X        <- data.frame(rnorm(100))
#' thetahat <- 0.1
#' thetavar <- 0.2
#' Xmean    <- data.frame(mean(X[,1]))
#' Xvar     <- var(X[,1])
#' rr      <- function(X,theta){exp(X*theta)}
#' risk.ratio.approximate.confidence(Xmean, Xvar, thetahat, rr, thetavar)
#' 
#' #We can force RR'.s CI to be >= 1.
#' #This should be done with extra methodological (philosophical) care as 
#' #RR>= 1 should only be assumed with absolute mathematical certainty
#' risk.ratio.approximate.confidence(Xmean, Xvar, thetahat, rr, thetavar, force.min = TRUE)
#'
#' #Example 2: Multivariate Relative Risk
#' #--------------------------------------------
#' set.seed(18427)
#' X1        <- rnorm(1000)
#' X2        <- runif(1000)
#' X         <- data.frame(t(colMeans(cbind(X1,X2))))
#' Xvar      <- cov(cbind(X1,X2))
#' thetahat  <- c(0.02, 0.01)
#' thetavar  <- matrix(c(0.1, 0, 0, 0.4), byrow = TRUE, nrow = 2)
#' rr        <- function(X, theta){exp(theta[1]*X[,1] + theta[2]*X[,2])}
#' risk.ratio.approximate.confidence(X, Xvar, thetahat, rr, thetavar) 
#'}
#' 
#'@importFrom MASS mvrnorm
#'@importFrom numDeriv grad
#'@importFrom stats qnorm
#'@keywords internal
#'@export

risk.ratio.approximate.confidence <- function(X, Xvar, thetahat, rr, thetavar,
                                  nsim = 1000, confidence = 95, 
                                  deriv.method.args = list(), 
                                  deriv.method = c("Richardson", "complex"), 
                                  check_thetas = TRUE,
                                  force.min = FALSE){
  
  #Get confidence
  check.confidence(confidence)
  .alpha <- max(0, 1 - confidence/100)
  .Z     <- qnorm(1-.alpha/2)
  
  #Check 
  .thetavar <- as.matrix(thetavar)
  if(check_thetas){ check.thetas(.thetavar, thetahat, NA, NA, "risk.ratio") }
  
  #Get number of simulations
  .nsim  <- max(10, ceiling(nsim))
  
  #To matrix
  .X      <- as.matrix(X)
  .Xvar   <- matrix(Xvar, ncol = sqrt(length(Xvar)))
  
  #Calculate the conditional expected value as a function of theta
  .Risk  <- function(.theta){
    .paf  <- pif.approximate(X =  .X, Xvar = .Xvar, thetahat = .theta, rr = rr,
                             cft=NA, deriv.method = deriv.method,
                             deriv.method.args = deriv.method.args,
                             check_exposure = FALSE, 
                             check_rr = FALSE, check_integrals = FALSE, 
                             is_paf = TRUE)
    return( 1/(1-.paf))
  }
  
  #Calculate the conditional variance as a function of theta
  .Variance <- function(.theta){
    rr.fun.x <- function(X){
      rr(X,.theta)
    }
    
    .dR0   <- as.matrix(grad(rr.fun.x, .X, method = deriv.method[1],
                             method.args = deriv.method.args))
    .var  <- t(.dR0)%*%.Xvar%*%(.dR0)
    return(.var)
  }
  
  #Get expected value and variance of that
  .meanvec   <- rep(NA, .nsim)
  .varvec    <- rep(NA, .nsim)
  .thetasim  <- mvrnorm(.nsim, thetahat, thetavar, empirical = TRUE)
  for (i in 1:.nsim){
    .meanvec[i]  <- .Risk(.thetasim[i,])
    .varvec[i]   <- .Variance(.thetasim[i,])
  }
  
  #Get variance of that
  .inversevarpaf <- var(.meanvec) + mean(.varvec)
  
  #Create the confidence intervals
  .squareroot <- .Z*sqrt(.inversevarpaf)
  .ciup       <- .Risk(thetahat) + .squareroot
  .cilow      <- (.Risk(thetahat)^2)/.ciup #Force rr > 0
  
  #If minimum is forced to 1 correct CI
  if (force.min){
    .cilow      <- ((.Risk(thetahat) - 1)^2)/(.ciup-1) + 1
  }
  
  #Compute the Risk Ratio intervals
  .cirisk        <- c("Lower_CI"       = .cilow, 
                      "Point_Estimate" =  .Risk(thetahat) , 
                      "Upper_CI"       = .ciup )
  
  #Return variance
  return(.cirisk)
  
}
