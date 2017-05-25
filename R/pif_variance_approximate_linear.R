#'@title Approximate Variance for the Potential Impact Fraction using the 
#'  approximate method
#'  
#'@description Function that calculates approximate variance to the potential 
#'  impact fraction.
#'  
#'@param X      Mean value of exposure levels from a cross-sectional random 
#'  sample.
#'  
#'@param Xvar      Variance of exposure levels.
#'  
#'@param thetahat  Estimator (vector or matrix) of \code{theta} for the Relative
#'  Risk function.
#'  
#'@param thetavar   Estimator of variance of \code{thetahat}
#'  
#'@param rr        Function for Relative Risk which uses parameter \code{theta}.
#'  The order of the parameters shound be \code{rr(X, theta)}.
#'  
#'  
#'  **Optional**
#'  
#'@param cft       Function \code{cft(X)} for counterfactual. Leave empty for 
#'  the Population Attributable Fraction \code{\link{paf}} where counterfactual 
#'  is 0 exposure.
#'  
#'@param nsim      Number of simulations for estimation of variance
#'  
#'@param check_thetas Checks that theta parameters are correctly inputed
#'  
#'@param check_cft  Check if counterfactual function \code{cft} reduces 
#'  exposure.
#'  
#'@param check_xvar Check if it is covariance matrix.
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
#'@param is_paf Force evaluation as paf
#'  
#' @author Rodrigo Zepeda-Tello \email{rzepeda17@gmail.com}
#' @author Dalia Camacho-García-Formentí \email{daliaf172@gmail.com}
#'  
#'@seealso \code{\link{pif.variance.linear}} for \code{\link{pif}} variance and 
#'  \code{\link{pif.confidence}} for confidence intervals of \code{\link{pif}}
#'  
#' @examples 
#' \dontrun{
#' #Example 1: Exponential Relative risk
#' #--------------------------------------------
#' set.seed(46987)
#' rr       <- function(X,theta){exp(X*theta)}
#' cft      <- function(X){0.5*X}
#' X        <- rbeta(100, 2, 3)
#' Xmean    <- data.frame(mean(X))
#' Xvar     <- var(X)
#' theta    <- 1.2
#' thetavar <- 0.15
#' pif.variance.approximate.linear(Xmean, theta, rr, thetavar, Xvar, cft) 
#'
#' #Example 2: Multivariate example
#' #--------------------------------------------
#' X1       <- rnorm(1000,3,.5)
#' X2       <- rnorm(1000,4,1)
#' X        <- data.frame(cbind(X1,X2))
#' Xmean    <- matrix(colMeans(X), ncol = 2)
#' Xvar     <- cov(X)
#' theta    <- c(0.12, 0.17)
#' thetavar  <- matrix(c(0.001, 0.00001, 0.00001, 0.004), byrow = TRUE, nrow = 2)
#' rr       <- function(X, theta){exp(theta[1]*X[,1] + theta[2]*X[,2])}
#' pif.variance.approximate.linear(Xmean, theta, rr, thetavar, Xvar,
#' cft = function(X){cbind(0.5*X[,1],0.4*X[,2])}, check_integrals = FALSE)
#' }
#'@importFrom MASS mvrnorm
#'@importFrom stats weighted.mean
#'@importFrom numDeriv grad
#'@keywords internal
#'@export


pif.variance.approximate.linear <- function(X, thetahat, rr, thetavar, Xvar,
                                     cft = NA, 
                                     check_thetas = TRUE, check_cft = TRUE, check_xvar = TRUE, check_rr = TRUE, 
                                     check_integrals = TRUE, check_exposure = TRUE, deriv.method.args = list(), 
                                     deriv.method = c("Richardson", "complex"), nsim = 1000,
                                     is_paf = FALSE){
  
  #Function for checking that thetas are correctly inputed
  if(check_thetas){ check.thetas(thetavar, thetahat, NA, NA, "approximate") }
  
  #Check counterfactual function
  if(!is.function(cft)){ is_paf <- TRUE }
  
  #Set X as matrix
  if(check_xvar) {Xvar   <- check.xvar(Xvar)}
  .Xvar   <- Xvar
  .X      <- as.matrix(X)
  
  #Set a minimum for nsim
  .nsim        <- max(nsim,10)
  
  #Get the expected pif
  .pifexp <- function(theta){
    pif.approximate(X = .X, Xvar = .Xvar, thetahat = theta, rr = rr, cft = cft, 
                    deriv.method.args = deriv.method.args, deriv.method = deriv.method, 
                    check_exposure = check_exposure, check_rr = check_rr, 
                    check_integrals = check_integrals, is_paf = is_paf)
  }
  
  #Get the variance of pif
  .pifvar <- function(theta){
    
    #Rewrite functions as functions of X
    rr.fun.x <- function(X){
      rr(X,theta)
    }
    
    rr.fun.cft <- function(X){
      cftX  <- cft(X)
      rr(cftX, theta)
    }
    
    #Estimate relative risks 
    dR0 <- as.matrix(grad(rr.fun.x, .X, method = deriv.method,
                          method.args = deriv.method.args))
    R0  <- rr.fun.x(.X)
    
    #Estimate cft
    if (is_paf){
      dR1 <- 0
      R1  <- 1
    } else {
      dR1 <- as.matrix(grad(rr.fun.cft, .X, method = deriv.method,
                            method.args = deriv.method.args))  
      R1  <- rr.fun.cft(.X)
    }
    
    #Calculate Taylor this way
    aux <- ((dR1*R0 - dR0*R1)/(R0^2))
    vr  <- t(aux)%*%.Xvar%*%(aux)
    
    return(vr)
  }
  
  #Get expected value and variance of that
  .meanvec   <- rep(NA, .nsim)
  .varvec    <- rep(NA, .nsim)
  .thetasim  <- mvrnorm(.nsim, thetahat, thetavar, empirical = TRUE)
  for (i in 1:.nsim){
    .meanvec[i]  <- .pifexp(.thetasim[i,])
    .varvec[i]   <- .pifvar(.thetasim[i,])
  }
  
  #Get variance of that
  .varpif <- var(.meanvec) + mean(.varvec)
  
  #Return variance
  return(.varpif)
  
}


