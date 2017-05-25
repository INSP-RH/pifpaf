#'@title Confidence Intervals for the Potential Impact Fraction when only mean
#'  and variance of exposure values is available using loglinear method
#'  
#'@description Confidence intervals for the Population Attributable Fraction for
#'  the approximate method where only mean and variance from a previous study is
#'  available.For relative risk inyective functions, the pif is inyective, and
#'  intervals can be calculated for log(pif), and then transformed to pif CI.
#'  
#'@param Xmean  Mean value of exposure levels from a cross-sectional.
#'  
#'@param Xvar   Variance of the exposure levels.
#'  
#'@param thetahat  Estimator (vector or matrix) of \code{theta} for the Relative
#'  Risk function \code{rr}
#'  
#'@param thetavar   Estimator of variance of \code{thetahat}
#'  
#'@param rr        Function for Relative Risk which uses parameter \code{theta}.
#'  The order of the parameters shound be \code{rr(X, theta)}.
#'  
#'  
#' \strong{**Optional**}
#'@param cft       Differentiable function \code{cft(X)} for counterfactual.
#'  Leave empty for the Population Attributable Fraction \code{\link{paf}} where
#'  counterfactual is 0 exposure.
#'  
#'@param confidence Concidence level (0 to 100) default = \code{95} \%
#'  
#'@param nsim      Number of simulations for estimation of variance
#'  
#'@param check_thetas Checks that theta parameters are correctly inputed
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
#'@param check_exposure  Check that exposure \code{X} is positive and numeric.
#'  
#'@param check_rr        Check that Relative Risk function \code{rr} equals 
#'  \code{1} when evaluated at \code{0}.
#'  
#'@param is_paf Boolean forcing evaluation of \code{paf}
#'  
#' @author Rodrigo Zepeda-Tello \email{rzepeda17@gmail.com}
#' @author Dalia Camacho-García-Formentí \email{daliaf172@gmail.com}
#'  
#' @examples 
#' 
#' #Example 1: Exponential Relative Risk
#' #--------------------------------------------
#' set.seed(46987)
#' rr      <- function(X,theta){exp(X*theta)}
#' cft     <- function(X){0.4*X}
#' Xmean   <- data.frame(3)
#' Xvar    <- 1
#' theta   <- 0.4
#' thetavar <- 0.001
#' pif.confidence.approximate.loglinear(Xmean, Xvar, theta, thetavar, rr, cft,
#' nsim = 1000)
#'
#' #Example 2: Multivariate Relative Risk
#' #--------------------------------------------
#'X1       <- rnorm(100,3,.5)
#'X2       <- rnorm(100,4,1)
#'X        <- data.frame(cbind(X1,X2))
#'Xmean    <- t(as.matrix(colMeans(X)))
#'Xvar     <- cov(X)
#'thetahat <- c(0.12, 0.17)
#'thetavar  <- matrix(c(0.001, 0.00001, 0.00001, 0.004), byrow = TRUE, nrow = 2)
#'rr       <- function(X, theta){exp(theta[1]*X[,1] + theta[2]*X[,2])}
#'pif.confidence.approximate.loglinear(Xmean, Xvar, thetahat, thetavar, 
#'rr, cft = function(X){0.8*X}, nsim = 100)
#' 
#'@importFrom MASS mvrnorm
#'@importFrom numDeriv hessian grad
#'@keywords internal
#'@export

pif.confidence.approximate.loglinear <- function(Xmean, Xvar, thetahat, thetavar, rr,
                                                 cft = NA, 
                                                 deriv.method.args = list(), 
                                                 deriv.method = c("Richardson", "complex"),
                                                 check_exposure = TRUE, check_rr = TRUE, check_integrals = TRUE,
                                                 nsim = 1000, confidence = 95, check_thetas = TRUE,
                                                 is_paf = FALSE){
  
  #Get confidence
  check.confidence(confidence)
  .alpha <- max(0, 1 - confidence/100)
  .Z     <- qnorm(1-.alpha/2)
  
  #Get method
  .method <- deriv.method[1]
  
  #Get number of simulations
  .nsim  <- max(10, ceiling(nsim))
  
  #Check 
  .thetavar <- as.matrix(thetavar)
  if(check_thetas){ check.thetas(.thetavar, thetahat, NA, NA, "log") }
  
  #Set X as matrix
  .Xmean  <- as.matrix(Xmean)
  .Xvar   <- as.matrix(Xvar)
  
  if(!is.function(cft)){ is_paf <- TRUE }
  
  #Calculate the conditional expected value as a function of theta
  .logpifexp <- function(.theta){
    
    .rr_cft_fun <- function(X){
      Xcft.value  <- cft(X)
      rr(Xcft.value,.theta)
    } 
    
    .rr_fun_x <- function(X){
      rr(X,.theta)
    }
    
    #Calculate the RR
    .hrr   <- hessian(.rr_fun_x, .Xmean, method = .method, method.args = deriv.method.args)
    .R0    <- .rr_fun_x(.Xmean)   + 0.5*sum(.hrr*.Xvar)
    
    #Estimate counterfactual
    if (is_paf){
      .R1    <- 1
    } else {
      .hcft  <- hessian(.rr_cft_fun, .Xmean, method = .method, method.args = deriv.method.args)
      .R1    <- .rr_cft_fun(.Xmean) + 0.5*sum(.hcft*.Xvar)
    }
    
    return( log(.R1)-log(.R0) )
  }
  
  #Inverse
  .pif      <- pif.approximate(X = .Xmean, Xvar = .Xvar, thetahat = thetahat, rr = rr, cft = cft,
                                deriv.method.args = deriv.method.args, deriv.method = deriv.method,
                               check_exposure = check_exposure, check_rr = check_rr, 
                               check_integrals = check_integrals, is_paf = is_paf)
  .inverse  <- 1 - .pif
  
  #Calculate the conditional variance as a function of theta
  
  .logpifvar <- function(.theta){
    .rr_fun_x <- function(X){
      rr(X,.theta)
    }
    
    .rr_cft_fun <- function(X){
      Xcft.value  <- cft(X)
      rr(Xcft.value,.theta)
    } 
    
    #Estimate RR part
    dR0    <- as.matrix(grad(.rr_fun_x, .Xmean,  method = .method, method.args = deriv.method.args))
    R0     <- rr(.Xmean, .theta)
    .var0  <- t(1/R0*dR0)%*%.Xvar%*%(1/R0*dR0)
    
    #Estimate counterfactual part
    if (is_paf){
      .var1  <- 0
    } else {
      dR1    <- as.matrix(grad(.rr_cft_fun, .Xmean, method = .method, method.args = deriv.method.args))
      R1     <- rr(cft(.Xmean), .theta)
      .var1  <- t(1/R1*dR1)%*%.Xvar%*%(1/R1*dR1)
    }
    
    .var   <- .var0 +.var1 
    
    return(.var)
  }
  
  #Get expected value and variance of that
  .logmeanvec   <- rep(NA, .nsim)
  .logvarvec    <- rep(NA, .nsim)
  .thetasim     <- mvrnorm(.nsim, thetahat, .thetavar, empirical = TRUE)
  for (i in 1:.nsim){
    .logmeanvec[i]  <- .logpifexp(.thetasim[i,])
    .logvarvec[i]   <- .logpifvar(.thetasim[i,])
  }
  
  #Get variance of that
  .logvarpif <- var(.logmeanvec) + mean(.logvarvec)
  
  #Create the confidence intervals
  .zqrt       <- .Z*sqrt(.logvarpif)
  
  
  #Compute the pif intervals
  .cipif         <- 1-c("Lower_CI" = .inverse*exp(.zqrt), "Point_Estimate" =  .inverse, "Upper_CI" = .inverse*exp(-.zqrt), "Estimated_Variance_log(PIF)" = .logvarpif)
  
  #Return variance
  return(.cipif)
  
}