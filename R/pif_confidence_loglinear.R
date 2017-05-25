#' @title Confidence intervals for the Potencial Impact Fraction, using the
#'   loglinear method
#'   
#' @description Confidence intervals for the Potencial Impact Fraction for
#'   relative risk inyective functions, the pif is inyective, and intervals can
#'   be calculated for log(pif), and then transformed to pif CI.
#'   
#' @param X         Random sample (can be vector or matrix) which includes
#'   exposure and covariates.
#'   
#' @param thetahat  Estimative of \code{theta} for the Relative Risk function
#'   
#' @param thetavar   Estimator of variance of thetahat
#'   
#' @param rr        Function for Relative Risk which uses parameter 
#'   \code{theta}. The order of the parameters shound be \code{rr(X, theta)}.
#'   
#'   
#'   **Optional**
#' @param cft       Function \code{cft(X)} for counterfactual. Leave empty for 
#'   the Population Attributable Fraction \code{PAF} where counterfactual is 0
#'   exposure.
#'   
#' @param weights    Survey \code{weights} for the random sample \code{X}.
#'   
#' @param nsim      Number of simulations for estimation of variance.
#'   
#' @param confidence Confidence level \% (default 95)
#'   
#' @param check_thetas Check that theta parameters are correctly inputed
#'   
#' @param check_exposure  Check that exposure \code{X} is positive and numeric
#'   
#' @param check_cft  Check if counterfactual function \code{cft} reduces
#'   exposure.
#'   
#' @param is_paf Force evaluation of paf
#'   
#' @author Rodrigo Zepeda-Tello \email{rzepeda17@gmail.com}
#' @author Dalia Camacho-García-Formentí \email{daliaf172@gmail.com}
#'  
#'   
#' @examples 
#' 
#' #Example 1: Exponential Relative risk
#' #--------------------------------------------
#' set.seed(18427)
#' X        <- rnorm(100,5,1)
#' thetahat <- 0.4
#' thetavar <- 0.1
#' cft      <- function(X){sqrt(X)}
#' pif.confidence.loglinear(X, thetahat, thetavar, rr = function(X, theta){exp(theta*X)}, cft)
#' 
#' #Same example with linear counterfactual
#' a    <- 0.5
#' cft  <- function(X){a*X}
#' pif.confidence.loglinear(X, thetahat, thetavar, function(X, theta){exp(theta*X)})
#' 
#' #' #Example 2: Multivariate Relative Risk
#' #--------------------------------------------
#' set.seed(18427)
#' X1        <- rnorm(100, 4,0.01)
#' X2        <- runif(100,0.4,2)
#' X         <- as.matrix(cbind(X1,X2))
#' thetahat  <- c(0.12, 0.03)
#' thetavar  <- matrix(c(0.01, 0, 0, 0.04), byrow = TRUE, nrow = 2)
#' rr        <- function(X, theta){
#'            .X <- as.matrix(X, ncol = 2)
#'            exp(theta[1]*.X[,1] + theta[2]*.X[,2])
#'            }
#' pif.confidence.loglinear(X, thetahat, thetavar, rr) 
#' 
#' @import MASS
#' @keywords internal
#' @export

pif.confidence.loglinear <- function(X, thetahat, thetavar, rr, 
                                     cft = NA,
                                     weights =  rep(1/nrow(as.matrix(X)),nrow(as.matrix(X))),
                                     nsim = 100, confidence = 95, check_thetas = TRUE, check_exposure = TRUE,
                                     check_cft = TRUE, is_paf = FALSE){
  
  #To matrix
  .X     <- as.data.frame(X)
  n      <- nrow(.X)
  
  #Get confidence
  check.confidence(confidence)
  .alpha <- max(0, 1 - confidence/100)
  .Z     <- qnorm(1-.alpha/2)

  #Get number of simulations
  .nsim  <- max(10, ceiling(nsim))
  
  #Check 
  .thetavar <- as.matrix(thetavar)
  if(check_thetas){ check.thetas(.thetavar, thetahat, NA, NA, "log") }
  
  #Check exposure levels
  if(check_exposure){check.exposure(.X)}
  
  #Check counterfactual exists
  if (!is.function(cft)){ is_paf <- TRUE}
  
  #Check expected values are finite
  infinite <- FALSE
  
  #Calculate the conditional expected value as a function of theta
  .logpifexp <- function(.theta){
    .pif <- pif(X = X, thetahat = .theta, rr = rr, cft = cft, method = "empirical", 
                weights = weights, Xvar=NA,
                check_exposure = FALSE, check_rr = FALSE,
                check_integrals = FALSE, is_paf = is_paf)
    return(1 - .pif)
  }
  
  #Get inverse for multiplying
  .inverse   <- 1 - pif(X = X, thetahat = thetahat, rr = rr, cft = cft, method = "empirical",
                        weights = weights, Xvar=NA,
                         check_exposure = FALSE, check_rr = FALSE,
                        check_integrals = FALSE, is_paf = is_paf)
  
  
  #Calculate the conditional variance as a function of theta
  .logpifvar <- function(.theta){
    s        <- sum(weights)
    s2       <- sum(weights^2)
    
    #Calculate rr
    .RO      <- weighted.mean(as.matrix(rr(.X,.theta)), weights)
    if (is.infinite(.RO)){
      warning("Expected value of Relative Risk is not finite")
      infinite <- TRUE
    } else {
      .varRO   <- (1/.RO^2)*s2*( s / (s^2 - s2) ) * 
        weighted.mean(as.matrix((rr(.X,.theta) - .RO)^2), weights)
    }
    
    
    if (is_paf){
      .RC      <- 1
      .varRC   <- 0
      .covRORC <- 0
    } else {
      
      .cft.X <- as.data.frame(cft(X))
      
      #Check counterfactual
      if(check_cft){check.cft(cft, X)}
      
      .RC      <- weighted.mean(as.matrix(rr(.cft.X,.theta)), weights)
      if (is.infinite(.RC)){
        warning("Expected value of Relative Risk under counterfactual is not finite")
        infinite <- TRUE
      } else {
        .varRC   <- (1/.RC^2)*s2*( s / (s^2 - s2) ) * 
          weighted.mean(as.matrix((rr(.cft.X,.theta) - .RC)^2), weights)
        .covRORC <- (1/(.RO*.RC))*s2*s/(s^2 - s2)   * 
          (weighted.mean(as.matrix((rr(.X, .theta))*(rr(.cft.X, .theta))), weights)-.RO*.RC)
      }  
    }
      
    if (!infinite){
      .var     <-  .varRO + .varRC - 2*.covRORC  
    } else {
      .var     <- Inf
    }
    
    return(.var)
  }
  
  #Get expected value and variance of that
  .logmeanvec   <- rep(NA, .nsim)
  .logvarvec    <- rep(NA, .nsim)
  .thetasim     <- mvrnorm(.nsim, thetahat, thetavar, empirical = TRUE)
  for (i in 1:.nsim){
    .logmeanvec[i]  <- .logpifexp(.thetasim[i,])
    .logvarvec[i]   <- .logpifvar(.thetasim[i,])
  }
  
  #Get variance of that
  .logvarpif <- var(.logmeanvec) + mean(.logvarvec)
  
  #Create the confidence intervals
  .zqrt      <- .Z*sqrt(.logvarpif)
  
  #Compute the pif intervals
  .cipif         <- 1 - c("Lower_CI" = .inverse*exp(.zqrt), 
                          "Point_Estimate" =  .inverse, 
                          "Upper_CI" = .inverse*exp(-.zqrt), 
                          "Variance_Estimate_log(PIF)" = .logvarpif)
  
  #Correct if .zqrt is infinite
  if (is.infinite(.zqrt)){.cipif["Lower_CI"] <- -Inf}
  
  #Return variance
  return(.cipif)
  
}