#' @title Approximate Variance for the Potential Impact Fraction
#' 
#' @description Function that calculates approximate variance of the potential impact fraction (linearization).
#' 
#' @param X         Random sample (can be vector or matrix) which includes exposure and covariates.
#' 
#' @param thetahat  Estimative of \code{theta} for the Relative Risk function
#' 
#' @param thetasd   Estimator of standard error of thetahat (usually standard error) 
#' 
#' @param rr        Function for relative risk
#' 
#' 
#' **Optional**
#'@param cft       Function \code{cft(X)} for counterfactual. Leave empty for 
#'                  the Population Attributable Fraction \code{PAF} where counterfactual is 0 exposure
#' 
#' @param nsim      Number of simulations for estimation of variance
#' 
#' @param weights    Survey \code{weights} for the random sample \code{X}
#' 
#' @param check_thetas Checks that theta parameters are correctly inputed
#' 
#' @import MASS stats
#' 
#' @author Rodrigo Zepeda Tello \email{rodrigo.zepeda@insp.mx}
#' @author Dalia Camacho García Formentí 
#' 
#' @examples 
#' 
#' #Example with risk given by HR (PAF)
#' set.seed(18427)
#' X <- rnorm(100,3,.5)
#' thetahat <- 0.12
#' thetavar <- 0.1
#' pif.variance.linear(X, thetahat, thetavar, function(X, theta){exp(theta*X)})
#' paf.variance.linear(X, thetahat, thetavar, function(X, theta){exp(theta*X)})
#' 
#' #Example with linear counterfactual
#' cft      <- function(X){0.3*X}
#' pif.variance.linear(X, thetahat, thetavar, function(X, theta){exp(theta*X)}, cft)
#' 
#' #Example with theta and X multivariate
#' set.seed(18427)
#' X1 <- rnorm(2000, 3,.5)
#' X2 <- rnorm(2000,3,.5)
#' X  <- as.matrix(cbind(X1,X2))
#' thetahat <- c(0.1, 0.03)
#' thetasd <- matrix(c(0.1, 0, 0, 0.05), byrow = TRUE, nrow = 2)
#' rr        <- function(X, theta){
#'            .X <- matrix(X, ncol = 2)
#'            exp(theta[1]*.X[,1] + theta[2]*.X[,2])
#'            }#' cft <- function(X){0.5*X}
#' pif.variance.linear(X, thetahat, thetasd, rr, cft) 
#' 
#' 
#' @export

pif.variance.linear <- function(X, thetahat, thetasd, rr, 
                                cft = function(Varx){matrix(0,ncol = ncol(as.matrix(Varx)), nrow = nrow(as.matrix(Varx)))}, 
                                weights =  rep(1/nrow(as.matrix(X)),nrow(as.matrix(X))), check_thetas = FALSE,  nsim = 100){
  #Set X as matrix
  .X  <- as.matrix(X)
  
  #Set a minimum for nsim
  .nsim        <- max(nsim,10)
  
  #Check 
  .thetasd <- as.matrix(thetasd)
  if(check_thetas){ check.thetas(.thetasd, thetahat, NA, NA, "linear") }
  
  #Get the expected pif
  .pifexp <- function(theta){
    R0 <- weighted.mean(rr(.X, theta), weights)
    RC <- weighted.mean(rr(cft(X), theta), weights)
    return(1-RC/R0)
  }
  
  #Get the variance of pif
  .pifvar <- function(theta){
    vr <- pif.conditional.variance.linear(X, theta, rr, cft, weights)
    return(vr)
  }
  
  #Get expected value and variance of that
  .meanvec   <- rep(NA, .nsim)
  .varvec    <- rep(NA, .nsim)
  .thetasim  <- mvrnorm(.nsim, thetahat, thetasd)
  for (i in 1:.nsim){
    .meanvec[i]  <- .pifexp(.thetasim[i,])
    .varvec[i]   <- .pifvar(.thetasim[i,])
  }
  
  #Get variance of that
  .varpif <- var(.meanvec) + mean(.varvec)
  
  #Return variance
  return(.varpif)
  
}