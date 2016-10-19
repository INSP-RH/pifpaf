#' @title Point Estimate of the Potential Impact Fraction via the Empirical Method
#' 
#' @description Function that calculates the potential impact fraction via the empirical method
#' 
#' @param X         Random sample (can be vector or matrix) which includes exposure and covariates.
#' 
#' @param thetahat  Estimative of \code{theta} for the Relative Risk function
#' 
#' @param rr        Function for relative risk
#' 
#' 
#' **Optional**
#' 
#' @param cft       Function \code{cft(X)} for counterfactual. Leave empty for \code{PAF} where counterfactual is 0 exposure
#' 
#' @param weights   Survey \code{weights} for the random sample \code{X}
#' 
#' @author Rodrigo Zepeda Tello \email{rodrigo.zepeda@insp.mx}
#' @author Dalia Camacho García Formentí
#' 
#' @examples 
#' 
#' #Example with risk given by HR
#' set.seed(18427)
#' X <- rnorm(100)
#' thetahat <- 0.12
#' pif.empirical(X, thetahat, function(X, theta){exp(theta*X)})
#' 
#' #Same example considering counterfactual of halfing exposure
#' pif.empirical(X, thetahat, function(X, theta){exp(theta*X)}, cft = function(X){ 0.5*X })
#' 
#' #Example with linear relative risk
#' pif.empirical(X, thetahat, function(X, theta){theta*X + 1}, cft = function(X){ 0.5*X })
#' 
#' @import matrixStats
#' 
#' @export


pif.empirical <- function(X, thetahat, rr, 
                          cft = function(Varx){matrix(0,ncol = ncol(as.matrix(Varx)), nrow = nrow(as.matrix(Varx)))}, 
                          weights =  rep(1/nrow(as.matrix(X)),nrow(as.matrix(X)))){
  
  #Set X as matrix
  .X  <- as.matrix(X)
  
  #Check that rr is 1 when X = 0
  check.rr(.X, thetahat, rr)
  
  #Estimate weighted sums
  .mux   <- sum(rr(.X,thetahat) * weights)
  .mucft <- sum(rr(cft(.X),thetahat) * weights)
  
  #Check that integrals make sense
  check.integrals(.mux, .mucft)
  
  #Calculate PIF
  .pif   <- 1 - .mucft/.mux

  return(.pif)
  
}


