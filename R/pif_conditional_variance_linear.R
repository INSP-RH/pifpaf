#' @title Conditional Variance for the Potential Impact Fraction
#' 
#' @description Function that calculates the conditional variance of the potential impact fraction (linearization).
#' 
#' @param X         Random sample (can be vector or matrix) which includes exposure and covariates.
#' 
#' @param thetahat  Estimative of \code{theta} for the Relative Risk function
#' 
#' @param rr        Function for relative risk
#' 
#' 
#' **Optional**
#' @param cft       Function \code{cft(X)} for counterfactual. Leave empty for 
#'                  the Population Attributable Fraction \code{PAF} where counterfactual is 0 exposure
#' 
#' @param weights    Survey \code{weights} for the random sample \code{X}
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
#' pif.conditional.variance.linear(X, thetahat,  rr = function(X, theta){exp(theta*X)})
#' 
#' #Example with linear counterfactual
#' cft      <- function(X){0.3*X}
#' pif.conditional.variance.linear(X, thetahat,  rr = function(X, theta){exp(theta*X)}, cft)
#' 
#' #Example with theta and X multivariate
#' set.seed(18427)
#' X1 <- rnorm(2000, 3,.5)
#' X2 <- rnorm(2000,3,.5)
#' X  <- as.matrix(cbind(X1,X2))
#' thetahat <- c(0.1, 0.03)
#' rr        <- function(X, theta){
#'            .X <- matrix(X, ncol = 2)
#'            exp(theta[1]*.X[,1] + theta[2]*.X[,2])
#'            }#' cft <- function(X){0.5*X}
#' pif.conditional.variance.linear(X, thetahat, rr, cft) 
#' 
#' @export

pif.conditional.variance.linear <- function(X, thetahat, rr, 
                                cft = function(Varx){matrix(0,ncol = ncol(as.matrix(Varx)), nrow = nrow(as.matrix(Varx)))}, 
                                weights =  rep(1/nrow(as.matrix(X)),nrow(as.matrix(X)))){
  
  check.cft(cft,X)
  .X         <- as.matrix(X)
  rr.fun     <- function(X){rr(X, thetahat)} 
  rr.cft.fun <- function(X){rr(cft(X), thetahat)}
  
  s        <- sum(weights)
  s2       <- sum(weights^2)
  .RO      <- weighted.mean(rr(.X,.theta), weights)
  .RC      <- weighted.mean(rr(.cft.X,.theta), weights)
  
  .varRO   <- s2*(s/(s^2 - s2)) * weighted.mean((rr(.X,.theta) - .RO)^2, weights)
  .varRC   <- s2*(s/(s^2 - s2)) * weighted.mean((rr(.cft.X,.theta) - .RC)^2, weights)
  .covRORC <- s2*s/(s^2 - s2)   * (weighted.mean((rr(.X, .theta))*(rr(.cft.X,theta)), weights)-.RO*.RC)

  Var       <- (1/.RO)^2*((.RC/.RO)^2*.varRO+.varRC-2*(.RC/.RO)*.covRORC)
  return(Var)
}
  
