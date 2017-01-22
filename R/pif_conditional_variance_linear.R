#' @title Conditional Variance for the Potential Impact Fraction
#' 
#' @description Function that calculates the conditional variance of the potential impact fraction (linearization).
#' 
#' @param X         Random sample (can be vector or matrix) which includes exposure and covariates.
#' 
#' @param thetahat  Estimative of \code{theta} for the Relative Risk function
#' 
#' @param rr        Function for Relative Risk which uses parameter 
#'   \code{theta}. The order of the parameters shound be \code{rr(X, theta)}.
#' 
#' 
#' **Optional**
#' @param cft       Function \code{cft(X)} for counterfactual. Leave empty for 
#'                  the Population Attributable Fraction \code{PAF} where counterfactual is 0 exposure.
#' 
#' @param weights    Survey \code{weights} for the random sample \code{X}.
#' 
#' @param check_cft  Check if counterfactual function \code{cft} reduces exposure.
#' 
#' @param is_paf     Force evaluation of paf
#' 
#' @author Rodrigo Zepeda Tello \email{rodrigo.zepeda@insp.mx}
#' @author Dalia Camacho García Formentí \email{daliaf172@gmail.com}
#' 
#' @seealso \code{pif.variance.linear} for \code{pif} variance without conditioning on \code{theta}  and \code{pif.variance.loglinear} for variance of \code{log(pif)}
#' 
#' @examples 
#' 
#' #Example 1: Exponential Relative risk
#' #--------------------------------------------
#' set.seed(18427)
#' X <- rnorm(100,3,.5)
#' thetahat <- 0.12
#' 
#' #When no counterfactual is specified calculates PAF
#' pif.conditional.variance.linear(X, thetahat,  rr = function(X, theta){exp(theta*X)})
#' 
#' #Example with linear counterfactual
#' cft      <- function(X){0.3*X}
#' pif.conditional.variance.linear(X, thetahat,  rr = function(X, theta){exp(theta*X)}, cft)
#' 
#' #Example 2: Multivariate case
#' #--------------------------------------------
#' set.seed(18427)
#' X1 <- rnorm(100, 3,.5)
#' X2 <- runif(100, 1, 1.5)
#' X  <- as.matrix(cbind(X1,X2))
#' thetahat <- c(0.1, 0.03)
#' rr        <- function(X, theta){
#'            .X <- matrix(X, ncol = 2)
#'            exp(theta[1]*.X[,1] + theta[2]*.X[,2])
#'            }
#' cft <- function(X){0.5*X}
#' pif.conditional.variance.linear(X, thetahat, rr, cft) 
#'
#' @importFrom stats weighted.mean
#' 
#' @export


  
pif.conditional.variance.linear <- function(X, thetahat, rr, 
                                            cft = function(Varx){matrix(0,ncol = ncol(as.matrix(Varx)), nrow = nrow(as.matrix(Varx)))}, 
                                            weights =  rep(1/nrow(as.matrix(X)),nrow(as.matrix(X))), check_cft = TRUE,
                                            is_paf = FALSE){
  
    #Check counterfactual  
    if (check_cft){check.cft(cft,X)}
  
    #Set as counterfactual
    .X     <- as.matrix(X)
    .cft.X <- as.matrix(cft(.X))
    
    #Sum of weights
    s        <- sum(weights)
    s2       <- sum(weights^2)
    
    #Estimate RR part
    .RO      <- weighted.mean(rr(.X, thetahat), weights)
    .varRO   <- s2*(s/(s^2 - s2)) * weighted.mean((rr(.X, thetahat) - .RO)^2, weights)
    
    #Estimate counterfactual part
    if (is_paf){
      .RC      <- 1
      .varRC   <- 0
      .covRORC <- 0
    } else {
      .RC      <- weighted.mean(rr(.cft.X, thetahat), weights)  
      .varRC   <- s2*(s/(s^2 - s2)) * weighted.mean((rr(.cft.X, thetahat) - .RC)^2, weights)
      .covRORC <- s2*s/(s^2 - s2)   * (weighted.mean((rr(.X, thetahat))*(rr(.cft.X, thetahat)), weights)-.RO*.RC)
    }
    
    
    #Calculate variance
    .Var       <- (1/.RO)^2*((.RC/.RO)^2*.varRO+.varRC-2*(.RC/.RO)*.covRORC)
    
    return(.Var)
  
}
  