#' @title Approximate Variance for the Potential Impact Fraction
#' 
#' @description Function that calculates approximate variance of the potential impact fraction (linearization).
#' 
#' @param X         Random sample (can be vector or matrix) which includes exposure and covariates.
#' 
#' @param thetahat  Maximum likelihood estimator of \code{theta} for the Relative Risk function
#'
#' @param thetavar  Estimator of variance of thetahat
#'
#' @param rr        Function for Relative Risk which uses parameter 
#'   \code{theta}. The order of the parameters shound be \code{rr(X, theta)}.
#' 
#' 
#' **Optional**
#' 
#' @param cft       Function \code{cft(X)} for counterfactual. Leave empty for 
#'   the Population Attributable Fraction \code{\link{paf}} where counterfactual
#'   is 0 exposure.
#'   
#' @param weights   Normalized survey \code{weights} for the sample \code{X}.
#' 
#' @param nsim      Number of simulations for estimation of variance
#' 
#' @param check_thetas Checks that theta parameters are correctly inputed
#' 
#' @param check_cft  Check if counterfactual function \code{cft} reduces exposure.
#' 
#' @param check_exposure Check if exposure > 0
#' 
#' @param is_paf Boolean to force paf evaluation.
#' 
#' @importFrom MASS mvrnorm 
#' @importFrom stats weighted.mean
#' 
#' @author Rodrigo Zepeda Tello \email{rodrigo.zepeda@insp.mx}
#' @author Dalia Camacho García Formentí \email{daliaf172@gmail.com}
#' 
#' @seealso \code{pif.variance.linear} for \code{pif} variance and \code{pif.variance.loglinear} for variance of \code{log(pif)}
#' and \code{pif.confidence} for confidence intervals of \code{pif}
#'
#' @examples 
#' 
#' #Example 1: Exponential Relative risk
#' #--------------------------------------------
#' set.seed(18427)
#' X <- rnorm(100,3,.5)
#' thetahat <- 0.12
#' thetavar <- 0.1
#' pif.variance.linear(X, thetahat, thetavar, function(X, theta){exp(theta*X)})
#' 
#' #Same example with linear counterfactual
#' cft      <- function(X){0.3*X}
#' pif.variance.linear(X, thetahat, thetavar, function(X, theta){exp(theta*X)}, cft)
#' 
#' #Example 2: Multivariate case
#' #--------------------------------------------
#' set.seed(18427)
#' X1 <- rnorm(100, 3,.5)
#' X2 <- runif(100, 1, 1.5)
#' X  <- as.matrix(cbind(X1,X2))
#' thetahat <- c(0.1, 0.03)
#' thetavar <- matrix(c(0.1, 0, 0, 0.05), byrow = TRUE, nrow = 2)
#' rr        <- function(X, theta){
#'            .X <- matrix(X, ncol = 2)
#'            exp(theta[1]*.X[,1] + theta[2]*.X[,2])
#'            }
#' cft <- function(X){0.5*X}
#' pif.variance.linear(X, thetahat, thetavar, rr, cft) 
#' 
#' @export

pif.variance.linear <- function(X, thetahat, thetavar, rr, 
                                cft = function(Varx){matrix(0,ncol = ncol(as.matrix(Varx)), nrow = nrow(as.matrix(Varx)))}, 
                                weights =  rep(1/nrow(as.matrix(X)),nrow(as.matrix(X))), 
                                check_thetas = TRUE,  check_cft = TRUE, 
                                check_exposure = TRUE, nsim = 1000, is_paf = FALSE){
  #Set X as matrix
  .X    <- as.matrix(X)
  
  #Check exposure values are greater than zero
  if(check_exposure){ check.exposure(.X) }
  
  #Set a minimum for nsim
  .nsim <- max(nsim, 10)
  
  #Check 
  .thetavar <- as.matrix(thetavar)
  if(check_thetas){ check.thetas(.thetavar, thetahat, NA, NA, "linear") }
  
  #Get the conditional expected pif
  .pifexp <- function(theta){
    pif(X = X, thetahat = theta, rr = rr, cft = cft, weights = weights,
        method = "empirical",  check_exposure = FALSE, 
        check_rr = FALSE, check_integrals = FALSE, is_paf = is_paf)
  }
  
  #Get the conditional variance of pif
  .pifvar <- function(theta){
    vr <- pif.conditional.variance.linear(X = .X, thetahat = theta, rr = rr, cft = cft, weights = weights, 
                                          check_cft = check_cft, is_paf = is_paf)
    return(vr)
  }
  
  #Get expected value and variance of that
  .meanvec   <- rep(NA, .nsim)
  .varvec    <- rep(NA, .nsim)
  .thetasim  <- mvrnorm(.nsim, thetahat, .thetavar, empirical = TRUE)
  for (i in 1:.nsim){
    .meanvec[i]  <- .pifexp(.thetasim[i,])
    .varvec[i]   <- .pifvar(.thetasim[i,])
  }
  
  #Get variance of that
  meanvec <<- .meanvec
  varvec <<- .varvec
  .varpif <- var(.meanvec) + mean(.varvec)
  
  #Return variance
  return(.varpif)
  
}