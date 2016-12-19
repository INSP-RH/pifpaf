#' @title Confidence intervals for the Population Attributable Fraction, using the inverse method
#' 
#' @description Confidence intervals for the Population Attributable Fraction for relative risk inyective functions, the PAF is inyective, and intervals can be calculated for the relative risk, and then transformed to PAF CI. This function works for both the empirical method and the approximate method.
#' 
#' @param X         Random sample (can be vector or matrix) which includes exposure and covariates. Or mean exposure from a previous study if no sample is available.
#' 
#' @param thetahat  Estimative of \code{theta} for the Relative Risk function
#' 
#' @param thetavar   Estimator of variance of thetahat 
#' 
#' @param rr        Function for relative risk
#' 
#' 
#' **Optional**
#' 
#' @param weights   Survey \code{weights} for the random sample \code{X}
#' 
#' @param method    Either \code{empirical} (default) or \code{approximate}. 
#' 
#' @param Xvar      Variance of exposure levels.
#' 
#' @param nsim      Number of simulations
#' 
#' @param confidence Confidence level \% (default 95)
#' 
#' @param force.min Boolean indicating whether to force the PAF to have a 
#'                  minimum value of 0 instead of allowing negative values (not recommended).
#'                  
#' @param check_thetas Checks that theta parameters are correctly inputed                
#' 
#' @author Rodrigo Zepeda Tello \email{rodrigo.zepeda@insp.mx}
#' @author Dalia Camacho García Formentí 
#' 
#' @examples 
#' 
#' #Example with risk given by HR
#' set.seed(18427)
#' X <- rnorm(100,3,.5)
#' thetahat <- 0.4
#' thetavar <- 0.1
#' paf.confidence.inverse(X, thetahat, thetavar, function(X, theta){exp(theta*X)})
#' 
#' #With larger sample the variance reduces
#' set.seed(18427)
#' X <- rnorm(10000,3,.5)
#' thetahat <- 0.12
#' thetavar <- 0.1
#' paf.confidence.inverse(X, thetahat, thetavar, function(X, theta){exp(theta*X)})
#' 
#' #We can force PAF's CI to be >= 0 
#' paf.confidence.inverse(X, thetahat, thetavar, function(X, theta){exp(theta*X)}, force.min = TRUE)
#' 
#' #Example with theta and X multivariate
#' set.seed(18427)
#' X1 <- rnorm(1000,3,.5)
#' X2 <- rnorm(1000,3,.5)
#' X  <- as.matrix(cbind(X1,X2))
#' thetahat <- c(0.12, 0.03)
#' thetavar <- matrix(c(0.1, 0, 0, 0.4), byrow = TRUE, nrow = 2)
#' rr <- function(X, theta){exp(theta[1]*X[,1] + theta[2]*X[,2])}
#' paf.confidence.inverse(X, thetahat, thetavar, rr) 
#' 
#' #Example: Approximate method
#'  set.seed(46987)
#'  rr      <- function(X,theta){exp(X*theta)}
#'  X       <- rnorm(100,3.2,1)
#'  Xmean   <- mean(X)
#'  Xvar    <- var(X)
#'  theta   <- 0.4
#'  thetavar <- 0.001
#'  paf.confidence.inverse(Xmean, thetahat = theta, thetavar = thetavar,
#'   rr=rr, method = "approximate", Xvar = Xvar)
#'  paf.confidence.inverse(X, thetahat = theta, thetavar = thetavar, 
#'  rr=rr, method = "empirical", Xvar = Xvar)
#'
#'#Example: Multidimensional example using approximate method
#'X1       <- rnorm(1000,3,.5)
#'X2       <- rnorm(1000,4,1)
#'X        <- as.matrix(cbind(X1,X2))
#'Xmean    <- colMeans(X)
#'Xvar     <- cov(X)
#'theta    <- c(0.12, 0.17)
#'thetavar  <- matrix(c(0.001, 0.00001, 0.00001, 0.004), byrow = TRUE, nrow = 2)
#'rr       <- function(X, theta){exp(theta[1]*X[,1] + theta[2]*X[,2])}
#'paf.confidence.inverse(Xmean, thetahat = theta, thetavar = thetavar, 
#'rr=rr, method = "approximate", Xvar = Xvar)
#' 
#' @import MASS
#' @export

paf.confidence.inverse <- function(X, thetahat, thetavar, rr, weights =  rep(1/nrow(as.matrix(X)),nrow(as.matrix(X))),
                                   nsim = 1000, confidence = 95, force.min = FALSE, check_thetas = TRUE, method = c("empirical", "approximate"), Xvar = NA){
  
  #Get method from vector
  .method <- as.vector(method)[1]
  
  #Change variance to matrix
  .thetavar <- as.matrix(thetavar)
  
  #Function checking confidence is correctly specified
  check.confidence(confidence)
  
  #Function for checking that thetas are correctly inputed
  if(check_thetas){ check.thetas(.thetavar, thetahat, NA, NA, "inverse") }
  
 rr.CI<- switch (.method,
          empirical = {
            risk.ratio.confidence(X = X, thetahat = thetahat, 
                                              thetavar = .thetavar, rr = rr, 
                                              weights =  weights, nsim = nsim, 
                                              confidence = confidence, force.min = force.min)
            },
          approximate = {
            Xvar <- check.xvar(Xvar)
            risk.ratio.approximate.confidence(Xmean = X, Xvar = Xvar, thetahat = thetahat, 
                                                            thetavar = .thetavar, rr = rr, nsim = nsim, 
                                                            confidence = confidence, force.min = force.min)
          },
          risk.ratio.confidence(X = X, thetahat = thetahat, 
                                thetavar = .thetavar, rr = rr, 
                                weights =  weights, nsim = nsim, 
                                confidence = confidence, force.min = force.min)
  )
  #Compute the PAF intervals
  .cipaf         <- 1-1/rr.CI
    
    #Return variance
    return(.cipaf)
  
}