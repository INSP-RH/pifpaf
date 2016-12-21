#' @title Confidence Intervals for the Population Attributable Fraction using loglinear method
#' 
#' @description Confidence intervals for the Population Attributable Fraction.
#' 
#' @param X     Mean value of exposure levels from a previous study.
#' 
#' @param Xvar      Variance of exposure levels from a previous study.
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
#' @param nsim      Number of simulations
#' 
#' @param weights   Survey \code{weights} for the random sample \code{X}
#' 
#' @param method    Either \code{empirical} (default), \code{kernel} or \code{approximate}. 

#' @param confidence Confidence level \% (default 95)
#' 
#' @param check_thetas Check that theta parameters are correctly inputed
#' 
#' @param ktype    \code{kernel} type from  \code{gaussian}, \code{epanechnikov}, \code{rectangular},
#'                  \code{triangular}, \code{biweight}, \code{cosine}, \code{optcosine}
#' 
#' @param bw        Smoothing bandwith parameter from density
#' 
#' @param cft.check Boolean. Check the counterfactual function?
#' 
#' @param adjust    Adjust bandwith parameter from density
#' 
#' @param npoints   Number of points
#' 
#' @author Rodrigo Zepeda Tello \email{rodrigo.zepeda@insp.mx}
#' @author Dalia Camacho García Formentí 
#' 
#' @examples 
#' #Example 1
#' set.seed(46987)
#' rr      <- function(X,theta){exp(X*theta)}
#' X       <- rnorm(100, 3, 0.1)
#' theta   <- 0.4
#' thetavar <- 0.001
#' paf.confidence.loglinear(X, theta, rr, thetavar, method = "empirical")
#' 
#' #Same example using approximate method
#' paf.confidence.loglinear(mean(X), theta, rr, thetavar, Xvar = var(X), method = "approximate")
#'
#'\dontrun{
#'#Example 2: Compare paf.variance.approximate with paf.variance.loglinear
#'X1       <- rnorm(100,3,.5)
#'X2       <- rnorm(100,4,1)
#'X        <- as.matrix(cbind(X1,X2))
#'Xmean    <- colMeans(X)
#'Xvar     <- cov(X)
#'theta    <- c(0.12, 0.17)
#'thetavar  <- matrix(c(0.001, 0.00001, 0.00001, 0.004), byrow = TRUE, nrow = 2)
#'rr       <- function(X, theta){exp(theta[1]*X[,1] + theta[2]*X[,2])}
#'
#'#Using empirical
#'paf.confidence.loglinear(X, theta, rr, thetavar)  
#'
#'#Using approximate
#'paf.confidence.loglinear(Xmean, theta, rr, thetavar, Xvar = Xvar, method = "approximate")
#'}
#'
#' @import MASS numDeriv
#' @export

paf.confidence.loglinear <- function(X, thetahat, rr, thetavar,
                                     weights = rep(1/nrow(as.matrix(X)), nrow(as.matrix(X))),
                                     confidence = 95, nsim = 1000, check_thetas = TRUE, 
                                     method = c("empirical", "kernel", "approximate"), Xvar = var(X),
                                     cft.check = TRUE, ktype = "epanechnikov", bw = "nrd0", adjust = 1,
                                     npoints = 1000){
  
  #Get confidence
  check.confidence(confidence)
  .alpha  <- max(0, 1 - confidence/100)
  .Z      <- qnorm(1-.alpha/2)
  
  #Get number of simulations
  .nsim   <- max(10, ceiling(nsim))
  
  #Get method
  .method <- method[1]
  
  #Check 
  .thetavar <- as.matrix(thetavar)
  if(check_thetas){ check.thetas(.thetavar, thetahat, NA, NA, "log") }
  
  #Set X as matrix
  check.exposure(X)
  .Xvar   <- check.xvar(Xvar)
  .X      <- matrix(X, ncol = ncol(.Xvar))
  
  #Calculate the conditional expected value as a function of theta
  .logpafexp <- function(.theta){
    .paf  <- pif(.X, thetahat = .theta, rr = rr, weights = weights, method = .method, Xvar = .Xvar,
                 cft.check = cft.check, ktype = ktype, bw = bw, adjust = adjust, npoints = npoints)
    .R0   <- 1/(1-.paf)
    return( -log(.R0) )
  }
  
  #Get 1 - paf for point estimate
  .paf     <- pif(.X, thetahat = thetahat, rr = rr, weights = weights, method = .method, Xvar = .Xvar,
                  cft.check = cft.check, ktype = ktype, bw = bw, adjust = adjust, npoints = npoints)
  .inverse <- 1 - .paf
  
  #Calculate the conditional variance as a function of theta
  switch(.method,
         approximate = {
           .logpafvar <- function(.theta){
             
             #Create function of RR  
             rr.fun.x <- function(X){
               rr(X, .theta)
             }
             
             #Work with linealization
             dR0   <- as.matrix(grad(rr.fun.x, .X))
             R0    <- rr(.X, .theta)
             .var  <- t(1/R0*dR0)%*%.Xvar%*%(1/R0*dR0)
             
             return(.var)
           }
         }, empirical = {
           
           #Estimate variance
           .logpafvar <- function(.theta){
             s  <- sum(weights)
             s2 <- sum(weights^2)
             .RO  <- weighted.mean(rr(.X,.theta), weights)
             .var <- (1/.RO^2)*( s / (s^2 - s2) ) * weighted.mean((rr(.X,.theta) - .RO)^2, weights)
             return(.var)
           }
         }
  )
  
  
  #Get expected value and variance of that
  .logmeanvec   <- rep(NA, .nsim)
  .logvarvec    <- rep(NA, .nsim)
  .thetasim     <- mvrnorm(.nsim, thetahat, thetavar, empirical = TRUE)
  for (i in 1:.nsim){
    .logmeanvec[i]  <- .logpafexp(.thetasim[i,])
    .logvarvec[i]   <- .logpafvar(.thetasim[i,])
  }
  
  #Get variance of that
  .logvarpaf <- var(.logmeanvec) + mean(.logvarvec)
  
  #Create the confidence intervals
  .zqrt      <- .Z*sqrt(.logvarpaf)

  
  #Compute the PAF intervals
  .cipaf     <- 1-c("Lower" = .inverse*exp(.zqrt), "Point_Estimate" =  .inverse, "Upper" = .inverse*exp(-.zqrt) )
  
  #Return variance
  return(.cipaf)
  
}