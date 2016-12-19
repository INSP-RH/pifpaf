#' @title Population Attributable Fraction with exponential relative risk function
#' 
#' @description Function that calculates the population attributable fraction with exponential relative risk function
#' 
#' @param X         Random sample (can be vector or matrix) which includes exposure and covariates.
#' 
#' @param thetahat  Estimative of \code{theta} for the Relative Risk function
#' 
#' 
#' **Optional**
#' 
#' @param weights   Survey \code{weights} for the random sample \code{X}
#' 
#' @param method    Either \code{empirical} (default), \code{kernel} or \code{approximate}. 
#' 
#' @param Xvar      Variance of exposure levels.
#' 
#' @param ktype    \code{kernel} type from  \code{gaussian}, \code{epanechnikov}, \code{rectangular},
#'                  \code{triangular}, \code{biweight}, \code{cosine}, \code{optcosine}
#' 
#' @param bw        Smoothing bandwith parameter from density
#' 
#' @param adjust    Adjust bandwith parameter from density
#' 
#' @param npoints   Number of points
#' 
#' @return paf      Estimate of population attributable fraction
#' 
#' @author Rodrigo Zepeda Tello \email{rodrigo.zepeda@insp.mx}
#' @author Dalia Camacho García Formentí
#' 
#' @examples 
#' 
#' #Example with risk given by HR
#' set.seed(18427)
#' X <- rnorm(100,3,.5)
#' thetahat <- 0.12
#' paf.exponential(X, thetahat)
#' 
#' #Same example with kernel method
#' paf.exponential(X, thetahat,  method = "kernel")
#' 
#' #Same example with approximate method
#' Xmean <- mean(X)
#' Xvar  <- var(X)
#' paf.exponential(X, thetahat, method = "approximate", Xvar = Xvar)
#' 
#' #Multidimensional example
#' X     <- matrix(c(rnorm(100,2,.7),rnorm(100,4,1)),ncol=2)
#' theta <- c(0.3,0.1)
#' paf.exponential(X,theta)
#' @export

paf.exponential <- function(X, thetahat,
                weights =  rep(1/nrow(as.matrix(X)),nrow(as.matrix(X))), 
                method = c("empirical", "kernel", "approximate"),
                Xvar = var(X),
                ktype = "epanechnikov", bw = "nrd0", adjust = 1, npoints = 1000){
  .X    <- as.matrix(X)
  if(dim(.X)[2] != length(thetahat)){
    stop("The amount of parameters in theta must be equal to the number of exposure values and covariates in each observation")
  }
  rr   <- function(x,theta){
    x   <- as.matrix(x)
    n   <- dim(x)[1]
    sol <- c()
    for(i in 1:n){
      xi <- x[i,]
      sol <- c(sol, exp(theta%*%xi))
    }
    return(sol)
  }
  .paf <- pif(X = X, thetahat = thetahat, rr = rr, 
              weights = weights, method = method, 
              Xvar = Xvar,
              ktype = ktype, bw = bw, adjust = adjust, 
              npoints = npoints)
  
  return(.paf)
}


