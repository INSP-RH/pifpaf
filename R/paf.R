#' @title Population Attributable Fraction
#' 
#' @description Function that calculates the population attributable fraction
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
#' paf(X, thetahat, function(X, theta){exp(theta*X)})
#' 
#' #Same example with kernel method
#' paf(X, thetahat, function(X, theta){exp(theta*X)}, method = "kernel")
#' 
#' #Example with linear relative risk
#' paf(X, thetahat, function(X, theta){theta*X + 1})
#' 
#' #Multivariate example
#' set.seed(18427)
#' X1 <- rnorm(1000,3,.5)
#' X2 <- rnorm(1000,3,.5)
#' X  <- as.matrix(cbind(X1,X2))
#' thetahat <- c(0.12, 0.03)
#' rr <- function(X, theta){exp(theta[1]*X[,1] + theta[2]*X[,2])}
#' paf(X, thetahat, rr) 
#' 
#' \dontrun{
#' #Multivariate cases cannot be evaluated with kernel method
#' paf(X, thetahat, rr, method = "kernel") 
#' }
#' 
#' 
#' @export

paf <- function(X, thetahat, rr, 
                weights =  rep(1/nrow(as.matrix(X)),nrow(as.matrix(X))), 
                method = c("empirical", "kernel", "approximate"),
                Xvar = var(X),
                ktype = "epanechnikov", bw = "nrd0", adjust = 1, npoints = 1000){
  
  
  .paf <- pif(X = X, thetahat = thetahat, rr = rr, 
              weights = weights, method = method, 
              Xvar = Xvar,
              ktype = ktype, bw = bw, adjust = adjust, 
              npoints = npoints)
  
  return(.paf)
}


